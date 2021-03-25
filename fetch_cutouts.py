# system
import os, sys, traceback, signal
# utilities
import re, click, urllib3
import yaml as yml
from datetime import datetime
# threading
import threading, queue
# astropy
from astropy.io import fits
import astropy.units as u
# configuration & processing
from cli_config import CLIConfig
from core.survey_abc import processing_status as ProcStatus, SurveyABC

LOG_FILE = "OutLOG.txt"

#Global pool manager
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
http = urllib3.PoolManager(
    num_pools = 60,
    maxsize   = 15,
    timeout   = 30.0,
    retries   = 3,
    block     = True
)

def set_sig_handler(threads):
    def sig_handler(sig, frame):
        #signal.signal(signal.SIGINT, original_sigint)
        msg = (lambda s: f"\n{len(s)*'*'}\n{s}\n{len(s)*'*'}\n")("***  CTRL-C RECEIVED! KILLING THREADS... ***")
        print(" ".join(re.sub('\n','\n  ',msg)))
        for t in threads:
            t.die()
        sys.exit(0)
    signal.signal(signal.SIGINT,sig_handler)
    # return

# kills a thread when given into the queue
class PoisonPill:
    def __init__(self):
        pass

class WorkerThread(threading.Thread):
    def __init__(self, worker, input_q, output_q=None, *args, **kwargs):
        self.input_q = input_q
        self.output_q = output_q
        self.worker = worker
        self.kill_recieved = False
        super().__init__(*args, **kwargs)

    def run(self):
        while not self.kill_recieved:
            task = self.input_q.get()
            if type(task) is PoisonPill:
                self.input_q.task_done()
                break
            try:
                ret = self.worker(task)
                # self.input_q.task_done()

                if self.output_q:
                    self.output_q.put(item=ret)
            #catch running exception and kill
            except Exception as e:
                msg = task['survey'].sprint(f"{str(e)} killing thread\n\n")
                try:
                    date_mark = str(datetime.now()) + ": "
                    with open(LOG_FILE, "a") as logfile:
                        logfile.write("\n\n"+date_mark+msg)
                        logfile.close()
                except:
                    print("Unable to write error to log: " + msg)
                # print("ret", ret)
                print(msg)
                self.die()
                # self.input_q.task_done()
            self.input_q.task_done()

        if self.kill_recieved:
            try:
                task['survey'].print("Bye!")
            except:
                print("Bye!")
            self.input_q.task_done()

    def die(self):
        self.kill_recieved = True

# grab a FITS hdu from some survey
def get_cutout(target):
    # all fits is list of one or more dicts
    all_fits = target['survey'].set_pid(target['pid']).get_cutout(target['position'], target['size'], target['group_by'])
    return all_fits

def save_cutout(all_fits):
    originals_end="_ORIGINALS"
    try:
        saved_fits = SurveyABC.save_and_serialize(all_fits, originals_path_end=originals_end)
        msg = f"{all_fits[0]['survey']}({all_fits[0]['filter']}): [Position:{all_fits[0]['position']} at "\
            f"radius {all_fits[0]['radius']} arcmin]: All Output Files Successfully Saved!"
        date_mark = str(datetime.now()) + ": "
        with open(LOG_FILE, "a") as logfile:
            logfile.write("\n\n"+date_mark+msg)
            logfile.close()
        print(msg)
    except Exception as e:
        print("Unable to save ")
        date_mark = str(datetime.now()) + ": "
        with open(LOG_FILE, "a") as logfile:
            logfile.write("\n\n"+date_mark+str(e))
            logfile.close()



def read_in_config(yml_file):
    try:
        file_obj = open(yml_file,'r')
        file_data = yml.load(file_obj, Loader=yml.SafeLoader)
        file_obj.close()
        params = {}
        params['surveys'] = file_data['cutouts']['surveys']
        params['radius'] = file_data['cutouts']['radius']
        params['group_by'] = file_data['cutouts']['group_by']
        params['output'] = file_data['configuration']['output']
        params['overwrite'] = file_data['configuration']['overwrite']
        params['flush'] = file_data['configuration']['flush']
    except Exception as e:
        print("YAML file read error: " +str(e))
        return None
    return params

def check_batch_csv(batch_files_string):
    batch_files = batch_files_string.split(',')
    good_files = [csv for csv in batch_files if csv[-4:]=='.csv']
    bad_files = [x for x in batch_files if x not in good_files]
    if len(bad_files)>0:
        print("only CSV batch file formats accepted! skipped "+ str(bad_files))
    return good_files

#cfg is a SURVEYABC object already configured
def process_requests(cfg):
    start = datetime.now()
    grabbers = 15
    # set up i/o queues
    in_q  = queue.Queue()
    out_q = queue.Queue()

    # toss all the targets into the queue, including for all surveys
    # i.e., some position in both NVSS and VLASS and SDSS, etc.
    for task in cfg.get_procssing_stack():
        task['survey'].attach_http_pool_manager(http)
        in_q.put(task)

    # need this for ctrl-c shutdown
    threads = list()

    # spin up a bunch of worker threads to process all the data
    # in principle these could be chained furprint(str(e), "killing thread")
    # targets -> hdus -> save to file -> process to jpg -> save to file
    for _ in range(grabbers):
        thread = WorkerThread(get_cutout, in_q, out_q)
        in_q.put(PoisonPill())
        threads.append(thread)
        thread.start()

    # save all from out queue as added there
    thread = WorkerThread(save_cutout, out_q)#.start()
    threads.append(thread)
    set_sig_handler(threads) # install ctrl-c handler
    thread.start()
    in_q.join()

    out_q.put(PoisonPill()) # add killmessage to end of queue
    out_q.join()

    print("time took: " +str(datetime.now()-start))

def check_group_by_string(group_by):
    case_match = group_by.upper()
    valid = ["MOSAIC", "NONE", "DATE-OBS"]
    if "DATE" in case_match:
        case_match = "DATE-OBS"
    if case_match in valid:
        return case_match
    else:
        raise Exception(f"group_by argument {group_by} is invalid!\n Valid group by options are: MOSAIC, NONE, DATE-OBS")

def parse_surveys_string(surveys):
    # regex to match if filters in brackets next to survey name
    # e.g. WISE[w1],SDSS(g,r,i)
    if surveys:
        return re.split(",(?![^[]*\\])(?![^(]*\\))",surveys)
    else:
        return []

##HANDLE COMMAND LINE INPUT
@click.group()
def cli():
    """\b
       Command Line Cutout Fetching program.
    """
    pass

@cli.command()
@click.option('--coords','-c', 'coords', required=False, default=None)
@click.option('--name','-n', 'name', required=False)
@click.option('--radius','-r', 'radius', required=False, type=int)
@click.option('--surveys','-s', 'surveys', required=False, type=str)
@click.option('--output','-o', 'data_out', required=False)
@click.option('--groupby','-g', 'group_by', required=False)
@click.option('--config', '-cf','config_file', required=False)
@click.option('--overwrite',is_flag=True, help='overwrite existing target files (default False)')
@click.option('--flush', is_flag=True, help='flush existing target files (supersedes --overwrite)')
def fetch(overwrite, flush, coords, name, radius=None, surveys=None, data_out=None, group_by='', config_file=''):
    """
    \b
    Single cutout fetching command.
    -c 'coords' for Source coordinates OR
    -n 'name' for Source name
        example accepted coordinate formats:
        > RA,DEC or 'RA, DEC' in degrees
        > '00h42m30s', '+41d12m00s' or 00h42m30s,+41d12m00s
        > '00 42 30 +41 12 00'
        > '00:42.5 +41:12'
        if name:
        > The name of the object to get coordinates for, e.g. 'M42'
    \b
    -r 'radius' is the Integer search radius around the specified source location in arcmin.
        The cutouts will be of maximum width and height of 2*radius
    \n
    \b
    -s 'surveys' is one or several surveys comma separated without spaces between.
        Implemented surveys include:
         - VLASS
         - GLEAM
            frequencies: f1 (072-103 MHz), f2 (103-034 MHz), f3 (139-170 MHz), f4 (170-231 MHz default)
         - FIRST
         - NVSS
         - WISE
            wavelengths: W1 (3.4μm default),  W2 (4.6μm),  W3 (12μm),  W4 (22μm)
         - PANSTARRS
            filters: g, r, i (default), z, y
         - SDSS-I/II
            filters: g (default), r, i
     \b
        Filters/Frequencies/Wavelengths for each survey may be specified in the following formats:
         > "WISE(w2),SDSS[g,r]"
         > "WISE[w1],VLASS"
         > "GLEAM(f1,f3)"
         > WISE,VLASS
        If no filters are specified then the default filter is used for each.
        If surveys argument is not specified then will fetch from ALL implemented
        surveys with default filters for each survey.
    \n

    \b
    -o 'output' is the directory location to save output FITS images to.
        Output will be furthered separated into subfolders for the corresponding survey.
        Default location is a folder named 'data_out/' in this current directory.
    \n
    \b
    -g 'groupby' is an option to separate FITS results by "MOSAIC", "DATE-OBS", or "NONE" (default).
        > "MOSAIC": if the requested position and radius straddle boundaries in multiple
                    FITS images for a given survey a mosaicked FITS file will be generated
                    from all of these input images with each input image as an extension of
                    the corresponding mosaicked FITS. Mosaics are largely provided for visual
                    use only.
        > "DATE-OBS": For surveys VLASS, FIRST, NVSS, or PanSTARRS a Mosaicked FITS is made
                    (when needed) for every unique DATE-OBS.
        > "NONE" (default): All resulting FITS images in the requested survey are returned
                    without doing any mosaicking
    \n
    \b
    -cf 'config' is to specify a YAML config file for settings, ex."config.yml".
        *Note: Specified command line args will overwrite these settings.

    """
    target = coords
    is_name = False
    if not name and not coords:
        print("Invalid options: coords or name required!\n")
        return
    if name and coords:
        print('Invalid options: must enter ONE of coords or name \n')
        return
    if name:
        target = name
        is_name = True

    if not radius and not config_file:
        print("\n must specify search radius or to use config.yml (--config_file)")
        return

    if radius:
        size = radius*2

    if group_by:
        group_by = group_by.upper()

    if config_file:
        config_dict = read_in_config(config_file)
        if not config_dict:
            print(f"no valid params parsed from config file {config_file}")
            return # config file error abort
        # use YAML configs only if args not given
        if not radius:
            size = config_dict['radius']*2
        if surveys is None:
            surveys = config_dict['surveys']
        if data_out is None:
            data_out = config_dict['output']
        if not overwrite:
            overwrite = config_dict['overwrite']
        if not flush:
            flush = config_dict['flush']
        if not group_by:
            group_by = config_dict['group_by']

    if data_out is None:
        data_out = 'data_out'

    relative_path = os.path.dirname(os.path.abspath(__file__))+'/'
    out_path = os.path.join(relative_path,data_out)

    if isinstance(surveys, str):
        surveys = parse_surveys_string(surveys)

    if isinstance(group_by, str):
        try:
            group_by = check_group_by_string(group_by)
        except Exception as e:
            print(str(e))
            return
    print(f"Using args: \n image size {size} \n surveys {surveys} \n group by: {group_by}\n")

    # configuration
    cfg = CLIConfig(surveys, out_path, group_by)
    cfg.set_single_target_params(target, size, is_name)
    if flush:
        cfg.flush_old_survey_data()
    print(f"Overwrite Mode: {cfg.set_overwrite(overwrite)}")
    # MAIN CALL
    process_requests(cfg)

# Notes: http://click.palletsprojects.com/en/5.x/options/
#
@cli.command()
@click.option('--file','-f', 'batch_files_string', required=True, help='batch file(s) name(s)')
@click.option('--radius','-r', 'radius', required=False, type=int)
@click.option('--surveys','-s', 'surveys', required=False, type=str)
@click.option('--output','-o', 'data_out', required=False)
@click.option('--groupby','-g', 'group_by', required=False)
@click.option('--config', '-cf','config_file', required=False)
@click.option('--overwrite', 'overwrite', is_flag=True, help='overwrite existing duplicate target files (default True)')
@click.option('--flush', 'flush', is_flag=True, help='flush existing target files (supersedes --overwrite)')
def fetch_batch( overwrite, flush, batch_files_string, radius=None, surveys=None, data_out=None, group_by='', config_file=''):
    """
       Batch cutout fetching command.

       \b
       Fetch batch cutouts for either a single csv, or multiple csv file(s)
       of source coordinates or names.
       \b
       -f "file" The CSV file(s) must at least have separate columns named "RA" and "Dec"
          (or any of the variants below, but there can only be one variant of
          RA and one of Dec per file). A column labelled "Name" or "NAME" may also be used.
          For a given source, coordinates will be evaluated via "RA" and "Dec" if
          they are non-empty. If a line does not have a valid coordinate position,
          but does have a "Name" column value, the service will attempt to resolve
          the source name.
       \b
          Accepted variants of RA and Dec Column header names are:
          R.A.
          Right Ascension
          RA (J2000)
          R.A. (J2000)
          Right Ascension (J2000)
          RAJ2000
          DEC
          DEC.
          Declination
          DEC (J2000)
          DEC. (J2000)
          Declination (J2000)
          DecJ2000
        \b
          Source names will be resolved via the Sesame Name Resolver:
          http://vizier.u-strasbg.fr/viz-bin/Sesame
       \b
       -r 'radius' is the Integer search radius around the specified source location in arcmin.
          The cutouts will be of maximum width and height of 2*radius
       \n
       \b
       -s 'surveys' is one or several surveys comma separated without spaces between.
           - VLASS
           - GLEAM
              frequencies: f1 (072-103 MHz), f2 (103-034 MHz), f3 (139-170 MHz), f4 (170-231 MHz default)
           - FIRST
           - NVSS
           - WISE
              wavelengths: W1 (3.4μm default),  W2 (4.6μm),  W3 (12μm),  W4 (22μm)
           - PANSTARRS
              filters: g, r, i (default), z, y
           - SDSS-I/II
              filters: g (default), r, i
       \b
          Filters/Frequencies/Wavelengths for each survey may be specified in the following formats:
           > "WISE(w2),SDSS[g,r]"
           > "WISE[w1],VLASS"
           > "GLEAM(f1,f3)"
           > WISE,VLASS
       \b
          If no filters are specified then the default filter is used for each.
          If surveys argument is not specified then will fetch from ALL implemented
          surveys with default filters for each survey.
       \n

       \b
       -o 'output' is the directory location to save output FITS images to.
          Output will be furthered separated into subfolders for the corresponding survey.
          Default location is a folder named 'data_out/' in this current directory.
       \n
       \b
       -g 'groupby' is an option to group FITS results by "MOSAIC", "DATE-OBS", or "NONE" (default).
           > "MOSAIC": if the requested position and radius straddle boundaries in multiple
                       FITS images for a given survey a mosaicked FITS file will be generated
                       from all of these input images with each input image as an extension of
                       the corresponding mosaicked FITS. Mosaics are largely provided for visual
                       use only.
           > "DATE-OBS": For surveys VLASS, FIRST, NVSS, or PanSTARRS a Mosaicked FITS is made
                       (when needed) for every unique DATE-OBS.
           > "NONE" (default): All resulting FITS images in the requested survey are returned
                       without doing any mosaicking

      \n
      \b
      -cf 'config' is to specify a YAML config file for settings, ex."config.yml".
          *Note: Specified command line args will overwrite these settings.
    """
    if not radius and not config_file:
        print("\n must specify search radius or to use config.yml (--config_file)")
        return
    if radius:
        size = radius*2

    if group_by:
        group_by = group_by.upper()

    if config_file:
        config_dict = read_in_config(config_file)
        if not config_dict:
            return # config file error abort
        # use YAML configs only if args not given
        if not radius:
            size = config_dict['radius']*2
        if surveys is None:
            surveys = config_dict['surveys']
        if data_out is None:
            data_out = config_dict['output']
        if not overwrite:
            overwrite = config_dict['overwrite']
        if not flush:
            flush = config_dict['flush']
        if not group_by:
            group_by = config_dict['group_by']

    if isinstance(surveys, str):
        surveys = parse_surveys_string(surveys)

    if isinstance(group_by, str):
        try:
            group_by = check_group_by_string(group_by)
        except Exception as e:
            print(str(e))
            return
    print(f"Using args: \n image size {size} \n surveys {surveys} \n group by: {group_by}\n")

    accepted_batch_files = check_batch_csv(batch_files_string)
    if not accepted_batch_files:
        print("no valid CSV batch files specified!")
        return
    print(f"Using batch csv: {accepted_batch_files}")

    if data_out is None:
        # make output oflder resemble batch file name
        if len(accepted_batch_files) == 1:
            data_out = accepted_batch_files[0].split("/")[-1].replace('.csv', '_out')
        else:
            data_out = 'data_out'
    relative_path = os.path.dirname(os.path.abspath(__file__))+'/'
    out_path = os.path.join(relative_path,data_out)
    # configuration
    cfg = CLIConfig(surveys, out_path, group_by)
    cfg.set_batch_targets(accepted_batch_files, relative_path, size)
    if flush:
        cfg.flush_old_survey_data()
    print(f"Overwrite Mode: {cfg.set_overwrite(overwrite)}")
    process_requests(cfg)

if __name__ == "__main__":
    cli()
    print("hmm")
