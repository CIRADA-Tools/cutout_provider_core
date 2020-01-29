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
from core.survey_abc import processing_status as ProcStatus

#Global pool manager
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
http = urllib3.PoolManager(
    num_pools = 60,
    maxsize   = 15,
    timeout   = 30.0,
    retries   = 10,
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

# TODO (Issue #11): Need to vet this for thread-saftey -- tricky!
class WorkerThread(threading.Thread):
    def __init__(self, worker, input_q, output_q=None, group_by=None, *args, **kwargs):
        self.input_q = input_q
        self.output_q = output_q
        self.group_by = group_by
        self.worker = worker
        self.kill_recieved = False
        super().__init__(*args, **kwargs)

    def run(self):
        while not self.kill_recieved:
            task = self.input_q.get()
            # if it's swallowed
            if type(task) is PoisonPill:
                self.input_q.task_done()
                break
            ret = self.worker(task)
            self.input_q.task_done()
            if self.output_q:
                self.output_q.put(item=ret)

        if self.kill_recieved:
            try:
                task['survey'].print("Bye!")
            except:
                print("Bye!")

    def die(self):
        self.kill_recieved = True

# grab a FITS hdu from some survey
def get_cutout(target):
    fetched = target['survey'].set_pid(target['pid']).get_cutout(target['position'], target['size'])
    target['hdu'] = fetched['cutout']
    target['log'] = {
        'msg': fetched['message'],
        'sts': fetched['status']
    }
    target['originals'] = fetched['raw_tiles']
    return target

# save an HDU into a file
def save_cutout(target):
    def pack_message(msg,msg_log):
        p="\nLOG: "
        return re.sub(r"^\n","",((f"{p}"+f"{p}".join(msg_log.splitlines())) if msg_log != "" else "")+f"{p}{msg}")

    if target['hdu']:
        target['hdu'].writeto("{0}".format(target['filename']), overwrite=True)
        # add original raw tiles as extensions to mosaic
        if len(target['originals'])>1:
            fits_img = fits.open(target['filename'], mode='append')
            for raw in target['originals']:
                fits_img.append(raw[0])
            fits_img.close(output_verify="silentfix")
        msg = target['survey'].sprint(f"{target['filename']} SAVED!",buffer=False)
    else:
        msg = target['survey'].sprint(f"Cutout at (RA, Dec) of ({target['position'].ra.to(u.deg).value}, {target['position'].dec.to(u.deg).value}) degrees /w size={target['size']} returned None.",buffer=False)
        log_msg = pack_message(msg,target['log']['msg'])
        if not ProcStatus.touch_file(f"{target['filename']}",target['log']['sts'],log_msg):
            target['survey'].sprint("Failed to dump %s with log info...\n%s" % (
                re.sub(r"\.fits$",f".{target['log']['sts'].name}",target['filename']),
                log_msg
            ), buffer=False)
    print(msg)

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
    grabbers = 60
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
    # in principle these could be chained further, such that you could go
    # targets -> hdus -> save to file -> process to jpg -> save to file
    for _ in range(grabbers):
        thread = WorkerThread(get_cutout, in_q, out_q, cfg.group_by).start()
        in_q.put(PoisonPill())
        threads.append(thread)

    # save all from out queue as added there
    thread = WorkerThread(save_cutout, out_q).start()
    threads.append(thread)
    set_sig_handler(threads) # install ctrl-c handler
    in_q.join()

    out_q.put(PoisonPill()) # add killmessage to end of queue
    out_q.join()

    print("time took: " +str(datetime.now()-start))

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

@cli.command()
@click.option('--coords','-c', 'coords', required=False, default=None)
@click.option('--name','-n', 'name', required=False)
@click.option('--radius','-r', 'radius', required=False, type=int)
@click.option('--surveys','-s', 'surveys', required=False, type=str)
@click.option('--output','-o', 'data_out', required=False)
@click.option('--groupby','-g', 'group_by', required=False, type=str)
@click.option('--config', '-cf','config_file', required=False, help='Specify YAML config file for settings, ex. "config.yml". *Note: Specified command line args will overwrite these settings')
@click.option('--overwrite',is_flag=True, help='overwrite existing target files (default False)')
@click.option('--flush', is_flag=True, help='flush existing target files (supersedes --overwrite)')
def fetch(overwrite, flush, coords, name, radius=None, surveys=None, data_out=None, group_by=None config_file=''):
    """
    \b
    Single cutout fetching command.
    Coordinates -c 'coords' OR Source name, -n 'name'.
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
    Implemented surveys include: FIRST,VLASS,WISE,SDSS,PANSTARRS,NVSS
    \b
    Filters for each survey may be specified in the following formats:
        > "WISE(w2),SDSS[g,r]"
        > "WISE[w1],VLASS"
        > WISE,VLASS
    \b
    If no filters are specified then the default filter is used for each.
    If surveys argument is not specified then will fetch from ALL implemented
    surveys with default filters for each survey.
    \n

    \b
    -o 'output' is the directory location to save output FITS images to.
    Default location is a folder named 'data_out/' in this current directory.
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
            print("YAML surveys", surveys)
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
    print(f"Using args: \n image size {size} \n surveys {surveys}")

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
@click.option('--config', '-cf','config_file', required=False, help='Specify YAML config file for settings, ex. "config.yml". *Note: Specified command line args will overwrite these settings')
@click.option('--overwrite', 'overwrite', is_flag=True, help='overwrite existing target files (default False)')
@click.option('--flush', 'flush', is_flag=True, help='flush existing target files (supersedes --overwrite)')
def fetch_batch( overwrite, flush, batch_files_string, radius=None, surveys=None, data_out=None, config_file=''):
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
       Implemented surveys include: FIRST,VLASS,WISE,SDSS,PANSTARRS,NVSS
       \b
       Filters for each survey may be specified in the following formats:
           > "WISE(w2),SDSS[g,r]"
           > "WISE[w1],VLASS"
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
    """
    if not radius and not config_file:
        print("\n must specify search radius or to use config.yml (--config_file)")
        return
    if radius:
        size = radius*2

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

    if data_out is None:
        data_out = 'data_out'

    relative_path = os.path.dirname(os.path.abspath(__file__))+'/'
    out_path = os.path.join(relative_path,data_out)

    if isinstance(surveys, str):
        surveys = parse_surveys_string(surveys)
    print(f"Using args: \n image size {size} \n surveys {surveys}")

    accepted_batch_files = check_batch_csv(batch_files_string)
    if not accepted_batch_files:
        print("no valid CSV batch files specified!")
        return
    print(f"Using batch csv: {accepted_batch_files}")
    # configuration
    cfg = CLIConfig(surveys, out_path)
    cfg.set_batch_targets(accepted_batch_files, relative_path, size)
    if flush:
        cfg.flush_old_survey_data()
    print(f"Overwrite Mode: {cfg.set_overwrite(overwrite)}")
    process_requests(cfg)

if __name__ == "__main__":
    cli()
