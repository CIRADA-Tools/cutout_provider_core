# system
import os, sys, traceback, signal
# utilities
import re, click, urllib3
import yaml as yml
# threading
import threading, queue
# astropy
from astropy.io import fits
import astropy.units as u
# configuration & processing
from surveys.survey_config import SurveyConfig
from surveys.survey_abc import processing_status as ProcStatus

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
    def __init__(self, worker, input_q, output_q=None, *args, **kwargs):
        self.input_q = input_q
        self.output_q = output_q
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

    #print("SAVING CUTOUT ")
    if target['hdu']:
        #TODO add mosaicked in filename if mosaicked
        target['hdu'].writeto("{0}".format(target['filename']), overwrite=True)
        # add original raw tiles as extensions to mosaic
        if len(target['originals'])>1:
            fits_img = fits.open(target['filename'], mode='append')
            for raw in target['originals']:
                fits_img.append(raw[0])
            #fits_img.writeto("{0}".format(target['filename']))
            fits_img.close(output_verify="silentfix")

        msg = target['survey'].sprint(f"{target['filename']} done!",buffer=False)
    else:
        msg = target['survey'].sprint(f"Cutout at (RA, Dec) of ({target['position'].ra.to(u.deg).value}, {target['position'].dec.to(u.deg).value}) degrees /w size={target['size']} returned None.",buffer=False)
        log_msg = pack_message(msg,target['log']['msg'])
        if not ProcStatus.touch_file(f"{target['filename']}",target['log']['sts'],log_msg):
            target['survey'].print("Failed to dump %s with log info...\n%s" % (
                re.sub(r"\.fits$",f".{target['log']['sts'].name}",target['filename']),
                log_msg
            ), buffer=False)
    print(msg)

#cfg is a SURVEYABC object already configured
def process_requests(cfg):
    grabbers = 60
    savers = 1
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
        thread = WorkerThread(get_cutout, in_q, out_q)
        thread.start()
        in_q.put(PoisonPill())
        threads.append(thread)

    for _ in range(savers):
        thread = WorkerThread(save_cutout, out_q)
        thread.start()
        threads.append(thread)
    set_sig_handler(threads) # install ctrl-c handler
    in_q.join()

    # testing out 1 save to file threads (absolutely not necessary)
    for _ in range(savers):
        out_q.put(PoisonPill())
    out_q.join()

# # TODO: May not need this...
# # define the default config file with absolute path
# this_source_file_dir = re.sub(r"(.*/).*$",r"\1",os.path.realpath(__file__))
# default_config = this_source_file_dir + 'config.yml'

##HANDLE COMMAND LINE INPUT
@click.group()
def cli():
    """\b
       Survey cutout fetching script.
       Command help: <command> --help
    """
    ## TODO:  put args instructions here for the --help command
    pass

@cli.command()
@click.option('--coords','-c', 'coords', required=False, default=None)
@click.option('--name','-n', 'name', required=False)
@click.option('--radius','-r', 'radius', required=True, type=int)
@click.option('--surveys','-s', 'surveys', required=False, default=[])
@click.option('--overwrite',is_flag=True,default=None,help='overwrite existing target files')
@click.option('--flush',is_flag=True,default=None,help='flush existing target files (supersedes --overwrite)')
def fetch(coords, name, radius, surveys, overwrite, flush):
    """
    \b
    Fetch cutouts for either a single set of coordinates, 'coords' or Source name, 'name'.
    'radius' is the image cutout radius in arcmin.
    'surveys' is one or several surveys comma separated without spaces between.
    If surveys argument is not specified then will fetch from all implemented surveys.
    \n
    \b
    example accepted coordinate formats:
        > RA,DEC or 'RA, DEC' in degrees
        > '00h42m30s', '+41d12m00s' or 00h42m30s,+41d12m00s
        > '00 42 30 +41 12 00'
        > '00:42.5 +41:12'
    if name:
        > The name of the object to get coordinates for, e.g. 'M42'
    """
    target = coords
    size = radius*2
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
    if surveys:
        surveys=surveys.split(',')

    print("target", target)
    print("image size", size)
    print("surveys", surveys)
    # provide list of surveys rather than name of a config file to parse items from
    cfg = SurveyConfig(surveys)
    cfg.set_single_target_params(target, size, is_name)
    print(f"Overwrite Mode: {cfg.set_overwrite(overwrite)}")
    print(f'Flush Mode: {cfg.set_flush(flush)}')
    process_requests(cfg)

@cli.command()
@click.option('--flush',is_flag=True,default=None,help='flush existing target files')
@click.argument('config-file')
def maintenance(config_file,flush):
    """
       Cutout pipeline maintenance command.

          CONFIG_FILE = yaml configuration file

       \b
       For CONFIG_FILE formating help details see:
          python fetch_cutouts.py batch-process --help
    """
    if flush:
        # load the configuration
        click.echo(f"Using Configuration: {config_file}")
        try:
            cfg = SurveyConfig(config_file)
            cfg.force_flush()
        except Exception as e:
            click.echo("ERROR: %s\n> %s" % (e, "\n> ".join(traceback.format_exc().splitlines())))
    else:
        click.echo("Please enter option flag/s.")

# Notes: http://click.palletsprojects.com/en/5.x/options/
@cli.command()
@click.argument('config-file')
@click.option('--overwrite',is_flag=True,default=None,help='overwrite existing target files')
@click.option('--flush',is_flag=True,default=None,help='flush existing target files (supersedes --overwrite)')
def batch_process(config_file,overwrite,flush):
    """
       Batch cutout fetching command.

           CONFIG_FILE = yaml configuration file

       The following is an example of a configuration file.\n
       \b
          cutouts:
              ra_dec_deg_csv_files:
                  - targets.csv
              box_size_arcmin: 3
              surveys:
                  - VLASS
                  - WISE:
                        filters: [w1]
                  - PanSTARRS:
                        filters: [g,i]
          configuration:
              local_root: testing/data/out
              overwrite: False
              flush: True

       where targets.csv is a csvfile of source coordinates or names.\n
       \b
       The CSV file must at least have separate columns named "RA" and "Dec"
       (or any of the variants below, but there can only be one variant of
       RA and one of Dec per file). A column labelled "Name" may also be used.
       For a given source, coordinates will be evaluated via "RA" and "Dec" if
       they are non-empty. If a line does not have a valid position, but does
       have a "Name" value, the service will attempt to resolve the "Name".
       \b
        Accepted variants of RA and Dec are:
        \b
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
    """
    # load the configuration
    print(f"Using Configuration: {config_file}")
    cfg = SurveyConfig(config_file)
    print(f"Overwrite Mode: {cfg.set_overwrite(overwrite)}")
    print(f'Flush Mode: {cfg.set_flush(flush)}')
    process_requests(cfg)

if __name__ == "__main__":
    cli()
