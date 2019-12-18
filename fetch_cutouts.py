# system
import os, sys, traceback, signal
# thread-salf version of urllib
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
http = urllib3.PoolManager(
    num_pools = 60,
    maxsize   = 15,
    timeout   = 30.0,
    retries   = 10,
    block     = True
)
# utilities
import re, click
import yaml as yml
# astropy
from astropy.io import fits
import astropy.units as u
# threading
import threading, queue
# configuration
from surveys.survey_config import SurveyConfig
# processing
from surveys.survey_abc import processing_status as ProcStatus


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
    return target


# save an HDU into a file
def save_cutout(target):
    def pack_message(msg,msg_log):
        p="\nLOG: "
        return re.sub(r"^\n","",((f"{p}"+f"{p}".join(msg_log.splitlines())) if msg_log != "" else "")+f"{p}{msg}")

    #print("SAVING CUTOUT ")
    if target['hdu']:
        target['hdu'].writeto("{0}".format(target['filename']), overwrite=True)
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


# TODO: May not need this...
# define the default config file with absolute path
this_source_file_dir = re.sub(r"(.*/).*$",r"\1",os.path.realpath(__file__))
default_config = this_source_file_dir + 'config.yml'


@click.group()
def cli():
    """\b
       Survey cutout fetching script.
       Command help: <command> --help
    """
    pass

@cli.command()
@click.argument('target')
@click.argument('size')
@click.argument('surveys', required=False, nargs=-1)
def fetch(target, size, surveys):
    """
    fetch cutouts from single coordinate, 'target'
    'size' is the image width in arcmin as an integer
    specify one or several surveys comma separated
    if surveys argument not specified then fetch from all implemented surveys
    """
    print("target", target)
    print("size", size)
    print(surveys)
    surveys = list(surveys)
    print("surveys", surveys)
    # if not surveys:
    #     survey_list = []
    # else:
    #     survey_list = surveys.split(' ')
    # print("surveys parsed:", survey)

    # provide list of surveys rather than name of a config file to parse items from
    cfg = SurveyConfig(surveys)
    if not cfg:
        print("invalid or missing arguments! Request aborted")
        return
    cfg.set_single_target_params(target, size)
    
    print(cfg.size_arcmin)

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

@cli.command()
@click.argument('config-file')
def status(config_file):
    """
       Cutout pipeline status command.
       For CONFIG_FILE formating help details see:
          python fetch_cutouts.py batch-process --help
    """
    pass


@cli.command()
@click.option('--flush',is_flag=True,default=None,help='flush existing target files')
@click.argument('config-file')
def maintenance(
    config_file,
    flush
):
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
def batch_process(
    config_file,
    overwrite,
    flush
):
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

       where targets.csv must at least have columns of ra and dec, e.g.,\n
       \b
          id,RA,DEC
          6,162.33807373,-0.668059408665
          107,0.417553782463,1.09196579456
          191,132.387832642,1.7280052900299998
          777,176.34243774400002,15.4953451157
          986,208.63735961900002,28.2433853149
    """

    # load the configuration
    print(f"Using Configuration: {config_file}")
    cfg = SurveyConfig(config_file)
    print(f"Overwrite Mode: {cfg.set_overwrite(overwrite)}")
    print(f'Flush Mode: {cfg.set_flush(flush)}')

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


if __name__ == "__main__":
    #batch_process()
    cli()
