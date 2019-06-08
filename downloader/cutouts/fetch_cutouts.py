# system
import os
import sys
import signal

# utilities
import re
import click
import yaml as yml

# astropy
from astropy.io import fits

# threading
import threading
import queue

# configuration
from surveys.survey_config import SurveyConfig


def set_sig_handler(threads):
    def sig_handler(sig, frame):
        #signal.signal(signal.SIGINT, original_sigint)
        print(" ".join("  ***  CTRL-C  RECEIVED! KILLING THREADS... ***"))
        for t in threads:
            t.die()
        sys.exit(0)
    signal.signal(signal.SIGINT,sig_handler)
    

# kills a thread when given into the queue
class PoisonPill:
    def __init__(self):
        pass


class WorkerThread(threading.Thread):
    def __init__(self, work, input_q, output_q=None, *args, **kwargs):
        self.input_q = input_q
        self.output_q = output_q
        self.work = work
        self.kill_recieved = False
        super().__init__(*args, **kwargs)

    def run(self):
        while not self.kill_recieved:
            work_in = self.input_q.get()

            # if it's swallowed
            if type(work_in) is PoisonPill:
                self.input_q.task_done()
                return

            ret = self.work(work_in)
            self.input_q.task_done()

            if self.output_q:
                self.output_q.put(item=ret)

        if self.kill_recieved:
            try:
                survey = work_in['survey']
                survey.print(' '.join('SHUTTING DOWN GRACEFULLY...'))
                del survey
            except:
                print("Bye!")
            del self.input_q
            if self.output_q:
                del self.output_q

    def die(self):
        self.kill_recieved = True


# grab a FITS hdu from some survey
def get_cutout(target):
    target['hdu'] = target['survey'].get_cutout(target['coord'], target['size'])
    return target


# save an HDU into a file
def save_cutout(target):
    if target['hdu']:
        target['hdu'].writeto("{0}".format(target['filename']), overwrite=True)
    else:
        survey = type(target['survey']).__name__
        msg_str = f"cutout at {target['coord']} returned None"
        prefix_msg_str = "\n".join([f"{survey}: {s}" for s in msg_str.splitlines()])
        print(prefix_msg_str)


@click.command()
@click.option('--config-file',default='config_default.yml',help='yaml search parameters configuration file')
def batch_process(config_file="config_default.yml"):
    """Survey Cutout fetching script (cf., config_default.yml)"""

    # load yaml configuration file
    print(f"Using Configuration: {config_file}")
    cfg = SurveyConfig(config_file)

    grabbers = 10
    savers = 1

    # set up i/o queues
    in_q  = queue.Queue()
    out_q = queue.Queue()

    # toss all the targets into the queue, including for all surveys
    # i.e., some position in both NVSS and VLASS and SDSS, etc.
    for task in cfg.get_procssing_stack():
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
    batch_process()
