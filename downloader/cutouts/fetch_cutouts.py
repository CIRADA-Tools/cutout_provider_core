import re
import os

import yaml as yml
import click

# An example of how to use the survey cutout code
import csv
import sys

import threading
import queue

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u

from surveys.nvss import NVSS
from surveys.first import FIRST
from surveys.wise import WISE
from surveys.sdss import SDSS
from surveys.vlass import VLASS
from surveys.panstarrs import PanSTARRS


# kills a thread when given into the queue
class PoisonPill:
    def __init__(self):
        pass


# A thread that grabs data from a queue, processes, then optionally tosses into another queue
class WorkerThread(threading.Thread):
    def __init__(self, work, input_q, output_q=None, *args, **kwargs):

        self.input_q = input_q
        self.output_q = output_q
        self.work = work

        super().__init__(*args, **kwargs)

    def run(self):

        while True:

            work_in = self.input_q.get()

            # if it's swallowed
            if type(work_in) is PoisonPill:
                self.input_q.task_done()
                return

            ret = self.work(work_in)
            self.input_q.task_done()

            if self.output_q:
                self.output_q.put(item=ret)


def csv_to_dict(filename):
    entries = []

    with open(filename, 'r') as infile:
        c = csv.DictReader(infile)
        for entry in c:
            entries.append(entry)

    return entries


def get_target_list(config_file="config.yml"):
    config  = yml.load(open(config_file,'r'))['cutouts']

    targets = list()
    size = config['box_size_armin'] * u.arcmin
    for coord_csv_file in config['ra_deg_deg_csv_files']:
        sources = csv_to_dict(coord_csv_file)

        # make all keys lower case to effectively reference case-variants of RA and Dec.
        sources = [{k.lower(): v for k,v in s.items()} for s in sources]

        # extract position information
        targets.extend([
            {
                'coord': SkyCoord(x['ra'], x['dec'], unit=(u.deg, u.deg)),
                'size': size
            }
            for x in sources])

    return targets

def get_surveys(config_file="config.yml"):
    all_surveys = (
        FIRST.__name__,
        NVSS.__name__,
        VLASS.__name__,
        WISE.__name__,
        PanSTARRS.__name__,
        SDSS.__name__,
    )
        
    requested = yml.load(open(config_file))['cutouts']['surveys']

    surveys = list()
    for survey in all_surveys:
        for request in requested:
            if request.lower() == survey.lower():
                surveys.append(eval("%s()" % survey))
                break

    return set(surveys)


# grab a FITS hdu from some survey
def get_cutout(target):

    target['hdu'] = target['survey'].get_cutout(target['coord'], target['size'])
    return target


# save an HDU into a file
def save_cutout(target):

    if target['hdu']:
        target['hdu'].writeto("{0}".format(target['filename']), overwrite=True)
    else:
        #print("{0}: cutout at {1} returned None: {2}".format(survey, target['coord'], sys.stderr))
        survey = type(target['survey']).__name__
        msg_str = f"cutout at {target['coord']} returned None"
        prefix_msg_str = "\n".join([f"{survey}: {s}" for s in msg_str.splitlines()])
        print(prefix_msg_str)

@click.command()
@click.option('--config-file',default='config.yml',help='yaml search parameters configuration file')
def batch_process(config_file="config.yml"):
    """Suvery Cutout fetching script (cf., config.yml)"""

    print(f"Using Configuration: {config_file}")

    out_dir = "data/out"
    try:
        os.makedirs(out_dir)
    except FileExistsError:
        print(f"Using FITS output dir: {out_dir}")
    else:
        print(f"Created FITS output dir: {out_dir}")

    grabbers = 10
    savers = 1

    all_surveys = get_surveys(config_file)

    in_q = queue.Queue()
    out_q = queue.Queue()

    # toss all the targets into the queue, including for all surveys
    # i.e some position in both NVSS and VLASS and SDSS, etc.
    targets = get_target_list(config_file)
    zero_padding = len(f"{len(targets)}")
    for idx, target in enumerate(targets):

        for s in all_surveys:

            t = dict(target)

            t['survey'] = s
            size = re.sub(r"\.?0+$","","%f" % t['size'].value)
            t['filename'] = f"{out_dir}/J{s.get_sexy_string(t['coord'])}_s{size}arcmin_{type(t['survey']).__name__}.fits"

            in_q.put(t)

    # spin up a bunch of worker threads to process all the data
    # in principle these could be chained further, such that you could go
    # targets -> hdus -> save to file -> process to jpg -> save to file
    for _ in range(grabbers):

        WorkerThread(get_cutout, in_q, out_q).start()
        in_q.put(PoisonPill())

    # testing out 1 save to file threads (absolutely not necessary)
    for _ in range(savers):
        WorkerThread(save_cutout, out_q).start()
    in_q.join()

    for _ in range(savers):
        out_q.put(PoisonPill())
    out_q.join()


if __name__ == "__main__":
    batch_process()
