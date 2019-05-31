import re
import os

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


def get_target_list():

    # initial list of 17 galaxies
    #list1 = csv_to_dict('BH12_multicomponent_zlt0p05_targetlist2.csv')
    #list2 = csv_to_dict('BH12_120RA180_minus4DEC16_RC4_targetlist.csv')
    #sources = list1+list2

    sources =  csv_to_dict('targs_for_michelle.csv')

    # make all keys lower case to effectively reference case-variants of RA and Dec.
    sources = [{k.lower(): v for k,v in s.items()} for s in sources]

    # extract position information
    targets = [
        {
            'coord': SkyCoord(x['ra'], x['dec'], unit=(u.deg, u.deg)),
            'size': 5*u.arcmin
        }
        for x in sources]

    return targets


# grab a FITS hdu from some survey
def get_cutout(target):

    target['hdu'] = target['survey'].get_cutout(target['coord'], target['size'])
    return target


# save an HDU into a file
def save_cutout(target):

    if target['hdu']:
        target['hdu'].writeto("{0}".format(target['filename']), overwrite=True)
    else:
        print("{0} cutout at {1} returned None".format(target['survey'].__name__, target['coord']), sys.stderr)


def batch_process():

    out_dir = "data/out"
    try:
        os.makedirs(out_dir)
    except FileExistsError:
        print(f"Using FITS output dir: {out_dir}")
    else:
        print(f"Created FITS output dir: {out_dir}")

    grabbers = 10
    savers = 1

    all_surveys = (

        NVSS,
        FIRST,
        SDSS,
        VLASS,
        PanSTARRS

    )

    in_q = queue.Queue()
    out_q = queue.Queue()

    # toss all the targets into the queue, including for all surveys
    # i.e some position in both NVSS and VLASS and SDSS, etc.
    targets = get_target_list()
    zero_padding = len(f"{len(targets)}")
    for idx, target in enumerate(targets):

        for s in all_surveys:

            t = dict(target)

            t['survey'] = s
            sexadecimal = "%02d%02d%02.1f" % t['coord'].ra.hms+re.sub(r"([+-])\d",r"\1","%+d%02d%02d%02.0f" % t['coord'].dec.signed_dms)
            size = re.sub(r"\.?0+$","","%f" % t['size'].value)
            t['filename'] = f"{out_dir}/J{sexadecimal}_s{size}_{t['survey'].__name__}.fits"

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
