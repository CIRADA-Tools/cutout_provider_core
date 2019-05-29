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
    list1 = csv_to_dict('BH12_multicomponent_zlt0p05_targetlist2.csv')
    list2 = csv_to_dict('BH12_120RA180_minus4DEC16_RC4_targetlist.csv')

    targets = [
        {
            'coord': SkyCoord(x['RA'], x['Dec'], unit=(u.deg, u.deg)),
            'size': 5*u.arcmin
        }
        for x in list1+list2]

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
    for idx, target in enumerate(get_target_list()):

        for s in all_surveys:

            t = dict(target)

            t['survey'] = s
            t['filename'] = "data/out/{1}_{0}.fits".format(t['survey'].__name__, idx)

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
