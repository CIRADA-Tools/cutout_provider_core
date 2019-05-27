import csv

import pickle
import tempfile
import shutil
import os
import errno
import io
import sys
from pathlib import Path

import urllib

import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
from astropy.nddata.utils import NoOverlapError

import montage_wrapper

from .survey import Survey
from .fits_request import get_fits
from .fits_request import create_fits


def load_subtile_csv():

    script_location = Path(__file__).absolute().parent

    with open(script_location/'vlass_subtiles.csv', 'r') as infile:
        reader = csv.DictReader(infile)

        extracted = ((tile['ra_centre'], tile['dec_centre'], tile['file']) for tile in reader)
        ra, dec, file = zip(*extracted)

        ra = Angle(list(ra), unit=u.hourangle)
        dec = Angle(list(dec), unit=u.deg)

        coords = SkyCoord(ra=ra, dec=dec)

        subtiles = {
            'files': file,
            'catalog': coords
        }

        return subtiles


def get_subtiles():

    script_location = Path(__file__).absolute().parent

    try:
        with open(script_location/'vlass_tiles.pkl', 'rb') as infile:
            tile_list = pickle.load(infile)

    except IOError as e:
        if e.errno == errno.ENOENT:
            print("No sky coordinates")
            print("Generating sky coordinates (should only need to do this once)")
            tile_list = load_subtile_csv()
            print("Saving sky coordinates")
            with open(script_location/'vlass_tiles.pkl', 'wb+') as outfile:
                pickle.dump(tile_list, outfile)

        else:
            raise

    return tile_list


def closest_tiles(coord, tile_list, coord_size=5*u.arcmin, max_results=None):

    mask = coord.separation(tile_list['catalog']) < 0.5*u.deg + coord_size

    indices = (i for i, m in enumerate(mask) if m)

    results = [tile_list['files'][i].split('/')[-1] for i in indices][0:max_results]

    return results


def tile_query(tile_list):

    def query(coord, *args, **kwargs):
        return closest_tiles(coord, tile_list, *args, **kwargs)

    return query


intersecting_tiles = tile_query(get_subtiles())


def cutout(hdu, position, size):

    w = WCS(hdu.header)

    # trim to 2d from 4d
    w = w.dropaxis(2).dropaxis(2)
    img_data = np.squeeze(hdu.data)

    stamp = Cutout2D(img_data, position, size, wcs=w, mode='partial', copy=True)
    header = stamp.wcs.to_header()
    img = fits.PrimaryHDU(stamp.data, header=header)

    # writing to a pretend file in memory
    mem_file = io.BytesIO()
    img.writeto(mem_file)
    return mem_file.getvalue()


# make the directory structure if it doesn't exist
def make_dir(dirname):

    try:
        os.makedirs(dirname)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def mosaic(cutouts):

    td = tempfile.mkdtemp()
    input_dir = '{directory}/input'.format(directory=td)
    output_dir = '{directory}/output'.format(directory=td)
    make_dir(input_dir)

    try:
        for i, c in enumerate(cutouts):

            with open('{directory}/{name}.fits'.format(directory=input_dir, name=i), 'wb') as tmp:
                tmp.write(bytes(c))
        os.listdir(td)
        montage_wrapper.mosaic(input_dir, output_dir)

        with open('{outdir}/mosaic.fits'.format(outdir=output_dir), 'rb') as f:

            merged = f.read()

    finally:
        shutil.rmtree(output_dir)
        shutil.rmtree(input_dir)
        shutil.rmtree(td)

    return merged


def get_query_url(tilename, position, size):

    url = 'http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/cutout'

    pos = 'Circle ICRS {0} {1} {2}'.format(
        position.ra.to_string(unit=u.deg, decimal=True),
        position.dec.to_string(unit=u.deg, decimal=True),
        size.to(u.deg).value
    )

    query_dict = {
        'uri': "ad:VLASS/{0}".format(tilename),
        'cutout': pos
    }

    query = urllib.parse.urlencode(query_dict)

    return "{url}?{query_string}".format(url=url, query_string=query)


class VLASS(Survey):

    @staticmethod
    def get_cutout(position, size):

        if position.dec.value > 89:
            print("Warning: cutouts near the poles may give unexpected results/fail", file=sys.stderr)

        tiles = intersecting_tiles(position, size)  # the actual file names

        urls = [get_query_url(tile, position, size) for tile in tiles]

        if tiles:
            hdu_lists = [h for h in [get_fits(url) for url in urls] if h]
        else:
            print("Cannot find {}, perhaps this hasn't been covered by VLASS".format(position.to_string('hmsdms')), file=sys.stderr)
            return None

        cutouts = []

        for hdu_list in hdu_lists:
            try:
                img = cutout(hdu_list[0], position, size)
                cutouts.append(img)
            except NoOverlapError as e:
                pass

        if len(cutouts) > 1:
            try:
                c = mosaic(cutouts)
            except montage_wrapper.status.MontageError as e:
                print(e, file=sys.stderr)
                return None
        else:
            try:
                c = cutouts[0]
            except IndexError:
                return None

        return create_fits(c)

