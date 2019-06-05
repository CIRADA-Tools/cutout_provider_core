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

from astroquery.cadc import Cadc

# TODO: factor out
import montage_wrapper


def load_subtile_csv():

    script_location = Path(__file__).absolute().parent

    with open(script_location/'VLASS_subtiles.csv', 'r') as infile:
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

    if len(results) > 1:
         import re
         sexadecimal = "%02d%02d%02.1f" % coord.ra.hms+re.sub(r"([+-])\d",r"\1","%+d%02d%02d%02.0f" % coord.dec.signed_dms)
         print(f"VLASS: Mosaicing J{sexadecimal} of Size {coord_size} /w {len(results)} tiles...")
    elif len(results) == 1:
         print(f"VLASS: Only 1 tile to cut -- I love simplicity!")
    else:
         print(f"VLASS: Not tiles found for ({coord.ra.to(u.hms)}, {coord.dec})")

    return results


def tile_query(tile_list):

    def query(coord, *args, **kwargs):
        return closest_tiles(coord, tile_list, *args, **kwargs)

    return query


intersecting_tiles = tile_query(get_subtiles())

# Greg Sivalkoff patch...
def get_tile_urls(position,size):
    def construct_cadc_url(baseurl, position, radius):
        ICRS_position = position.transform_to('icrs')
        basefile = baseurl.split('pub/')[1].split('?')[0]
        if (basefile[-10:] == 'subim.fits' and basefile[:6] == 'VLASS/'):
            url = ( 'https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/cutout?uri=ad:' 
                   + urllib.parse.quote(basefile) 
                   + ('&cutout=Circle+ICRS+{}+{}+{}').format(ICRS_position.ra.degree,
                                                             ICRS_position.dec.degree,
                                                             radius.to(u.degree).value))
            return url
        else:
            print('CADC URL appears to be incorrect: {}'.format(basefile))
            return None
    cadc = Cadc()
    radius = size/np.sqrt(2.0)
    metadata = cadc.query_region(
        coordinates = position,
        radius      = radius.to(u.deg).value,
        collection  = 'VLASS'
    )
    if len(metadata) == 0:
        return list()
    base_urls = cadc.get_data_urls(metadata)
    urls = [construct_cadc_url(base_url, position, radius) for base_url in base_urls]
    return urls


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

def squeeze(hdu):

    w = WCS(hdu.header)

    # trim to 2d from 4d
    w = w.dropaxis(2).dropaxis(2)
    img_data = np.squeeze(hdu.data)

    img = fits.PrimaryHDU(img_data, header=hdu.header)

    # writing to a pretend file in memory
    mem_file = io.BytesIO()
    img.writeto(mem_file)
    return mem_file.getvalue()


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


from .survey import Survey
class VLASS(Survey):
    def __init__(self,is_cutout_server=True):
        super().__init__()
        self.is_cutout_server = is_cutout_server

        if self.is_cutout_server:
            print("=> Using CADC cutout server!")
        else:
            print("=> Using VLASS quick-look images!")


    def get_tile_cutouts(self, position, size):
        if position.dec.value > 89:
            print("Warning: cutouts near the poles may give unexpected results/fail", file=sys.stderr)

        if self.is_cutout_server:
            urls = get_tile_urls(position,size)
            if len(urls)>0:
                hdu_lists = [h for h in [self.get_fits(url) for url in urls] if h]
            else:
                print("Cannot find {}, perhaps this hasn't been covered by VLASS".format(position.to_string('hmsdms')), file=sys.stderr)
                return list() 
        else:
            tiles = intersecting_tiles(position, size)  # the actual file names
            urls = [get_query_url(tile, position, size) for tile in tiles]
            if tiles:
                hdu_lists = [h for h in [self.get_fits(url) for url in urls] if h]
            else:
                print("Cannot find {}, perhaps this hasn't been covered by VLASS".format(position.to_string('hmsdms')), file=sys.stderr)
                return list()

        cutouts = list()
        for hdu_list in hdu_lists:
            try:
                img = squeeze(hdu_list[0]) if self.is_cutout_server else cutout(hdu_list[0], position, size)
                cutouts.append(img)
            except NoOverlapError as e:
                if self.is_cutout_server:
                    import re
                    sexadecimal = "%02d%02d%02.0f" % position.ra.hms+re.sub(r"([+-])\d",r"\1","%+d%02d%02d%02.0f" % position.dec.signed_dms)
                    print(f"Overlap Error: {sexadecimal} => {urls}")
                else:
                    # TODO: The original script produces a this error, but ignores it... prelminary exploration 
                    #       seems to indicate it's because the cutouts overlap with no excess... so perhaps its
                    #       OK, to trim the cutouts servers stuff to taste using cutout(hdu_list[0], position, size), 
                    #       instead of squeeze(hdu_list[0]) on self.is_cutout_server == True; effectly it would trim
                    #       fat caused by the mapping size -> radius = size/sqrt(2).
                    pass

        return cutouts


    def get_cutout(self,position, size):
        cutouts = self.get_tile_cutouts(position,size)

        if len(cutouts)==0:
            return None

        if len(cutouts) > 1:
            try:
                c = self.mosaic(cutouts)
            except montage_wrapper.status.MontageError as e:
                print(e, file=sys.stderr)
                return None
        else:
            try:
                c = cutouts[0]
            except IndexError:
                return None

        return self.create_fits(c)

