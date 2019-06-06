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

from astropy.nddata.utils import NoOverlapError

from astroquery.cadc import Cadc

# TODO: factor out
import montage_wrapper


# tile_query() static function contructor...
def tile_query():
    # get vlass quick-look image info list
    def get_subtiles():
        script_location = Path(__file__).absolute().parent
        try:
            # load the pickle quick-look image info
            with open(script_location/'vlass_tiles.pkl','rb') as infile:
                tile_list = pickle.load(infile)
        except IOError as e:
            # whoops!
            if e.errno == errno.ENOENT:
                # csv info loader
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
                # create the pickle from the vlass csv quick-look image info file
                print("No sky coordinates")
                print("Generating sky coordinates (should only need to do this once)")
                tile_list = load_subtile_csv()
                print("Saving sky coordinates")
                with open(script_location/'vlass_tiles.pkl', 'wb+') as outfile:
                    pickle.dump(tile_list, outfile)
            else:
                raise
        return tile_list

    # define closest_tiles file function
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

    # return a closer of tite_query() with loaded vlass title information
    def query(coord, *args, **kwargs):
        return closest_tiles(coord, get_subtiles(), *args, **kwargs)

    return query

# TODO: No longer thread safe -- fix!
# instatiate the title_query() function
intersecting_tiles = tile_query()


from .survey import Survey
class VLASS(Survey):
    def __init__(self,is_cutout_server=True):
        super().__init__()

        self.is_cutout_server = is_cutout_server
        self.intersecting_tiles = intersecting_tiles

        if self.is_cutout_server:
            print("=> Using CADC cutout server!")
        else:
            print("=> Using VLASS quick-look images!")


    # Greg Sivakoff's CADC cutout server url script...
    def __get_cadc_cutout_urls(self,position,size):
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


    # Micheal Ramsay's VLASS quick-look image url script...
    def __get_vlass_quick_look_image_urls(self,position,size):
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

        tiles = self.intersecting_tiles(position, size)  # the actual file names
        if len(tiles)==0:
            print("Cannot find {}, perhaps this hasn't been covered by VLASS".format(position.to_string('hmsdms')), file=sys.stderr)
            return list()
        urls = [get_query_url(tile, position, size) for tile in tiles]
        return urls


    def get_tile_urls(self,position,size):
        if self.is_cutout_server:
            urls = self.__get_cadc_cutout_urls(position,size)
        else:
            urls = self.__get_vlass_quick_look_image_urls(position,size)
        if len(urls)==0 and self.is_cutout_server:
            print("Cannot find {}, perhaps this hasn't been covered by VLASS".format(position.to_string('hmsdms')), file=sys.stderr)
        return urls


    #            * * * D E P R E C A T E D * * *
    #def get_tiles(self,position,size):
    #    urls = self.get_tile_urls(position,size)
    #    hdu_lists = [h for h in [self.get_fits(url) for url in urls] if h]
    #    return hdu_lists
    #
    #
    #def paste_tiles(self,hdu_tiles,position,size):
    #    img = None
    #    if len(hdu_tiles) > 1:
    #        try:
    #            imgs = [img for img in [self.get_image(tile) for tile in hdu_tiles]]
    #            img = self.create_fits(self.mosaic(imgs))
    #        except montage_wrapper.status.MontageError as e:
    #            print(f"Mosaicing Failed: {e}: file={sys.stderr}")
    #    elif len(hdu_tiles) == 1:
    #            img = hdu_tiles[0]
    #    #if img:
    #    #    try:
    #    #        img = self.trim_tile(img,position,size)
    #    #    except Exception as e:
    #    #        print(f"Trim Failed: {e}")
    #    return img
    #
    #
    ##            * * * D E P R E C A T E D * * *
    ##def get_tile_cutouts(self, position, size):
    ##    if position.dec.value > 89:
    ##        print("Warning: cutouts near the poles may give unexpected results/fail", file=sys.stderr)
    ##
    ##    hdu_lists = self.get_tiles(position,size)
    ##
    ##    cutouts = list()
    ##    for hdu_list in hdu_lists:
    ##        try:
    ##            img = self.get_image(hdu_list) if self.is_cutout_server else self.cutout(hdu_list[0], position, size)
    ##            cutouts.append(img)
    ##        except NoOverlapError as e:
    ##            if self.is_cutout_server:
    ##                import re
    ##                #sexadecimal = "%02d%02d%02.0f" % position.ra.hms+re.sub(r"([+-])\d",r"\1","%+d%02d%02d%02.0f" % position.dec.signed_dms)
    ##                print(f"Overlap Error: {self.get_sexy_string(position)} => {urls}")
    ##            else:
    ##                # TODO: The original script produces a this error, but ignores it... prelminary exploration 
    ##                #       seems to indicate it's because the cutouts overlap with no excess... so perhaps its
    ##                #       OK, to trim the cutouts servers stuff to taste using cutout(hdu_list[0], position, size), 
    ##                #       instead of self.get_image(hdu_list[0]) on self.is_cutout_server == True; effectly it would
    ##                #       trim fat caused by the mapping size -> radius = size/sqrt(2).
    ##                #
    ##                #       Need to consult with Michael Ramsay (original author).
    ##                #
    ##                #       Update: See (*) TODO. (In consultaion with M Ramsay Jun 5, 19.)
    ##                pass
    ##
    ##    return cutouts
    ##
    ##
    ##def get_cutout(self,position, size):
    ##    cutouts = self.get_tile_cutouts(position,size)
    ##
    ##    if len(cutouts)==0:
    ##        return None
    ##
    ##    if len(cutouts) > 1:
    ##        try:
    ##            c = self.mosaic(cutouts)
    ##        except montage_wrapper.status.MontageError as e:
    ##            print(e, file=sys.stderr)
    ##            return None
    ##    else:
    ##        try:
    ##            c = cutouts[0]
    ##        except IndexError:
    ##            return None
    ##
    ##    return self.create_fits(c)
    ## (*) TODO: replacement process... that is, first mosiac then cut...
    #def get_cutout(self,position, size):
    #    return self.paste_tiles(self.get_tiles(position,size),position,size)

