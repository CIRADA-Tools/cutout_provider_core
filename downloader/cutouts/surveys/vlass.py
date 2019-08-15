import os
import io
import sys
#import shutil
#import errno
#import tempfile
from pathlib import Path

import urllib

# deprecated 
#import csv
#import pickle

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy import units as u

#from astropy.nddata.utils import NoOverlapError

from astroquery.cadc import Cadc

# deprecated 
## tile_query() static function contructor...
#def tile_query():
#    # get vlass quick-look image info list
#    def get_subtiles():
#        script_location = Path(__file__).absolute().parent
#        try:
#            # load the pickle quick-look image info
#            with open(script_location/'vlass_tiles.pkl','rb') as infile:
#                tile_list = pickle.load(infile)
#        except IOError as e:
#            # whoops!
#            if e.errno == errno.ENOENT:
#                # csv info loader
#                def load_subtile_csv():
#                    script_location = Path(__file__).absolute().parent
#                    with open(script_location/'VLASS_subtiles.csv', 'r') as infile:
#                        reader = csv.DictReader(infile)
#                        extracted = ((tile['ra_centre'], tile['dec_centre'], tile['file']) for tile in reader)
#                        ra, dec, file = zip(*extracted)
#                        ra = Angle(list(ra), unit=u.hourangle)
#                        dec = Angle(list(dec), unit=u.deg)
#                        coords = SkyCoord(ra=ra, dec=dec)
#                        subtiles = {
#                            'files': file,
#                            'catalog': coords
#                        }
#                        return subtiles
#                # create the pickle from the vlass csv quick-look image info file
#                print("No sky coordinates")
#                print("Generating sky coordinates (should only need to do this once)")
#                tile_list = load_subtile_csv()
#                print("Saving sky coordinates")
#                with open(script_location/'vlass_tiles.pkl', 'wb+') as outfile:
#                    pickle.dump(tile_list, outfile)
#            else:
#                raise
#        return tile_list
#
#    # define closest_tiles file function
#    def closest_tiles(coord, tile_list, coord_size=5*u.arcmin, max_results=None):
#        mask = coord.separation(tile_list['catalog']) < 0.5*u.deg + coord_size
#        indices = (i for i, m in enumerate(mask) if m)
#        results = [tile_list['files'][i].split('/')[-1] for i in indices][0:max_results]
#        if len(results) > 1:
#             import re
#             sexadecimal = "%02d%02d%02.1f" % coord.ra.hms+re.sub(r"([+-])\d",r"\1","%+d%02d%02d%02.0f" % coord.dec.signed_dms)
#             print(f"VLASS: Mosaicking J{sexadecimal} of Size {coord_size} /w {len(results)} tiles...")
#        elif len(results) == 1:
#             print(f"VLASS: Only 1 tile to cut -- I love simplicity!")
#        else:
#             print(f"VLASS: Not tiles found for ({coord.ra.to(u.hms)}, {coord.dec})")
#        return results
#
#    # return a closer of tite_query() with loaded vlass title information
#    def query(coord, *args, **kwargs):
#        return closest_tiles(coord, get_subtiles(), *args, **kwargs)
#
#    return query


# deprecated 
## TODO (Issue #11): No longer thread safe -- fix!
## instatiate the title_query() function
#intersecting_tiles = tile_query()


from .survey_abc import SurveyABC
class VLASS(SurveyABC):
    # deprecated 
    #def __init__(
    #    self,
    #    is_cutout_server=True, # True use Cadc I/F; False Use VLASS Quick-Look Images
    #):
    def __init__(self):
        super().__init__()

        # deprecated 
        #self.is_cutout_server   = is_cutout_server
        #
        #self.intersecting_tiles = intersecting_tiles

        # TODO (Issue #11): This gets printed for each instance... mabye add a method and call it once... 
        #       hmmm... get_init_message_string()...? ... can use in survey_config... 
        #       avec try ... except ...
        #if self.is_cutout_server:
        #    self.print("=> Using CADC cutout server!")
        #else:
        #    self.print("=> Using VLASS quick-look images!")


    # deprecated 
    ## Greg Sivakoff's CADC cutout server url script...
    #def __get_cadc_cutout_urls(self,position,size):
    #    def construct_cadc_url(baseurl, position, radius):
    #        ICRS_position = position.transform_to('icrs')
    #        basefile = baseurl.split('pub/')[1].split('?')[0]
    #        if (basefile[-10:] == 'subim.fits' and basefile[:6] == 'VLASS/'):
    #            url = ( 'https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/cutout?uri=ad:' 
    #                   + urllib.parse.quote(basefile) 
    #                   + ('&cutout=Circle+ICRS+{}+{}+{}').format(ICRS_position.ra.degree,
    #                                                             ICRS_position.dec.degree,
    #                                                             radius.to(u.degree).value))
    #            return url
    #        else:
    #            self.print('CADC URL appears to be incorrect: {}'.format(basefile))
    #            return None
    #    cadc = Cadc()
    #    radius = size/np.sqrt(2.0)
    #    metadata = cadc.query_region(
    #        coordinates = position,
    #        radius      = radius.to(u.deg).value,
    #        collection  = 'VLASS'
    #    )
    #    if len(metadata) == 0:
    #        return list()
    #    base_urls = cadc.get_data_urls(metadata)
    #    urls = [construct_cadc_url(base_url, position, radius) for base_url in base_urls]
    #    return urls
    #
    #
    ## Michael Ramsay's VLASS quick-look image url script...
    #def __get_vlass_quick_look_image_urls(self,position,size):
    #    def get_query_url(tilename, position, size):
    #        url = 'http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/cutout'
    #        pos = 'Circle ICRS {0} {1} {2}'.format(
    #            position.ra.to_string(unit=u.deg, decimal=True),
    #            position.dec.to_string(unit=u.deg, decimal=True),
    #            size.to(u.deg).value
    #        )
    #        query_dict = {
    #            'uri': "ad:VLASS/{0}".format(tilename),
    #            'cutout': pos
    #        }
    #        query = urllib.parse.urlencode(query_dict)
    #        return "{url}?{query_string}".format(url=url, query_string=query)
    #
    #    tiles = self.intersecting_tiles(position, size)  # the actual file names
    #    if len(tiles)==0:
    #        #self.print("Cannot find {}, perhaps this hasn't been covered by VLASS".format(position.to_string('hmsdms')), file=sys.stderr)
    #        self.print("Cannot find {position.to_string('hmsdms')}, perhaps this hasn't been covered by VLASS.")
    #        return list()
    #    urls = [get_query_url(tile, position, size) for tile in tiles]
    #    return urls


    @staticmethod
    def get_supported_filters():
        return None


    def get_filter_setting(self):
        return None


    # deprecated 
    #def get_tile_urls(self,position,size):
    #    if self.is_cutout_server:
    #        urls = self.__get_cadc_cutout_urls(position,size)
    #    else:
    #        urls = self.__get_vlass_quick_look_image_urls(position,size)
    #    if len(urls)==0 and self.is_cutout_server:
    #        #self.print("Cannot find {}, perhaps this hasn't been covered by VLASS".format(position.to_string('hmsdms')), file=sys.stderr)
    #        self.print("Cannot find {position.to_string('hmsdms')}, perhaps this hasn't been covered by VLASS.")
    #    return urls
    def get_tile_urls(self,position,size):
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
                self.print('CADC URL appears to be incorrect: {}'.format(basefile))
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
        if len(urls)==0:
            self.print("Cannot find {position.to_string('hmsdms')}, perhaps this hasn't been covered by VLASS.")
        return urls


    def get_fits_header_updates(self,header,position,size):
        ###complex file name - extract from header info
        fpartkeys = [f'FILNAM{i+1:02}' for i in range(12)]
        nameparts = [header[key] for key in fpartkeys]
        ###create single string - FILNAM12 goes after a constant
        vfile = nameparts[0]
        for i in range(len(nameparts)-2):
            vfile = vfile + '.' + nameparts[i+1]
        vfile = vfile + '.pbcor.' + nameparts[len(nameparts)-1] + '.subim.fits'
        header_updates = {
            'BAND': ('2-4 GHz', 'Frequency coverage of observation'),
            'RADESYS':  (header['RADESYS'], 'Coordinate system used'),
            'DATE-OBS': (header['DATE-OBS'], 'Obs. date'),
            # TODO (Issue #6): same as MJD-OB... depricate.
            #'MJD': (Time(header['DATE-OBS']).mjd, 'MJD of the observation date'),
            'BUNIT': ('Jy/beam', 'Pixel flux unit'),
            'BMAJ':  (header['BMAJ'], 'Beam major axis [deg]'),
            'BMIN':  (header['BMIN'], 'Beam minor axis [deg]'),
            'BPA':   (header['BPA'], 'Beam position angle'),
            # TODO (Issue #6): might be already in wcs part of header...
            'STOKES': (header['BTYPE'], 'Stokes polarisation'),
            # TODO (Issue #6): Tiling issue and based on quick-look images -- I think...
            'IMFILE': (vfile, 'VLASS image file'),
        }
        return header_updates

