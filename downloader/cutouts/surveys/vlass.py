import os
import io
import sys
from pathlib import Path

import urllib

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy import units as u

from astroquery.cadc import Cadc


from .survey_abc import SurveyABC
class VLASS(SurveyABC):
    def __init__(self):
        super().__init__()

    @staticmethod
    def get_supported_filters():
        return None


    def get_filter_setting(self):
        return None


    def get_tile_urls(self,position,size):
        # Notes: upgrade to more robust method using JJ's recipe...
        #    1) pip install https://github.com/astropy/astroquery/archive/master.zip
        #    2) results = cadc.query_region(coords,
        #                 radius=value * u.deg,
        #                 collection='VLASS')
        #    3) cadc.get_image_list(query_result, coordinates, radius)
        # where,
        #    cadc = astroquery.cadc.Cadc()
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

