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

    # this will work for ANY collection from CADC
    def get_tile_urls(self,position,size):
        cadc = Cadc()
        radius = (size/2.0).to(u.deg)
        results = cadc.query_region(
            coordinates = position,
            radius      = radius,
            collection  = 'VLASS'
        )
        ### If adding any filters in then this is where would do it!!!#####
        #### e.g. filtered_results = results[results['time_exposure'] > 120.0] #####
        if len(results) == 0:
            return list()
        urls = cadc.get_image_list(results, position, radius)
        if len(urls)==0:
            self.print("Cannot find {position.to_string('hmsdms')}, perhaps this hasn't been covered by VLASS.")
        return urls

    def get_fits_header_updates(self,header, all_headers=None):
        self.print("header updates")
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
            'BTYPE': (header['BTYPE'], 'Stokes polarisation'),
            # TODO (Issue #6): Tiling issue and based on quick-look images -- I think...
            # 'IMFILE': (vfile, 'VLASS image file'),
            'COMMENT': "Quick Look images do not fully sample the PSF, and are cleaned to a threshold " \
                       " of ~5 sigma (details can be found in the weblogs for individual images). " \
                       "They are used for Quality Assurance and for transient searches, but should not " \
                       "be used for any other purpose. In addition to imaging artifacts, source " \
                       "positions can be off by up to 1-arcsec, and the flux density uncertainties " \
                       "are ~10-20%.",
        }
        ### ONLY FOR MOSAICKED. complex file name list all originals gone into mosaic
        if all_headers:
            for num,head in enumerate(all_headers):
                fpartkeys = [f'FILNAM{i+1:02}' for i in range(12)]
                nameparts = [head[key] for key in fpartkeys]
                header_updates['IMFILE'+str(num+1).zfill(2)]='.'.join(nameparts) + '.subim.fits'
                ###create single string - FILNAM12 goes after a constant
        return header_updates
