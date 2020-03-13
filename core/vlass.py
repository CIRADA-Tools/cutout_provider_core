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
from .toolbox import pad_string_lines


class VLASS(SurveyABC):
    def __init__(self):
        super().__init__()
        self.needs_trimming = False

    @staticmethod
    def get_supported_filters():
        return []

    @staticmethod
    def get_epoch(fileOrURL):
        # works for VLASS quicklook images
        # extracts '1.1' from url e.g. https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/sync?ID=ad%3AVLASS%2FVLASS1.1.ql.T11t36.J235800%2B003000.10.2048.v1.I.iter1.image.pbcor.tt0.subim.fits&RUNID=x9l8m9ed6s139yir&POS=CIRCLE+6.361111354985913e-06+2.5277771640709656e-06+0.08333333333333333
        # or filename e.g. VLASS__1.1.ql.T11t01.J000000+000000.10.2048.v1.I.iter1.image.pbcor.tt0.subim_s3.0arcmin.fits
        return fileOrURL.split('.ql')[0][-3:]

    def get_filter_setting(self):
        return None

    def add_cutout_service_comment(self, hdu):
        hdu.header.add_comment(pad_string_lines("Quick Look images do not fully sample the PSF, and are cleaned to a " \
                   "threshold of ~5 sigma (details can be found in the weblogs for " \
                   "individual images). They are used for Quality Assurance and for " \
                   "transient searches, but should not be used for any other purpose. In " \
                   "addition to imaging artifacts, source positions can be off by up to " \
                   "1-arcsec, and the flux density uncertainties are ~10-20%. \
                   "),after=-1)
        hdu.header.add_comment(pad_string_lines("The direct data service at CADC was used to provide this cutout: " \
                                " (https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/doc/data/) \
                                "), after=-1)

    # this will work for ANY collection from CADC
    def get_tile_urls(self,position,size):
        cadc = Cadc()
        radius = (size/2.0).to(u.deg)
        urls = cadc.get_images(
            coordinates = position,
            radius      = radius,
            collection  = 'VLASS',
            get_url_list= True
        )
        ### If adding any filters in then this is where would do it!!!#####
        #### e.g. filtered_results = results[results['time_exposure'] > 120.0] #####
        if len(urls) == 0:
            self.print("Cannot find {position.to_string('hmsdms')}, perhaps this hasn't been covered by VLASS.")
            return list()
        # if len(results) == 0:
        #     return list()
        # print("or this one?")
        # urls = cadc.get_image_list(results, position, radius)
        # if len(urls)==0:
        #     self.print("Cannot find {position.to_string('hmsdms')}, perhaps this hasn't been covered by VLASS.")
        return urls

    def get_fits_header_updates(self,header, all_headers=None):
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
        }
        ### ONLY FOR MOSAICKED. complex file name list all originals gone into mosaic
        if all_headers:
            for num,head in enumerate(all_headers):
                fpartkeys = [f'FILNAM{i+1:02}' for i in range(12)]
                nameparts = [head[key] for key in fpartkeys]
                header_updates['IMFILE'+str(num+1).zfill(2)]='.'.join(nameparts) + '.subim.fits'
                ###create single string - FILNAM12 goes after a constant
        return header_updates
