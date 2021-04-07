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
from .survey_filters import vlass_epoch


class VLASS(SurveyABC):
    def __init__(self, filter=None):
        super().__init__()
        self.needs_trimming = False
        self.filter = filter

    @staticmethod
    def get_supported_filters():
        return vlass_epoch

    @staticmethod
    def get_epoch(fileOrURL):
        # works for VLASS quicklook images
        # extracts '1.1' from url e.g. https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/sync?ID=ad%3AVLASS%2FVLASS1.1.ql.T11t36.J235800%2B003000.10.2048.v1.I.iter1.image.pbcor.tt0.subim.fits&RUNID=x9l8m9ed6s139yir&POS=CIRCLE+6.361111354985913e-06+2.5277771640709656e-06+0.08333333333333333
        # or filename e.g. VLASS__1.1.ql.T11t01.J000000+000000.10.2048.v1.I.iter1.image.pbcor.tt0.subim_s3.0arcmin.fits
        return fileOrURL.split('.ql')[0][-3:]

    @staticmethod
    # based on larger QL url
    def get_cutout_url(ql_url,coords, radius):
        standard_front = 'https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/sync?ID=ad%3AVLASS%2F'
        # '?RUNID' messes everything up, don't need it
        encoded_ql = urllib.parse.quote(ql_url.split("/")[-1])
        #encoded_ql= urllib.parse.quote(ql_url.split("/")[-1])
        encoded_ql = encoded_ql.replace('%3F','&').replace('?','&')
        #safe_coords = SkyCoord(ra=self.component.RA, dec=self.component.DEC, unit='deg')
        cutout_end = f"&CIRCLE={coords.ra.value}+{coords.dec.value}+{radius.value}"
        return standard_front+encoded_ql+cutout_end

    def get_filter_setting(self):
        return self.filter

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
        urls = []
        # urls = cadc.get_images(
        #     coordinates = position,
        #     radius      = radius,
        #     collection  = 'VLASS',
        #     get_url_list= True
        # )

        #change to
        # Observation.obsID = VLASS2.1.T08t18.J112603-083000
        all_rows = cadc.exec_sync(f"Plane.publisherID, Observation.requirements_flag FROM caom2.Plane AS Plane JOIN caom2.Observation AS Observation \
                                    ON Plane.obsID = Observation.obsID WHERE  ( Observation.collection = 'VLASS' \
                                    AND INTERSECTS( CIRCLE('ICRS', {position.ra.value}, {position.dec.value},  {radius.value}), \
	                                      Plane.position_bounds ) = 1 \
                                          AND  ( Plane.quality_flag IS NULL OR Plane.quality_flag != 'junk' \
                                          OR Plane.quality_flag != 'junk') AND Observation.requirements_flag != 'fail' )")

        if len(all_rows)>0:
            for row in all_rows:
                ql_urls = self.cadc.get_data_urls(all_rows)
                if ql_urls:
                    for url in ql_urls:
                        cutout_url = VLASS.get_cutout_url(url[0], position, radius)
                        urls.append(cutout_url)
        ### If adding any filters in then this is where would do it!!!#####
        #### e.g. filtered_results = results[results['time_exposure'] > 120.0] #####

        # if len(results) == 0:
        #     return list()
        # print("or this one?")
        # urls = cadc.get_image_list(results, position, radius)
        # if len(urls)==0:
        #     self.print("Cannot find {position.to_string('hmsdms')}, perhaps this hasn't been covered by VLASS.")
        if self.filter:
            final_urls = []
            for url in urls:
                epoch = VLASS.get_epoch(url)
                if epoch==self.filter.value:
                    final_urls.append(url)
            urls = final_urls
        if len(urls) == 0:
            self.print("Cannot find {position.to_string('hmsdms')}, perhaps this hasn't been covered by VLASS.")
            return list()
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
