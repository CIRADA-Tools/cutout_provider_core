import re

from astropy.table import Table
from astropy import units as u
from astroquery.ibe import IbeClass


from .survey_abc import SurveyABC
from .survey_filters import wise_filters
class WISE(SurveyABC):
    def __init__(self,filter=wise_filters.w1):
        super().__init__()
        self.filter = filter
        self.metadata_root = 'p3am_cdd'
        self.url_root = f"https://irsa.ipac.caltech.edu/ibe/data/wise/allwise/{self.metadata_root}"

    @staticmethod
    def get_supported_filters():
        return wise_filters

    def get_filter_setting(self):
        return self.filter

    def __get_coadd_ids(self,metadata):
        if len(metadata) > 0:
           coadd_ids = metadata[metadata['band']==self.filter.value]['coadd_id']
           if len(coadd_ids) > 0:
               return coadd_ids
        #return Table()
        return list()

    def __get_fits_urls(self,coadd_ids):
        urls = list()
        for coadd_id in coadd_ids:
            coaddgrp = coadd_id[:2]
            coadd_ra = coadd_id[:4]
            urls.append(f"{self.url_root}/{coaddgrp:s}/{coadd_ra:s}/{coadd_id:s}/{coadd_id}-{self.filter.name}-int-3.fits")
        return urls

    def get_tile_urls(self,position,size):
        #status = list()
        wise = IbeClass()
        edge = size.to(u.deg)
        metadata = wise.query_region(
            coordinate = position,
            mission    = 'wise',
            dataset    = 'allwise',
            table      = self.metadata_root,
            columns    = 'band,coadd_id',
            width      = str(edge),
            # height     = edge, # makes square by default
            intersect  = 'COVERS'
        )
        if len(metadata)==0:
            self.print(f"Position ({position.ra}, {position.dec}) has overlapping tiles.")
        coadd_ids = self.__get_coadd_ids(metadata)
        fits_urls = self.__get_fits_urls(coadd_ids)
        return fits_urls

    def get_fits_header_updates(self,header, all_headers=None):
        header_updates = {
            'BAND':     (f'{self.filter.name.upper()}', 'Filter used in observation'),
            'DATE-OBS': (header['MIDOBS'], 'Median observation date of stack'),
             # TODO (Issue #6): 'IMFILE' inadequate for mosaics...
            'IMFILE':   (header['COADDID'], 'ATLAS image identifier')
        }
        return header_updates
