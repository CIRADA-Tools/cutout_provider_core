from astropy import units as u
from astropy.table import Table

from astroquery.ibe import IbeClass

from .survey import Survey
from .fits_request import get_fits

from enum import Enum
class wiseFilters(Enum):
    w1 = 1
    w2 = 2
    w3 = 3
    w4 = 4

class WISE(Survey):
    def __init__(self,filter=wiseFilters.w1):
        self.filter = filter
        self.metadata_root = 'p3am_cdd'
        self.url_root = f"https://irsa.ipac.caltech.edu/ibe/data/wise/allwise/{self.metadata_root}"

    def __get_coadd_ids(self,metadata):
        if len(metadata) > 0:
           coadd_ids = metadata[metadata['band']==self.filter.value]['coadd_id']
           if len(coadd_ids) > 0:
               return coadd_ids
        return Table()

    def __get_fits_urls(self,coadd_ids):
        urls = list()
        for coadd_id in coadd_ids:
            coaddgrp = coadd_id[:2]
            coadd_ra = coadd_id[:4]
            urls.append(f"{self.url_root}/{coaddgrp:s}/{coadd_ra:s}/{coadd_id:s}/{coadd_id}-{self.filter.name}-int-3.fits")
        return urls

    def get_cutout(self,position, size):

        wise = IbeClass()

        edge = size.to(u.deg)
        metadata = wise.query_region(
            coordinate = position,
            mission    = 'wise',
            dataset    = 'allwise',
            table      = self.metadata_root,
            width      = edge,
            height     = edge,
            intersect  = 'COVERS',
            most_centered = True
        )

        coadd_ids = self.__get_coadd_ids(metadata)
        fits_urls = self.__get_fits_urls(coadd_ids)
        for i in fits_urls:
            print(f"> {i}")
        #if len(fits_urls)==0:
        #    print("> https://WHOOPS!")
        #    print(metadata)
        #    print(position)

        return get_fits(fits_urls[0]) if len(fits_urls) > 0 else None


