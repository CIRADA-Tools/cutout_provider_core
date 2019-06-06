import re

from astropy import units as u
from astropy.table import Table

from astroquery.ibe import IbeClass

from enum import Enum
class wiseFilters(Enum):
    w1 = 1
    w2 = 2
    w3 = 3
    w4 = 4

from .survey import Survey
class WISE(Survey):
    def __init__(self,filter=wiseFilters.w1):
        super().__init__()
        self.filter = filter
        self.metadata_root = 'p3am_cdd'
        self.url_root = f"https://irsa.ipac.caltech.edu/ibe/data/wise/allwise/{self.metadata_root}"

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
            width      = edge,
            height     = edge,
            intersect  = 'COVERS'
        )

        #if len(metadata)==0:
        #    # TODO: handle later
        #    #metadata = wise.query_region(
        #    #    coordinate = position,
        #    #    mission    = 'wise',
        #    #    dataset    = 'allwise',
        #    #    table      = self.metadata_root,
        #    #    columns    = 'band,coadd_id',
        #    #    width      = edge,
        #    #    height     = edge,
        #    #    intersect  = 'OVERLAPS'
        #    #)
        #    #status.append(f"> WISE: Position ({position.ra}, {position.dec}) has {len(metadata)} overlapping tiles.")
        #    ##status.append(f"{metadata[metadata['band']==1]}")
        #    status.append(f"> Position ({position.ra}, {position.dec}) has overlapping tiles.")
        #else:
        #    status.append("> ...")
        if len(metadata)==0:
            self.print(f"Position ({position.ra}, {position.dec}) has overlapping tiles.")

        coadd_ids = self.__get_coadd_ids(metadata)
        fits_urls = self.__get_fits_urls(coadd_ids)

        return fits_urls

    #            * * * D E P R E C A T E D * * *
    #def get_cutout(self,position, size):
    #    status = list()
    #
    #    wise = IbeClass()
    #
    #    edge = size.to(u.deg)
    #    metadata = wise.query_region(
    #        coordinate = position,
    #        mission    = 'wise',
    #        dataset    = 'allwise',
    #        table      = self.metadata_root,
    #        columns    = 'band,coadd_id',
    #        width      = edge,
    #        height     = edge,
    #        intersect  = 'COVERS'
    #    )
    #
    #    if len(metadata)==0:
    #        # TODO: handle later
    #        #metadata = wise.query_region(
    #        #    coordinate = position,
    #        #    mission    = 'wise',
    #        #    dataset    = 'allwise',
    #        #    table      = self.metadata_root,
    #        #    columns    = 'band,coadd_id',
    #        #    width      = edge,
    #        #    height     = edge,
    #        #    intersect  = 'OVERLAPS'
    #        #)
    #        #status.append(f"> WISE: Position ({position.ra}, {position.dec}) has {len(metadata)} overlapping tiles.")
    #        ##status.append(f"{metadata[metadata['band']==1]}")
    #        status.append(f"> WISE: Position ({position.ra}, {position.dec}) has overlapping tiles.")
    #    else:
    #        status.append("> WISE: ...")
    #
    #    coadd_ids = self.__get_coadd_ids(metadata)
    #    fits_urls = self.__get_fits_urls(coadd_ids)
    #
    #    # temp prints
    #    if len(fits_urls)==0:
    #        status.append("> https://WHOOPS!")
    #    else:
    #        for fits_url in fits_urls:
    #            status.append(f"> {fits_url}")
    #    print("{0}".format('\n'.join(status)))
    #
    #    return self.get_fits(fits_urls[0]) if len(fits_urls) > 0 else None


