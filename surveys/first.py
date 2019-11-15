from astropy import units as u
from astropy.time import Time

from .survey_abc import SurveyABC
class FIRST(SurveyABC):
    def __init__(self):
        super().__init__()


    @staticmethod
    def get_supported_filters():
        return None


    def get_filter_setting(self):
        return None


    def get_tile_urls(self,position,size):
        url = 'https://third.ucllnl.org/cgi-bin/firstimage'

        position = position.to_string('hmsdms', sep=' ')

        post_values = {
            'RA': position,  # note that this includes DEC
            'Equinox': 'J2000',
            'ImageSize': size.to(u.arcmin).value,
            'ImageType': 'FITS Image',
            'Download': 1,
            'FITS': 1
        }
        # TODO (Issue #9): FIRST mosaicking work in progress... consider using SIAP...
        # MAST I/F notes:
        # [1] https://third.ucllnl.org/cgi-bin/firstcutout -- FIRST bread crumbs
        #     http://sundog.stsci.edu/first/images.html
        #     http://archive.stsci.edu/vlafirst/search.php
        #     http://archive.stsci.edu/searches.html
        #     http://archive.stsci.edu/vo/mast_services.html (SCS vs SIAP)
        # [2] https://astroquery.readthedocs.io/en/latest/mast/mast.html -- Python MAST I/F (appears to have not support for VLA-FIRST)
        # [3] http://docs.astropy.org/en/stable/io/votable/ -- python VOTable I/F
        # [4] http://docs.g-vo.org/DaCHS/tutorial.pdf -- ?
        #str = " ***\n"
        #str += "%s\n" % url
        #str += "post_values = %s" % (" {\n>"+",\n> ".join([f"{k}: {post_values[k]}" for k in post_values.keys()])+"\n> }\n")
        #str += f" ======> {self.pack(url, post_values)}"
        #print(str)

        return [self.pack(url, post_values)]


    def get_fits_header_updates(self,header):
        header_updates = {
            'BAND':     ('1.4 GHz', 'Frequency of observation'),
            'DATE-OBS': (header['DATE-OBS'], 'Obs. date (yearmonth)'),
            # TODO (Issue #6): same as MJD-OB... depricate.
            #'MJD': (Time(header['DATE-OBS']).mjd, 'Median MJD of obs month (00:00:00 on 15th)'),
            # TODO (Issue #6): this is probably already in the wcs part of the header
            'BUNIT': ('Jy/beam', 'Pixel flux unit'),
            'BMAJ':  (header['BMAJ'], 'Beam major axis [deg]'),
            'BMIN':  (header['BMIN'], 'Beam minor axis [deg]'),
            'BPA':   (header['BPA'],  'Beam position angle'),
             # TODO (Issue #6): 'FNAME' inadequate for mosaics...
            'IMFILE': (header['FIELDNAM'], 'FIRST coadded image')
        }
        return header_updates
