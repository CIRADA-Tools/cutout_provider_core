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

        return [self.pack(url, post_values)]


    def get_fits_header_updates(self,header,position,size):
        header_updates = {
            'BAND':     ('1.4 GHz', 'Frequency of observation'),
            'DATE-OBS': (header['DATE-OBS'], 'Obs. date (yearmonth)'),
            # TODO: same as MJD-OB... depricate.
            #'MJD': (Time(header['DATE-OBS']).mjd, 'Median MJD of obs month (00:00:00 on 15th)'),
            # TODO: this is probably already in the wcs part of the header
            'BUNIT': ('Jy/beam', 'Pixel flux unit'),
            'BMAJ':  (header['BMAJ'], 'Beam major axis [deg]'),
            'BMIN':  (header['BMIN'], 'Beam minor axis [deg]'),
            'BPA':   (header['BPA'],  'Beam position angle'),
             # TODO: 'FNAME' inadequate for mosaics...
            'IMFILE': (header['FIELDNAM'], 'FIRST coadded image')
        }
        return header_updates

