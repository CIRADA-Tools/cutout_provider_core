from astropy import units as u

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


    def format_fits_header(self,hdu,position,size):
        return hdu

