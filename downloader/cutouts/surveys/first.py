from astropy import units as u

from .survey import Survey
from .fits_request import get_fits


class FIRST(Survey):

    @staticmethod
    def get_cutout(position, size):

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

        return get_fits(url, post_values)
