from math import log10, floor
from astropy import units as u
import numpy as np

from .survey_abc import SurveyABC
class NVSS(SurveyABC):
    def __init__(self):
        super().__init__()
        self.needs_trimming = False

    @staticmethod
    def get_supported_filters():
        return []

    def get_filter_setting(self):
        return None

    def add_cutout_service_comment(self, hdu):
        hdu.header.add_comment('This cutout was provided by the CIRADA project '\
                                '(www.cirada.ca) using the Postage Stamp Server '\
                                'hosted at NRAO: (https://www.cv.nrao.edu/nvss/postage.shtml) \
                                ', after=-1)

    def get_tile_urls(self,position,size):
        def round_sig(x, sig=3):
            """Rounds x to nearest sigificant figure, sig > 2."""
            return round(x.value,(sig if sig > 2 else 3)-int(floor(log10(abs(x.value))))-1)*x.unit
        # base url
        url = 'https://www.cv.nrao.edu/cgi-bin/postage.pl'
        # images can't be bigger than this
        max_image_pixels = 510 * u.pix # max is 512x512 but down two to be safe
        resolution = 2.25*u.arcsec/u.pix
        # # min allow pixel scale use this as base resolution
        min_scale = 2.25*u.arcsec/u.pix #0.001*u.arcsec/u.pix
        # # max pixel scale without information loss as multiple of min ~15
        max_scale = 15*u.arcsec/u.pix
        # # calculate/desire pixel scale
        desired_scale = round_sig((size/max_image_pixels).to(u.arcsec/u.pix),4)
        # # TODO (Issue #21): NVSS has a maxium practical size, for the most minimum pixel_scale at 15"/pix, so we should mosaick
        # #       if desired_pixel_scale > max_pixel_scale... if we want to keep the highest resolution.
        scale = min(max_scale,max(desired_scale,min_scale))
        pixel_scale = np.ceil(scale/resolution)*resolution # as int multiple of resolution

        # convert position to url query string
        position_components = position.to_string('hmsdms', sep=' ').split(' ')
        #position_components = position.to_string('hmsdms', precision=2, sep=' ').split(' ')
        ra  = " ".join(position_components[0: 3])
        dec = " ".join(position_components[3: 6])

        # what data to send in POST request
        post_values = {
            'Equinox': 'J2000',
            'PolType': 'I',
            'RA': ra,
            'Dec': dec,
            'Size': '{0:.5f} {0:.5f}'.format(size.to(u.deg).value),
            'Cells': '{0:.5f} {0:.5f}'.format(pixel_scale.value),  # pixel size in arc seconds
            'MAPROJ': 'SIN',
            'Type': 'application/octet-stream'
        }
        self.print(f"({position.ra.to(u.deg).value},{position.dec.to(u.deg).value}) => ({post_values['RA']},{post_values['Dec']})")
        #self.print(f"URL: {self.pack(url, post_values)}")
        return [self.pack(url, post_values)]

    def get_fits_header_updates(self,header, all_headers=None):
        return None
