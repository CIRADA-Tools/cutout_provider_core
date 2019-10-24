from astropy import units as u

from .survey_abc import SurveyABC
class NVSS(SurveyABC):
    def __init__(self):
        super().__init__()


    @staticmethod
    def get_supported_filters():
        return None


    def get_filter_setting(self):
        return None


    def get_tile_urls(self,position,size):
        url = 'https://www.cv.nrao.edu/cgi-bin/postage.pl'

        # images can't be bigger than this
        #max_image_pixels = 510 * u.pix
        max_image_pixels = 512 * u.pix
        # this is as low a resolution (higher scales being lower resolution) as you can get
        # without losing any info
        max_pix_scale = 15 * (u.arcsec/u.pix)

        # this could
        speculative_pix_scale = (size / max_image_pixels).to(u.arcsec / u.pix)
        # TODO: NVSS has a maxium practical size, for the most minimum pix_scale at 15"/pix, so we should mosaick 
        #       if speculative_pix_scale > max_pix_scale... if we want to keep the highest resolution.
        pix_scale = max(max_pix_scale, speculative_pix_scale)

        position_components = position.to_string('hmsdms', sep=' ').split(' ')
        #position_components = position.to_string('hmsdms', precision=2, sep=' ').split(' ')

        ra = " ".join(position_components[0: 3])
        dec = " ".join(position_components[3: 6])

        # what data to send in POST request
        post_values = {
            'Equinox': 'J2000',
            'PolType': 'I',
            'RA': ra,
            'Dec': dec,
            'Size': '{0:.5f} {0:.5f}'.format(size.to(u.deg).value),
            'Cells': '{0:.5f} {0:.5f}'.format(pix_scale.value),  # pixel size in arc seconds
            'MAPROJ': 'SIN',
            'Type': 'application/octet-stream'
        }

        self.print(f"({position.ra.to(u.deg).value},{position.dec.to(u.deg).value}) => ({post_values['RA']},{post_values['Dec']})")
        #self.print(f"URL: {self.pack(url, post_values)}")

        return [self.pack(url, post_values)]


    def get_fits_header_updates(self,header,position,size):
        return None

