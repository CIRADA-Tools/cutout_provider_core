from astropy import units as u

from .survey import Survey
class NVSS(Survey):

    def get_cutout(self,position, size):

        url = 'https://www.cv.nrao.edu/cgi-bin/postage.pl'

        # images can't be bigger than this
        max_image_pixels = 510 * u.pix
        # this is as low a resolution (higher scales being lower resolution) as you can get
        # without losing any info
        max_pix_scale = 15 * (u.arcsec/u.pix)

        # this could
        speculative_pix_scale = (size / max_image_pixels).to(u.arcsec / u.pix)
        pix_scale = max(max_pix_scale, speculative_pix_scale)

        position_components = position.to_string('hmsdms', sep=' ').split(' ')

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

        return self.get_fits(url, post_values)
