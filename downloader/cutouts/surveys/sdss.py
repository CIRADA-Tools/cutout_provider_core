from astropy import units as u
import urllib

from .survey import Survey
class SDSS(Survey):

    def get_cutout(self,position,size):

        pix_scale = 0.262 * (u.arcsec/u.pix)
        pixels = (size / pix_scale).to(u.pix)

        query_dict = {
            'ra': position.ra.value,
            'dec': position.dec.value,
            'layer': 'sdss2',
            'bands': 'grz',
            'pixscale': pix_scale.value,
            'size': int(pixels.value),
        }

        query_string = urllib.parse.urlencode(query_dict)

        url = "http://legacysurvey.org/viewer/fits-cutout?{query}".format(query=query_string)

        return self.get_fits(url)
