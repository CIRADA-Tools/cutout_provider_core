import urllib
from astropy import units as u


from .survey import Survey
from .survey_filters import grizy_filters
class SDSS(Survey):
    def __init__(self,filter=grizy_filters.g,trimming_on=True):
        super().__init__(trimming_on)
        self.filter = filter


    def get_tile_urls(self,position,size):
        pix_scale = 0.262 * (u.arcsec/u.pix)
        pixels = (size / pix_scale).to(u.pix)

        query_dict = {
            'ra': position.ra.value,
            'dec': position.dec.value,
            'layer': 'sdss2',
            'bands': self.filter.name,
            'pixscale': pix_scale.value,
            'size': int(pixels.value),
        }

        query_string = urllib.parse.urlencode(query_dict)

        url = "http://legacysurvey.org/viewer/fits-cutout?{query}".format(query=query_string)

        return [url]

