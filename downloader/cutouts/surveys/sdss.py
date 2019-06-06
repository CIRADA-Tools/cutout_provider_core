from astropy import units as u
import urllib

from enum import Enum
class sdssFilters(Enum):
    g = 1 
    r = 2
    i = 3
    z = 4
    y = 5

from .survey import Survey
class SDSS(Survey):
    def __init__(self,filter=sdssFilters.g):
        super().__init__()
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


    #            * * * D E P R E C A T E D * * *
    #def get_cutout(self,position,size):
    #
    #    pix_scale = 0.262 * (u.arcsec/u.pix)
    #    pixels = (size / pix_scale).to(u.pix)
    #
    #    query_dict = {
    #        'ra': position.ra.value,
    #        'dec': position.dec.value,
    #        'layer': 'sdss2',
    #        'bands': 'grz',
    #        'pixscale': pix_scale.value,
    #        'size': int(pixels.value),
    #    }
    #
    #    query_string = urllib.parse.urlencode(query_dict)
    #
    #    url = "http://legacysurvey.org/viewer/fits-cutout?{query}".format(query=query_string)
    #
    #    return self.get_fits(url)
