import urllib
from astropy import units as u
from astroquery.sdss import SDSS as astroSDSS

from .toolbox import get_sexadecimal_string
from .survey_abc import SurveyABC
from .survey_filters import ugriz_filters
class SDSS(SurveyABC):
    def __init__(self,filter=ugriz_filters.g):
        super().__init__()
        self.filter = filter
        self.needs_trimming = True
        self.data_release = 14 #default is 14

    @staticmethod
    def get_supported_filters():
        return ugriz_filters

    def get_filter_setting(self):
        return self.filter

    # def get_tiles(self, position, size):
    #     self.print(f"getting tiles for {str(position)}\n" )
    #     url_placeholder = f"SDSS_J{get_sexadecimal_string(position)}_s{(size/2).to(u.arcmin).value}_arcmin_{self.filter.name}"
    #     try:
    #         hdul_lists = astroSDSS.get_images(coordinates=position, radius = 5.0*u.arcmin, band=self.filter.name, timeout=self.http_read_timeout, data_release=self.data_release)
    #     except Exception as e:
    #         print(str(e))
    #         hdul_lists = []
    #     if not hdul_lists:
    #         raise Exception(f"no valid {type(self).__name__} tiles found for Position: {position.ra.degree}, {position.dec.degree} filter: {self.filter.name}")
    #     print(hdul_lists[0].info())
    #     # print("here?", hdul_lists)
    #     hdul_tups = [(hdul_l[0], f"{url_placeholder}-{num}") for num,hdul_l in enumerate(hdul_lists) if hdul_l[0]]
    #     print(hdul_tups)
    #     return hdul_tups

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

    def get_fits_header_updates(self,header, all_headers=None):
        return None
