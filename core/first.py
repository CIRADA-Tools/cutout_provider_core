from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
import requests, json

from .survey_abc import SurveyABC
class FIRST(SurveyABC):
    def __init__(self):
        super().__init__()
        self.needs_trimming = False


    @staticmethod
    def get_supported_filters():
        return []

    def get_filter_setting(self):
        return None

    # def find_all_sources(self, position,size):
    #     url = 'https://archive.stsci.edu/vlafirst/search.php'
    #     # position = position.to_string('hmsdms', sep=' ')
    #     post_values = {
    #         'ra': '10 50 07.270', #str(position.ra.to(u.deg).value),
    #         'dec': '+30 40 37.52', #str(position.dec.to(u.deg).value),
    #         'equinox': 'J2000',
    #         'radius': (size/2).to(u.arcmin).value,
    #         'selectedColumnsCsv': 'v_ra,v_dec',
    #         'coordformat': 'dec',
    #         'outputformat': 'JSON',
    #         'action': 'Search'
    #     }
    #     query = self.pack(url, post_values)
    #     if self.http is None:
    #         #response = urllib.request.urlopen(request)
    #         matches = requests.get(query, verify=False, timeout=self.http_read_timeout)
    #     else:
    #         matches = self.http.request('GET',query, timeout=self.http_read_timeout)
    #     return matches

    # def get_tile_urls(self,position,size):
    #     matches = self.find_all_sources(position,size)
    #     print("here?", matches, matches.content)
    #     print("size", size.to(u.arcmin).value)
    #     if matches.status_code==200:
    #         all_urls = []
    #         try:
    #             results = json.loads(matches.content)
    #         except Exception as e:
    #             if "no rows found" in str(matches.content):
    #                 raise Exception(f"No FIRST sources found at position: {position.to_string('decimal')} within radius: {(size/2).to(u.arcmin).value}")
    #             return all_urls
    #         url = 'https://third.ucllnl.org/cgi-bin/firstimage'
    #         for coords in results:
    #             print(coords)
    #             sk_coords = SkyCoord(coords['RA (J2000)'], coords['Dec (J2000)'], unit="deg")
    #             post_values = {
    #                 'RA': sk_coords.to_string('hmsdms', sep=' '),  # note that RA, DEC combined into RA here
    #                 'Equinox': 'J2000',
    #                 'ImageSize': size.to(u.arcmin).value,
    #                 'ImageType': 'FITS Image',
    #                 'Download': 1,
    #                 'FITS': 1
    #             }
    #             all_urls.append(self.pack(url, post_values))
    #         print(all_urls)
    #         return all_urls
    #     else:
    #         return []

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


    def get_fits_header_updates(self,header, all_headers=None):
        header_updates = {
            'BAND':     ('1.4 GHz', 'Frequency of observation'),
            'DATE-OBS': (header['DATE-OBS'], 'Obs. date (yearmonth)'),
            # TODO (Issue #6): same as MJD-OB... depricate.
            #'MJD': (Time(header['DATE-OBS']).mjd, 'Median MJD of obs month (00:00:00 on 15th)'),
            # TODO (Issue #6): this is probably already in the wcs part of the header
            'BUNIT': ('Jy/beam', 'Pixel flux unit'),
            'BMAJ':  (header['BMAJ'], 'Beam major axis [deg]'),
            'BMIN':  (header['BMIN'], 'Beam minor axis [deg]'),
            'BPA':   (header['BPA'],  'Beam position angle'),
             # TODO (Issue #6): 'FNAME' inadequate for mosaics...
            'IMFILE': (header['FIELDNAM'], 'FIRST coadded image')
        }
        return header_updates
