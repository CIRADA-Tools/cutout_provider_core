import urllib, requests, io
from astropy import units as u
from astropy.io.votable import parse#parse_single_table
# from astroquery.sdss import SDSS as astroSDSS

from .survey_abc import SurveyABC
from .survey_filters import gleam_frequency
from .toolbox import pad_string_lines
class GLEAM(SurveyABC):
    def __init__(self,filter=gleam_frequency.f4):
        # GLEAM 'filter' is actually frequency
        super().__init__()
        self.filter = filter
        self.needs_trimming = True
        self.projection = "SIN"
        self.base_url = "http://gleam-vo.icrar.org"

    @staticmethod
    def get_supported_filters():
        return gleam_frequency

    def get_filter_setting(self):
        return self.filter

    def add_cutout_service_comment(self, hdu):
        hdu.header.add_comment(pad_string_lines('This cutout was provided using the GLEAM Postage Stamp Service ' \
                                'through the GLEAM VO vlient (http://gleam-vo.icrar.org/gleam_postage/q/info) \
                                '), after=-1)

    def get_fits_matches(self,position,size):
        deg_size = float(size.to_value(u.degree))
        query_dict = {
            'freq': self.filter.value,
            'pos': f"{position.ra.value},{position.dec.value}",
            'proj_opt': self.projection,
            'size': deg_size,
        }
        query_string = urllib.parse.urlencode(query_dict)
        url = f"{self.base_url}/gleam_postage/q/siap.xml?{query_string}"
        if self.http is None:
            matches = requests.get(url, verify=False, timeout=self.http_read_timeout)
        else:
            matches = self.http.request('GET',url, timeout=self.http_read_timeout)
        return matches

    # code referenced and adapted from https://github.com/ICRAR/gleamvo-client/blob/master/gleam_client.py
    #file_id=mosaic_Week3_162-170MHz.fits&regrid=1&projection=SIN&fits_format=1
    def get_tile_urls(self,position,size):
        #GLEAMCUTOUT?
        matches = self.get_fits_matches(position,size)
        good=False
        try: # different response object for cutout server and CLI
            if matches.status_code==200: # webserver
                good =True
        except:
            if matches.status==200: # CLI
                good =True
        all_urls = []
        if good:
            try:
                if self.http is None:
                    response_data = bytearray(matches.content)
                else:
                    response_data = bytearray(matches.data)
                # a pretend file in memory
                xml = io.BytesIO(response_data)
                xml.seek(0)
                # urls = parse_single_table(xml).array['accref']
                tables = parse(xml)
                for t in tables.iter_tables():
                    url = t.array['accref']
                    if url:
                        url = url[0].decode("utf-8")
                        all_urls.append(url)
            except Exception as e:
                print(str(e))
                if "no rows found" in str(matches.content):
                    raise Exception(f"No GLEAM sources found at position: {position.to_string('decimal')} within radius: {(size/2).to(u.arcmin).value}")
                return all_urls

        return all_urls

    def get_fits_header_updates(self,header, all_headers=None):
        return None
