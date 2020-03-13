import os
import re
import urllib

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from astropy import units as u
from astropy.table import Table
from astropy.table import vstack
from astropy.units import Quantity
from astropy.coordinates import SkyCoord

from .survey_abc import SurveyABC
from .survey_filters import grizy_filters
from .toolbox import pad_string_lines

class PS1SkyTessellationPatterns:
    # cf., https://outerspace.stsci.edu/display/PANSTARRS/PS1+Sky+tessellation+patterns
    def __init__(self):
        # load the projection cell table
        this_source_file_dir = re.sub(r"(.*/).*$",r"\1",os.path.realpath(__file__))
        ps1grid = fits.open(this_source_file_dir+'ps1grid.fits')[1].data
        self.zones = list(ps1grid['ZONE'])
        self.proj_cells = list(ps1grid['PROJCELL'])
        self.row_cells = list(ps1grid['NBAND'])
        self.decs  = [d*u.deg for d in ps1grid['DEC']]
        self.x_sizes = [x*u.pixel for x in ps1grid['XCELL']]
        self.y_sizes = [y*u.pixel for y in ps1grid['YCELL']]
        self.x_subs  = [x for x in ps1grid['CRPIX1']]
        self.y_subs  = [y for y in ps1grid['CRPIX2']]
        self.min_decs = [d*u.deg for d in ps1grid['DEC_MIN']]
        self.max_decs = [d*u.deg for d in ps1grid['DEC_MAX']]

        self.pix_scale = 0.25 * (u.arcsec/u.pix)
        self.sub_cells = 10

    def __sanitize_ra(self,ra):
        if not isinstance(ra,Quantity):
            ra = (ra % 360.0) * u.deg
        return ra

    def __sanitize_dec(self,dec):
        if not isinstance(dec,Quantity):
            dec = dec * u.deg
        return dec

    def dec(self, zone):
        if not zone in self.zones:
            return None
        return self.decs[self.zones.index(zone)]

    def zone(self,dec):
        dec = self.__sanitize_dec(dec)
        for i,zone in enumerate(self.zones):
            if self.min_decs[i] <= dec and dec < self.max_decs[i]:
                 return zone
        return None

    def projcell(self,ra,dec):
        zone = self.zone(dec)
        if zone is None:
            return None
        ra = self.__sanitize_ra(ra)
        ra_cells = self.row_cells[self.zones.index(zone)]
        d_ra = 360.0*u.deg/ra_cells
        for i in range(ra_cells):
            ra_i = ((i * d_ra - d_ra/2.0).value % 360.0) * u.deg
            if (i == 0 and ((ra_i <= ra and ra < 0*u.deg) or (0*u.deg <= ra and ra < d_ra/2.0))) or (ra_i <= ra and ra < ra_i+d_ra):
                return i+self.proj_cells[self.zones.index(zone)]
        return None

    def projcell_center(self,ra,dec):
        projcell = self.projcell(ra,dec)
        if projcell is None:
            return None
        zone_i = self.zones.index(self.zone(dec))
        dec_center = self.decs[zone_i]
        ra_center = (projcell-self.proj_cells[zone_i])*360.0*u.deg/self.row_cells[zone_i]
        return SkyCoord(ra_center,dec_center)

    # debug
    def __get_tile_hieght(self,dec):
        zone_i = self.zones.index(self.zone(dec))
        return self.max_decs[zone_i]-self.min_decs[zone_i]
    def __gvstacket_tangent_plane_dec(self,ra,dec):
        return np.arctan(np.tan((dec+self.__get_tile_hieght(dec)/2.0).to(u.rad).value)*np.cos(ra.to(u.rad).value))*180.0*u.deg/np.pi

    def skycell(self,ra,dec):
        projcell = self.projcell(ra,dec)
        if projcell is None:
            return None
        #r = SkyCoord(self.__sanitize_ra(ra),self.__sanitize_dec(dec))
        #r_0 = self.projcell_center(ra,dec)
        ra = self.__sanitize_ra(ra)
        dec = self.__sanitize_dec(dec)
        zone_i = self.zones.index(self.zone(dec))
        d_dec = (self.max_decs[zone_i]-self.min_decs[zone_i])/(self.sub_cells-1)
        min_dec = self.min_decs[zone_i]
        si_dec = None
        # TODO: cat seem to get xy to make table query
        for i in range(self.sub_cells):
            dec_i = i * d_dec + min_dec
            #print(f"{dec}: ({dec_i},{dec_i+d_dec})")
            if dec_i <= dec and dec < dec_i+d_dec:
                 si_dec = i
                 break
        print(f"si_dec={si_dec}")
        return projcell


class PANSTARRS(SurveyABC):
    def __init__(self,filter=grizy_filters.i):
        super().__init__()
        self.pixel_scale = 0.25 * (u.arcsec/u.pix)
        self.filter = filter
        self.needs_trimming = True # may not need to trim??

    @staticmethod
    def get_supported_filters():
        return grizy_filters

    def add_cutout_service_comment(self, hdu):
        pass
        # hdu.header.add_comment(pad_string_lines('The PS1 Image Cutout Service ' \
        #             'hosted by STScI was used to provide this cutout: ' \
        #             '(https://outerspace.stsci.edu/display/PANSTARRS/PS1+Image+Cutout+Service) \
        #             '), after=-1)

    def get_skycells(self,position,size):
        service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
        def make_url(a,d):
            return f"{service}?ra={a.to(u.deg).value}&dec={d.to(u.deg).value}&format=fits&filters={self.filter.name}"

        # extract ra/dec's
        ra = position.ra
        dec = position.dec

        # the skycell at input (ra,dec)
        url = make_url(ra,dec)
        skycells = Table.read(url, format='ascii', fast_reader=True)

        # TODO (Issue #8): kludge: this is a really bad way of finding the neigbhoring skeyscells...
        # Notes: https://outerspace.stsci.edu/display/PANSTARRS/PS1+Sky+tessellation+patterns
        urls = list()
        r = size/2.0
        urls.append(make_url(ra-r,dec-r))
        urls.append(make_url(ra-r,dec+r))
        urls.append(make_url(ra+r,dec+r))
        urls.append(make_url(ra+r,dec-r))
        urls.append(make_url(ra-r,dec))
        urls.append(make_url(ra+r,dec))
        urls.append(make_url(ra,dec-r))
        urls.append(make_url(ra,dec+r))
        # try:
        for url in urls:
            try:
                self.print("GETTING SKYCELL AT: ", url)
                sc = Table.read(url, format='ascii')
                is_in_skycells = False
                for skycell in skycells:
                    if sc['projcell'] == skycell['projcell'] and sc['subcell'] == skycell['subcell']:
                        is_in_skycells = True
                if not is_in_skycells:
                    skycells = vstack([skycells,sc])
            except Exception as e:
                self.print("prob doesn't exist?", str(e))
        return skycells

    def get_filter_setting(self):
        return self.filter

    def get_tile_urls(self,position,size):
        urls = list()
        size_pixels = int(np.ceil((size/self.pixel_scale).to(u.pix).value))
        for skycell in self.get_skycells(position,size):
            urls.append(f"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?ra={position.ra.to(u.deg).value}&dec={position.dec.to(u.deg).value}&size={size_pixels}&format=fits&red={skycell['filename']}")
        return urls

    def get_fits_header_updates(self,header, all_headers=None):
        survey = type(self).__name__
        header_updates = {
            'BAND':     (f'{self.filter.name}-band', 'Filter used in observation'),
            'DATE-OBS': (Time(header['MJD-OBS'], format='mjd').isot, 'Obs. date'),
            'STK_TYPE': (header['STK_TYPE'], f'{survey} image stack type'),
            'STK_ID':   (header['STK_ID'],   f'{survey} image stack ID'),
            'SKYCELL':  (header['SKYCELL'],  f'{survey} image sky cell'),
            'TESS_ID':  (header['TESS_ID'],  f'{survey} tesselation')
        }
        #TODO consider adding in records from all images gone into mosaic
        # using all_headers
        if all_headers:
            for num,head in enumerate(all_headers):
                fpartkeys = ['TESS_ID', 'SKYCELL','STK_ID']
                nameparts = [head[key] for key in fpartkeys]
                header_updates['IMFILE'+str(num+1).zfill(2)]='.'.join(nameparts) + '.stk.i.unconv.fits_sci.fits'
                ###create single string - FILNAM12 goes after a constant
        return header_updates
