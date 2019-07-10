import os
import re
import urllib

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from astropy.table import Table
from astropy import units as u
from astropy.units import Quantity

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
            # TODO (Issue #6): shift so ra_i is a centered:i.e., i=0 => (-2.5,2.5) instead of (0,5) for the zone 13 case.
            #ra_i = i * d_ra - d_ra/2.0
            ra_i = i * d_ra
            if ra_i <= ra and ra < ra_i+d_ra:
                return i+self.proj_cells[self.zones.index(zone)]
        return None

    def skycell(self,ra,dec):
        projcell = self.projcell(ra,dec)
        if projcell is None:
            return None
        zone_i = self.zones.index(self.zone(dec))
        projcell_0 = self.proj_cells[zone_i]
        ra_cells = self.row_cells[zone_i]
        ra = self.__sanitize_ra(ra)
        #ra_offset=ra-(projcell-projcell_0)*360.0*u.deg/ra_cells
        ra_offset=ra-(projcell-projcell_0+0.5)*360.0*u.deg/ra_cells
        print(f"zone_i: {zone_i}, projcell: {projcell}, projcell_0: {projcell_0}, ra_cells: {ra_cells}, ra: {ra}, ra_offset: {ra_offset}")
        return projcell

from .survey_abc import SurveyABC
from .survey_filters import grizy_filters
class PanSTARRS(SurveyABC):
    def __init__(self,filter=grizy_filters.i):
        super().__init__()

        self.filter = filter


    # TODO: Determine if useful.
    #   * * * D E P R E C A T E D * * *
    #def __combine_bands(self,bands):
    #    hdus = [b[0] for b in bands]
    #    wcs = WCS(hdus[0].header)
    #    header = wcs.to_header()
    #
    #    data = np.stack(reversed([h.data for h in hdus]))
    #
    #    img = fits.PrimaryHDU(data, header=header)
    #
    #    return img

    # The following two helper functions were taken essentially in full from:
    # https://ps1images.stsci.edu/ps1image.html

    def __getimages(self,ra, dec, size=240, filters="grizy"):
        """Query ps1filenames.py service to get a list of images

           ra, dec = position in degrees
           size = image size in pixels (0.25 arcsec/pixel)
           filters = string with filters to include
           Returns a table with the results
        """

        service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
        #url = (f"{service}?ra={ra}&dec={dec}&size={size}&format=fits&filters={filters}")
        url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
               "&filters={filters}").format(**locals())

        print(f"URL: {url}")

        # TODO: notes...
        # [1] https://outerspace.stsci.edu/display/PANSTARRS/PS1+Image+Cutout+Service
        # [2] https://outerspace.stsci.edu/display/PANSTARRS/PS1+Sky+tessellation+patterns#PS1Skytessellationpatterns-Skycells
        table = Table.read(url, format='ascii')
        print(f"TALBE: {table}")
        return table


    def __geturl(self,ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False):
        """Get URL for images in the table

        ra, dec = position in degrees
        size = extracted image size in pixels (0.25 arcsec/pixel)
        output_size = output (display) image size in pixels (default = size).
                      output_size has no effect for fits format images.
        filters = string with filters to include
        format = data format (options are "jpg", "png" or "fits")
        color = if True, creates a color image (only for jpg or png format).
                Default is return a list of URLs for single-filter grayscale images.
        Returns a string with the URL
        """

        if color and format == "fits":
            raise ValueError("color images are available only for jpg or png formats")
        if format not in ("jpg", "png", "fits"):
            raise ValueError("format must be one of jpg, png, fits")
        table = self.__getimages(ra, dec, size=size, filters=filters)
        #url = (f"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?ra={ra}&dec={dec}&size={size}&format={format}")
        url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
               "ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())

        if output_size:
            url = url + "&output_size={}".format(output_size)
        # sort filters from red to blue
        flist = ["yzirg".find(x) for x in table['filter']]
        table = table[np.argsort(flist)]
        if color:
            if len(table) > 3:
                # pick 3 filters
                table = table[[0, len(table) // 2, len(table) - 1]]
            for i, param in enumerate(["red", "green", "blue"]):
                url = url + "&{}={}".format(param, table['filename'][i])
        else:
            urlbase = url + "&red="
            url = []
            for filename in table['filename']:
                url.append(urlbase + filename)

        return url


    @staticmethod
    def get_supported_filters():
        return grizy_filters


    def get_filter_setting(self):
        return self.filter


    def get_tile_urls(self,position,size):
        pix_scale = 0.25 * (u.arcsec/u.pix)

        ra     = position.ra.to(u.deg).value
        dec    = position.dec.to(u.deg).value
        pixels = int(np.ceil((size / pix_scale).to(u.pix).value))

        url = self.__geturl(ra=ra, dec=dec, size=pixels, filters=self.filter.name, format='fits')

        return url


    def get_fits_header_updates(self,header,position,size):
        survey = type(self).__name__
        header_updates = {
            'BAND':     (f'{self.filter.name}-band', 'Filter used in observation'),
            'DATE-OBS': (Time(header['MJD-OBS'], format='mjd').isot, 'Obs. date'),
            'STK_TYPE': (header['STK_TYPE'], f'{survey} image stack type'),
            'STK_ID':   (header['STK_ID'],   f'{survey} image stack ID'),
            'SKYCELL':  (header['SKYCELL'],  f'{survey} image sky cell'),
            'TESS_ID':  (header['TESS_ID'],  f'{survey} tesselation')
        }
        return header_updates

