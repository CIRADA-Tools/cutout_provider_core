from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS

import urllib

import numpy
from astropy.table import Table

from .survey import Survey
from .fits_request import get_fits


class PanSTARRS(Survey):

    @staticmethod
    def get_cutout(position, size, filters='gri'):

        pix_scale = 0.25 * (u.arcsec/u.pix)
        pixels = (size / pix_scale).to(u.pix)

        ra = position.ra.to(u.deg).value
        dec = position.dec.to(u.deg).value
        pixels = int((size / pix_scale).to(u.pix).value)

        url = PanSTARRS.geturl(ra=ra, dec=dec, size=pixels, filters=filters, format='fits')

        if url:
            bands = [y for y in [get_fits(x) for x in url] if y]
            return PanSTARRS.combine_bands(bands)
        else:
            return None

    @staticmethod
    def combine_bands(bands):

        hdus = [b[0] for b in bands]
        wcs = WCS(hdus[0].header)
        header = wcs.to_header()

        data = numpy.stack(reversed([h.data for h in hdus]))

        img = fits.PrimaryHDU(data, header=header)

        return img

    # The following two helper functions were taken essentially in full from:
    # https://ps1images.stsci.edu/ps1image.html

    @staticmethod
    def getimages(ra, dec, size=240, filters="grizy"):

        """Query ps1filenames.py service to get a list of images

        ra, dec = position in degrees
        size = image size in pixels (0.25 arcsec/pixel)
        filters = string with filters to include
        Returns a table with the results
        """

        service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
        url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
               "&filters={filters}").format(**locals())

        table = Table.read(url, format='ascii')
        return table

    @staticmethod
    def geturl(ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False):

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
        table = PanSTARRS.getimages(ra, dec, size=size, filters=filters)
        url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
               "ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
        if output_size:
            url = url + "&output_size={}".format(output_size)
        # sort filters from red to blue
        flist = ["yzirg".find(x) for x in table['filter']]
        table = table[numpy.argsort(flist)]
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
