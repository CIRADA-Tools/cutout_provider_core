import io
import os
import sys
import tempfile
import shutil

import urllib.request
import urllib.parse
import urllib.error

import re

import numpy as np
from astropy.wcs import WCS
from astropy.io import fits

from astropy import units as u

import montage_wrapper
from astropy.nddata.utils import Cutout2D

# abstract class for a survey
from abc import ABC, abstractmethod
class SurveyABC(ABC):
    def __init__(self,
        pid = None        # 4 threads
    ):
        super().__init__()

        self.pid = None


    def set_pid(self,pid):
        self.pid = pid


    def get_sexadecimal_string(self,position):
        sexadecimal = "%02d%02d%02.0f" % position.ra.hms+re.sub(r"([+-])\d",r"\1","%+d%02d%02d%02.0f" % position.dec.signed_dms)
        return sexadecimal 


    def get_ra_dec_string(self,position):
        return "%f%+f degree" % (position.ra.to(u.deg).value,position.dec.to(u.deg).value)


    def print(self,msg,show_caller=False):
        my_name    = type(self).__name__ + (f"[{sys._getframe(1).f_code.co_name}]" if show_caller else "")
        my_pid     = "" if self.pid is None else f"pid={self.pid}"
        my_filter  = (lambda f: "" if f is None else f"filter='{f.name}'")(self.get_filter_setting())
        prefix = f"{my_name}({my_pid}{'' if my_pid=='' or my_filter=='' else ','}{my_filter})"
        prefixed_output = "\n".join([f"{prefix}: {s}" for s in msg.splitlines()])
        print(prefixed_output)


    def pack(self,url,payload=None):
        if payload:
            data = urllib.parse.urlencode(payload).encode('utf-8')
            request = urllib.request.Request(url, data)
        else:
            request = url
        return request


    # get data over http post
    def __send_request(self, url, payload=None):

        potential_retries = 5

        request = self.pack(url,payload)

        while potential_retries > 0:

            try:
                response = urllib.request.urlopen(request)
            except urllib.error.HTTPError:
                pass
            except ConnectionResetError:
                pass
            else:
                try:
                    response_data = bytearray(response.read())
                    return response_data
                except:
                    # E.g., One noted error was 'IncompleteRead'...
                    #
                    # Note: this apperently caused an FITSFixedWarning, e.g.,
                    #
                    #    WARNING: FITSFixedWarning: 'datfix' made the change 'Changed '' to '2013-05-16T18:02:48.842''. [astropy.wcs.wcs]
                    #
                    # in astropy.wcs (cf., http://docs.astropy.org/en/stable/wcs/); hence this try..except wrapping.
                    #
                    pass

            potential_retries -= 1

        self.print(f"WARNING: Bailed on fetch '{url}'")


    # make the directory structure if it doesn't exist
    def __make_dir(self, dirname):
        try:
            os.makedirs(dirname)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise


    def mosaic(self, cutouts):
        td = tempfile.mkdtemp()
        input_dir = '{directory}/input'.format(directory=td)
        output_dir = '{directory}/output'.format(directory=td)
        self.__make_dir(input_dir)

        try:
            for i, c in enumerate(cutouts):
                with open('{directory}/{name}.fits'.format(directory=input_dir, name=i), 'wb') as tmp:
                    tmp.write(bytes(c))

            os.listdir(td)
            montage_wrapper.mosaic(input_dir, output_dir)

            with open('{outdir}/mosaic.fits'.format(outdir=output_dir), 'rb') as f:
                merged = f.read()
        finally:
            shutil.rmtree(output_dir)
            shutil.rmtree(input_dir)
            shutil.rmtree(td)

        return merged


    def get_image(self,hdu):
        img_data = np.squeeze(hdu[0].data)
        img = fits.PrimaryHDU(img_data, header=hdu[0].header)

        # writing to a pretend file in memory
        mem_file = io.BytesIO()
        img.writeto(mem_file)
        return mem_file.getvalue()


    # tries to create a fits file from bytes.
    # on fail it returns None
    def create_fits(self, data, position):
        # a pretend file in memory
        fits_file = io.BytesIO(data)
        fits_file.seek(0)

        # open fits hdul
        try:
            hdul = fits.open(fits_file)
        except OSError as e:
            self.print("Badly formatted FITS file: {0}\n\treturning None".format(str(e)), file=sys.stderr)
            return None

        # get/check header field
        header = hdul[0].header
        if ('NAXIS'  in header and header['NAXIS']  <  2) or \
           ('NAXIS1' in header and header['NAXIS1'] == 0) or \
           ('NAXIS2' in header and header['NAXIS2'] == 0):
            self.print("WARINING: Ill-defined 'NAXIS/i': skipping...")
            return None

        # get/check data field
        data = hdul[0].data
        if data.min() == 0 and data.max() == 0:
            self.print("WARNING: Fits file contains no data: skipping...")
            return None

        #
        # OK, let's make sure there's some standard stuff in our headers!
        # We'll fix up the WCS later...
        # Originator: Yjan Gordon (JN1119).
        #

        # remove superfluous nesting
        # TODO: This seems to have no effect, once the fits file is reload... investigate!
        if data.ndim > 2:
             data = data[0]

        # define the NAXISi=1,2 unit fields
        if ('CUNIT1' in header) and ('CUNIT2' in header) and (header['CUNIT1'] != header['CUNIT2']):
            self.print(f"WARNING: ('CUNIT1'={header['CUNIT1']}) != ('CUNIT2'={header['CUNIT2']}) !")
        elif not (('CUNIT1' in header) or ('CUNIT2' in header)):
             pos_units = 'deg'
             if 'CTYPE1' in header:
                 header.insert('CTYPE1', ('CUNIT1', pos_units), after=True)
             else:
                 self.print(f"WARNING: 'CTYPE1' undefined!")
                 header['CUNIT1'] = pos_units
             if 'CTYPE2' in header:
                 header.insert('CTYPE2', ('CUNIT2', pos_units), after=True)
             else:
                 self.print(f"WARNING: 'CTYPE2' undefined!")
                 header['CUNIT2'] = pos_units

        # set LATPOLE
        if not ('LATPOLE' in header):
            header['LATPOLE'] = (position.dec.to(u.deg).value, 'Native latitude of celestial pole')

        # set LONPOLE
        if not ('LONPOLE' in header):
            ra = position.ra.to(u.deg).value
            pole_longitude = 180.0 if 90.0 < ra and ra <= 270.0 else 0.0
            header['LONPOLE'] = (pole_longitude, 'Native longitude of celestial pole')

        # define RADESYS field
        if not ('RADESYS' in header):
            header['RADESYS'] = ('FK5', 'assumed')

        # set SURVEY
        if not ('SURVEY' in header):
            header['SURVEY'] = (type(self).__name__, 'Survey image obtained from')

        # set BAND default
        if not ('BAND' in header):
            header['BAND'] = 'na'

        # set EPOCH
        if not ('EPOCH' in header):
            header['EPOCH'] = (2000.0, 'Julian epoch of observation')

        # define DATE-OBS field
        if not ('DATE-OBS' in header):
            header['DATE-OBS'] = 'na'

        return hdul


    def get_fits(self, url, position, payload=None):
        self.print(f"Fetching: {url}")
        response = self.__send_request(url, payload)
        # note that it returns None if the response isn't a valid fits
        return self.create_fits(response, position)


    def get_tiles(self,position,size):
        urls = self.get_tile_urls(position,size)
        if len(urls) > 0:
            hdu_lists = [h for h in [self.get_fits(url,position) for url in urls] if h]
        else:
            return None
        return hdu_lists


    def paste_tiles(self,hdu_tiles,position,size):
        if hdu_tiles is None:
            return None

        hdu = None
        if len(hdu_tiles) > 1:
            self.print(f"Pasting {len(hdu_tiles)} at J{self.get_sexadecimal_string(position)}")
            try:
                imgs = [img for img in [self.get_image(tile) for tile in hdu_tiles]]
                hdu  = self.create_fits(self.mosaic(imgs),position)
            except montage_wrapper.status.MontageError as e:
                self.print(f"Mosaicing Failed: {e}: file={sys.stderr}",True)
        elif len(hdu_tiles) == 1:
                hdu = hdu_tiles[0]
        return hdu


    def trim_tile(self, hdu, position, size):
        if hdu is None:
            return None
    
        w = WCS(hdu[0].header)
    
        # trim to 2d from nd
        naxis = w.naxis
        while naxis > 2:
            w = w.dropaxis(2)
            naxis -= 1
        img_data = np.squeeze(hdu[0].data)
    
        stamp = Cutout2D(img_data, position, size, wcs=w, mode='trim', copy=True)
        img = fits.PrimaryHDU(stamp.data, header=stamp.wcs.to_header())
    
        # writing to a pretend file in memory
        mem_file = io.BytesIO()
        img.writeto(mem_file)
        return img


    def get_cutout(self,position, size):
        try:
            tiles  = self.get_tiles(position,size)
            #tile   = self.paste_tiles(tiles,position,size)
            #cutout = self.trim_tile(tile,position,size)
            cutout = tiles[0] if tiles else None
        except Exception as e:
            self.print(f"{e}",True)
            cutout = None

        self.print(f"Finished processing J{self.get_sexadecimal_string(position)} ({self.get_ra_dec_string(position)}) cutout of size {size}.")

        return cutout


    @staticmethod
    @abstractmethod
    def get_supported_filters():
        pass


    @abstractmethod
    def get_filter_setting(self):
        pass


    @abstractmethod
    def get_tile_urls(self,position,size):
        pass


    #@abstrachmethod
    #def format_fits_header(self,hdu,position):
    #    pass

