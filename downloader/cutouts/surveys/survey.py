import io
import os
import sys
import tempfile
import shutil

import montage_wrapper

import urllib.request
import urllib.parse
import urllib.error

from astropy.io import fits

# abstract class for a survey
from abc import ABC, abstractmethod
class Survey(ABC):
    def __init__(self):
        super().__init__()


    # get data over http post
    def __send_request(self, url, payload=None):

        potential_retries = 5

        if payload:
            data = urllib.parse.urlencode(payload).encode('utf-8')
            request = urllib.request.Request(url, data)
        else:
            request = url

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

        print(f"WARNING: Bailed on fetch '{url}'")


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
        self.make_dir(input_dir)

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


    # tries to create a fits file from bytes.
    # on fail it returns None
    def create_fits(self, data):
        # a pretend file in memory
        fits_file = io.BytesIO(data)
        fits_file.seek(0)

        try:
            f = fits.open(fits_file)
        except OSError as e:
            print("Badly formatted FITS file: {0}\n\treturning None".format(str(e)), file=sys.stderr)
            return None

        head = f[0].header
        if head['NAXIS'] == 0 or head['NAXIS1'] == 0 or head['NAXIS2'] == 0:
            return None

        data = f[0].data

        # if data is all 0
        if data.min() == 0 and data.max() == 0:
            return None

        return f


    def get_fits(self, url, payload=None):
        response = self.__send_request(url, payload)
        # note that it returns None if the response isn't a valid fits
        return self.create_fits(response)


    # grab a cutout of size <size> centered on <position>
    @abstractmethod
    def get_cutout(self, position, size):
        pass
