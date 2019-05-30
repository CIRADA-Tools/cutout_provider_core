# downloading a FITS file from a URL

import io

import sys

import urllib.request
import urllib.parse
import urllib.error

from astropy.io import fits


# get data over http post
def send_request(url, payload=None):

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
                # Note: this apperently cause an FITSFixedWarning, e.g.,
                #
                #    WARNING: FITSFixedWarning: 'datfix' made the change 'Changed '' to '2013-05-16T18:02:48.842''. [astropy.wcs.wcs]
                #
                # in astropy.wcs (cf., http://docs.astropy.org/en/stable/wcs/); hence this try..except wrapping.
                #
                pass

        potential_retries -= 1

    print(f"WARNING: Bailed on fetch '{url}'")


# tries to create a fits file from bytes.
# on fail it returns None
def create_fits(data):

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


def get_fits(url, payload=None):

    response = send_request(url, payload)
    # note that it returns None if the response isn't a valid fits
    return create_fits(response)

