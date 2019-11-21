import io
import os
import sys
import traceback
import tempfile
import shutil

# *** IO_WRAPPER ***
# TODO (Issue #11): ...
# Notes: https://stackoverflow.com/questions/1218933/can-i-redirect-the-stdout-in-python-into-some-sort-of-string-buffer/33979942#33979942
from io import TextIOWrapper, BytesIO

import urllib.request
import urllib.parse
import urllib.error
import urllib3
import requests
from time import sleep

import re

import numpy as np
from astropy.wcs import WCS
from astropy.time import Time
from astropy.io import fits
from .survey_filters import HeaderFilter
from .survey_filters import get_header_pretty_string
from .survey_filters import sanitize_fits_date_fields

from astropy import units as u

import montage_wrapper as montage
from astropy.nddata.utils import Cutout2D


def get_sexadecimal_string(position):
    sexadecimal = "%02d%02d%02.0f" % position.ra.hms+re.sub(r"([+-])\d",r"\1","%+d%02d%02d%02.0f" % position.dec.signed_dms)
    return sexadecimal

from enum import Enum
class processing_status(Enum):
    idle      = "Waiting for fetching request"
    none      = "Cutout doesn't exist"
    fetching  = "Cutout fetching and processing"
    corrupted = "Header corrupted"
    error     = "Programming error"
    bailed    = "Fetching aborted"
    done      = "Cutout processed"

    @classmethod
    def __get_unreprocessable_list(self):
        return [self.none, self.corrupted, self.error]

    @classmethod
    def __get_unreprocessable_file_names(self,filename):
        prefix   = re.sub(r"\.\w+?$",".",filename)
        suffixes = [f"{s.name}" for s in self.__get_unreprocessable_list()]
        return [prefix+suffix for suffix in suffixes]

    @classmethod
    def touch_file(self, filename, status, msg):
        for unreprocessable in self.__get_unreprocessable_list():
            if status == unreprocessable:
                try:
                    handle = open(re.sub(r"\.fits$",f".{status.name}",filename),'w')
                    handle.write(unreprocessable.value+"...\n"+msg)
                    handle.close()
                except:
                    return False
                break
        return True

    @classmethod
    def is_processed(self,filename):
        files = [filename]+self.__get_unreprocessable_file_names(filename)
        for file in files:
            if os.path.isfile(file):
                return True
        return False

    @classmethod
    def get_file_listing(self,filename):
        file_listing = list()
        for file in [filename]+self.__get_unreprocessable_file_names(filename):
            if os.path.exists(file):
                file_listing.append(file)
        return file_listing

    @classmethod
    def flush(self,filename,is_all=False):
        msgs = list()
        for file in ([filename] if is_all else [])+self.__get_unreprocessable_file_names(filename):
            if os.path.exists(file):
                msgs.append(f"Flushed: {file}")
                os.remove(file)
        return "\n".join(msgs) if len(msgs) > 0 else None


# abstract class for a survey
from abc import ABC, abstractmethod
class SurveyABC(ABC):
    """
                      /\
                     /  \
                    / _o \
                   / <(\  \
                  /   />`A \
                 '----------`
                  Commenting
              Under Construction

        The following is a list of the class methods grouped by function.

            # class instance tracking number
            set_pid(pid)

            # http request settings
            set_http_request_retries(retries)
            set_http_wait_retry_s(wait_seconds)

            # for attatching http buffer pool manager (re., urllib3) -- default: urllib
            attach_http_pool_manager(http_pool_manager)

            # ra-dec output formating string
            get_sexadecimal_string(position)
            get_ra_dec_string(position)

            # std out printing function
            sprint(msg, diagnostic_msg=None, show_caller=False, is_traceback=True, buffer=True)
            print(msg, diagnostic_msg=None, show_caller=False, is_traceback=True, buffer=True)
            set_print_to_stdout()
            unset_print_to_stdout()

            # cutout fetching and processing stack
            pack(url, payload=None)
            standardize_fits_header_DATE_and_DATE_OBS_fields(date_obs_value)
            create_fits(data)
            get_fits(url, payload=None)
            get_tiles(position, size)
            mosaic(cutouts)
            paste_tiles(hdul_tiles, position, size)
            trim_tile(hdu, position, size)
            format_fits_hdu(hdu, position, size)

            # main cutout fetching rounting
            get_cutout(position, size)

            # abstract base class functions required by survey/child classes
            get_supported_filters()
            get_filter_setting(self)
            get_tile_urls(self,position,size)
            get_fits_header_updates(self,hdu,position,size)

        The main routine is get_cutout(position, size) which triggers the whole cutout
            o Fetching
            o Mosaicking
            o Trimming
            o Header Formating
        processing stack. The routine returns the following information.

               --------------------------------------
               | Dict Element | Contents            |
               |______________|_____________________|
               | cutout       | Processed Cutout    |
               | request_urls | Request URLS        |
               | raw_tiles    | Unprocessed Cutouts |
               | message      | Diagnotics Message  |
               | status       | Processing Status   |
               --------------------------------------


                      /\
                     /  \
                    / _o \
                   / <(\  \
                  /   />`A \
                 '----------`
                  Commenting
              Under Construction
    """
    def __init__(self,
        http_pool_manager = None, # thread safe pool manager (i.e., urllib3)
        pid = None # 4 threads
        ):
        ABC.__init__(self)

        self.processing_status = processing_status.idle

        self.pid = None
        self.http = None

        self.message_buffer = ""
        self.print_to_stdout = True

        self.request_urls_stack = list()
        self.mosaic_hdul_tiles_stack = list()

        self.http_request_retries = 5
        self.http_wait_retry_s = 5


    def __pop_processing_status(self):
        status = self.processing_status
        self.processing_status = processing_status.idle
        return status


    def set_pid(self, pid):
        self.pid = pid
        return self


    def set_http_request_retries(self,retries):
        self.http_request_retries = retries
        return self


    def set_http_wait_retry_s(self,wait_seconds):
        self.http_wait_retry_s = wait_seconds
        return self


    def attach_http_pool_manager(self,http_pool_manager):
        self.http = http_pool_manager
        return self


    def get_sexadecimal_string(self, position):
        return get_sexadecimal_string(position)


    def get_ra_dec_string(self, position):
        return "(%f,%+f) degrees" % (position.ra.to(u.deg).value,position.dec.to(u.deg).value)


    def __push_message_buffer(self,msg):
        self.message_buffer += msg+"\n"


    def __pop_message_buffer(self):
        msg = self.message_buffer
        self.message_buffer = ""
        return msg


    def sprint(self, msg, diagnostic_msg=None, show_caller=False, is_traceback=True, buffer=True):
        my_name   = type(self).__name__ + (f"[{sys._getframe(1).f_code.co_name}]" if show_caller else "")
        my_pid    = "" if self.pid is None else f"pid={self.pid}"
        my_filter = (lambda f: "" if f is None else f"filter='{f.name}'")(self.get_filter_setting())
        prefix = f"{my_name}({my_pid}{'' if my_pid=='' or my_filter=='' else ','}{my_filter})"
        if not (diagnostic_msg is None):
            if isinstance(diagnostic_msg,fits.header.Header):
                msg_str = msg + ("\nHEADER:\n>%s" % "\n>".join(get_header_pretty_string(diagnostic_msg).splitlines()))
            else:
                msg_str = msg + ("\nTRACEBACK:" if is_traceback else "") + "\n>%s" % "\n> ".join(diagnostic_msg.splitlines())
        elif isinstance(msg,fits.header.Header):
            msg_str = "HEADER:\n>%s" % "\n>".join(get_header_pretty_string(msg).splitlines())
        else:
            msg_str = msg
        prefixed_output = "\n".join([f"{prefix}: {s}" for s in msg_str.splitlines()])
        if buffer:
            self.__push_message_buffer(prefixed_output)
        return prefixed_output


    def set_print_to_stdout(self):
        self.print_to_stdout = True
        return self


    def unset_print_to_stdout(self):
        self.print_to_stdout = False
        return self


    def print(self, msg, diagnostic_msg=None, show_caller=False, is_traceback=True, buffer=True):
        #print(self.sprint(**{key: value for key, value in locals().items() if key not in 'self'}))
        message = self.sprint(**{key: value for key, value in locals().items() if key not in 'self'})
        if self.print_to_stdout:
            print(message)


    def pack(self, url, payload=None):
        if payload:
            data = urllib.parse.urlencode(payload).encode('utf-8')
            if self.http is None:
                request = urllib.request.Request(url, data)
            else:
                request = url+"/post?"+data.decode('utf-8')
        else:
            request = url
        return request


    # get data over http post
    def send_request(self, url, payload=None):

        #potential_retries = 5
        potential_retries = self.http_request_retries

        #TODO make sure cutoutserver didnt break with changes

        request = self.pack(url,payload)
        print("sending request for fits")
        while potential_retries > 0:
            try:
                if self.http is None:
                    #response = urllib.request.urlopen(request)
                    #webserver handles own process pool
                    response = requests.get(url, verify=False, timeout=5)
                else:
                    response = self.http.request('GET',request)

            except urllib.error.HTTPError as e:
                self.print(f"{e}",is_traceback=True)
            except ConnectionResetError as e:
                self.print(f"{e}",is_traceback=True)
            except urllib3.exceptions.MaxRetryError as e:
                self.print(f"{e}",is_traceback=True)
            except requests.exceptions.ReadTimeout as e:
                self.print(f"{e}", is_traceback=True)
            except requests.exceptions.ConnectionError as e:
                raise Exception(str(e))
            except Exception as e:
                print("OTHER exception" + str(e))
            else:
                try:
                    if self.http is None:
                        #response_data = bytearray(response.read())
                        #webserver way
                        response_data = bytearray(response.content)
                    else:
                        response_data = bytearray(response.data)
                    return response_data
                except Exception as e:
                    print(f"{e}",is_traceback=True)
            self.print(f"Taking a {self.http_wait_retry_s}s nap...")
            sleep(self.http_wait_retry_s)
            self.print("OK, lest trying fetching the cutout -- again!")
            potential_retries -= 1

        print(f"WARNING: Bailed on fetch '{url}'")
        self.processing_status = processing_status.bailed
        return None


    def standardize_fits_header_DATE_and_DATE_OBS_fields(self, date_obs_value):
        # standardize formating to 'yyyy-mm-ddTHH:MM:SS[.sss]': cf., 'https://heasarc.gsfc.nasa.gov/docs/fcg/standard_dict.html'.
        return re.sub(r"^(\d\d\d\d)(\d\d)(\d\d)$",r"\1-\2-\3T00:00:00.000",sanitize_fits_date_fields(date_obs_value))
        # TODO (Issue #6): consider replacing above with and replacing sanitize_fits_date_fields in this code with this (standardize_fits_header_DATE_and_DATE_OBS_fields)
        #return sanitize_fits_date_fields(date_obs_value)

    # tries to create a fits file from bytes.
    # on fail it returns None
    def create_fits(self, data):
        # a pretend file in memory
        fits_file = io.BytesIO(data)
        fits_file.seek(0)

        # open fits hdul
        try:
            hdul = fits.open(fits_file)
        except OSError as e:
            e_s = re.sub(r"\.$","",f"{e}")
            #self.print("Badly formatted FITS file: {0}\n\treturning None".format(str(e)), file=sys.stderr)
            self.print(f"OSError: {e_s}: Badly formatted FITS file: Cutout not found: Skipping...")
            self.processing_status = processing_status.corrupted
            return None

        # get/check header field
        header = hdul[0].header
        if ('NAXIS'  in header and header['NAXIS']  <  2) or \
           ('NAXIS1' in header and header['NAXIS1'] == 0) or \
           ('NAXIS2' in header and header['NAXIS2'] == 0):
            self.print(f"WARINING: Ill-defined 'NAXIS/i': {'NAXIS=%d => no cutout found:' % header['NAXIS'] if 'NAXIS' in header else ''} skipping...",header)
            self.processing_status = processing_status.corrupted
            return None

        # get/check data field
        data = hdul[0].data
        if data.min() == 0 and data.max() == 0:
            self.print("WARNING: Fits file contains no data: skipping...",header)
            self.processing_status = processing_status.corrupted
            return None

        # sanitize the date-obs field
        if 'DATE-OBS' in hdul[0].header:
            hdul[0].header['DATE-OBS'] = sanitize_fits_date_fields(hdul[0].header['DATE-OBS'])
        elif 'MJD-OBS' in hdul[0].header:
            hdul[0].header['DATE-OBS'] = Time(hdul[0].header['MJD-OBS'],format='mjd').isot

        # sanitize the rotation matrix
        # TODO (Issue #6): check for other antiquated stuff, i.e.,
        # http://tdc-www.harvard.edu/software/wcstools/cphead.wcs.html
        rotation_matrix_map = {
            'PC001001': 'PC1_1',
            'PC001002': 'PC1_2',
            'PC002001': 'PC2_1',
            'PC002002': 'PC2_2'
        }
        # for key in rotation_matrix_map.keys():
        #     if key in hdul[0].header:
        #         hdul[0].header.insert(key,(rotation_matrix_map[key],hdul[0].header[key]),after=True)
        #         hdul[0].header.remove(key)

        return hdul


    def get_fits(self, url, payload=None):
        self.print(f"Fetching: {url}")
        response = self.send_request(url, payload)
        return self.create_fits(response)


    def get_tiles(self, position, size):
        self.print(f"getting tile urls for {str(position)}\n\n\n\n" )
        self.request_urls_stack = self.get_tile_urls(position,size)
        if len(self.request_urls_stack) > 0:
            hdul_list = [hdul for hdul in [self.get_fits(url) for url in self.request_urls_stack] if hdul]
        else:
            self.processing_status = processing_status.none
            return None
        return hdul_list


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

        # TODO (Issue #11): So much for thread safe!
        ## *** IO_WRAPPER ***
        ## TODO (Issue #11): [1] Stream this
        ##       [2] Make a decorator
        ##       [3] Check if thread-safe
        ## setup the environment
        #old_stdout = sys.stdout
        #sys.stdout = TextIOWrapper(BytesIO(), sys.stdout.encoding)

        try:
            for i, c in enumerate(cutouts):
                with open('{directory}/{name}.fits'.format(directory=input_dir, name=i), 'wb') as tmp:
                    tmp.write(bytes(c))

            #os.listdir(td)

            # ok, let's mosaic!
            #montage.mosaic(input_dir, output_dir, bitpix=-64)
            montage.mosaic(input_dir, output_dir)

            #fits_metadata_file = f"{output_dir}/metadata.txt"
            #montage.commands.mImgtbl(input_dir,fits_metadata_file)
            #
            #header_template_file = f"{output_dir}/header.fits"
            #montage.mMakeHdr(fits_metadata_file,header_template_file)
            #with open(header_template_file,'r') as f:
            #    header_template = f.read()
            #self.print("HEADER_TEMPLATE:",header_template)

            with open('{outdir}/mosaic.fits'.format(outdir=output_dir), 'rb') as f:
                merged = f.read()
        finally:
            shutil.rmtree(output_dir)
            shutil.rmtree(input_dir)
            shutil.rmtree(td)

        ## *** IO_WRAPPER ***
        ## get output
        #sys.stdout.seek(0)      # jump to the start
        #out = sys.stdout.read() # read output

        ## *** IO_WRAPPER ***
        ## restore stdout
        #sys.stdout.close()
        #sys.stdout = old_stdout

        ## *** IO_WRAPPER ***
        ## Now we can print the output from montage
        #self.print(out)

        fits_file = io.BytesIO(merged)
        hdul = fits.open(fits_file)

        # debug
        #self.print("\n\n\n\nMOSAIKED:",hdul[0].header)

        # TODO (Issue #8): Still require investigation...
        #      but much better -- I think.
        mosiacked_header = hdul[0].header
        #wcs_header = HeaderFilter(hdul[0].header,is_add_wcs=True).get_header()
        wcs_header = WCS(hdul[0].header).to_header()
        hdul[0].header = wcs_header
        #if type(self).__name__ == 'PanSTARRS':
        #    self.print("MOSAICKED:",mosiacked_header)
        #    self.print("\n\n\n\nWCS:",wcs_header)

        # debug
        #self.print("\n\n\n\nWCS'd:",hdul[0].header)

        return hdul


    def get_image(self, hdu):
        # img_data = np.squeeze(hdu[0].data) # commenting this out doesn't change anything for Falon's code???
        img_data = hdu[0].data
        # we need to center the pixel ref's so as to reduce rotations during mosaicking
        # TODO (Issue #8):
        #      OK, so the scripette below removes the rotation,
        #         PC1_1   =     0.99982071950643 / Coordinate transformation matrix element
        #         PC1_2   =   -0.018934857951454 / Coordinate transformation matrix element
        #         PC2_1   =    0.018934857951454 / Coordinate transformation matrix element
        #         PC2_2   =     0.99982071950643 / Coordinate transformation matrix element
        #      but leaves a slight skew [1],
        #         PC1_2   = -2.4492935982947E-16 / Coordinate transformation matrix element
        #         PC2_1   =  2.4492935982947E-16 / Coordinate transformation matrix element
        #      for the PanSTARRS case (Ra=181.416667 deg, Dec=49.177778 deg, size=3').
        #      Investigate...
        #      --
        #      [1] Similar skewing was observed to be introduced for the WISE case (Ra=164.811020,
        #          dec=5.292027, size=3'), which did not have small rotation issues. This should,
        #          therefore, have very little impact on the other -- non-PanSTARRS --  mosaicking
        #          cases.

        # commenting these lines fixed Mosaick issue in webserver
        # hdu[0].header = HeaderFilter(hdu[0].header,is_add_wcs=True).update({
        #     'CRPIX1': (np.round(len(hdu[0].data[0])/2.0,1), 'Axis 1 reference pixel'),
        #     'CRPIX2': (np.round(len(hdu[0].data)/2.0,1), 'Axis 2 reference pixel')
        # }).get_header()

        # debug
        #self.print(WCS(hdu[0].header).to_header())

        img = fits.PrimaryHDU(img_data, header=hdu[0].header)

        mem_file = io.BytesIO()
        img.writeto(mem_file)
        mem_file.seek(0)
        return mem_file.getvalue()


    def paste_tiles(self, hdul_tiles, position):
        if hdul_tiles is None:
            return None

        hdu = None
        if len(hdul_tiles) > 1:
            self.print(f"Pasting {len(hdul_tiles)} at J{self.get_sexadecimal_string(position)}")
            try:
                imgs = [img for img in [self.get_image(tile) for tile in hdul_tiles]]
                # TODO (Issue #6): Need to handle multiple COADDID's...
                header_template = hdul_tiles[0][0].header
                hdu = self.mosaic(imgs)[0]
                # TODO: This requires vetting
                new_keys = set([k for k in hdu.header.keys()]+['PC1_1','PC1_2','PC2_1','PC2_2'])
                old_keys = set([k for k in header_template.keys()])
                for key in new_keys:
                    if key in old_keys:
                        header_template.remove(key)
                for key in hdu.header:
                    if key != 'COMMENT':
                        header_template[key] = hdu.header[key]
                # we need to have the header properly centered for trimming
                header_template = HeaderFilter(header_template,is_add_wcs=True).update({
                    'CRVAL1': (position.ra.to(u.deg).value, 'RA at reference pixel'),
                    'CRVAL2': (position.dec.to(u.deg).value, 'Dec at reference pixel')
                }).update({
                    'CRPIX1': (np.round(len(hdu.data[0])/2.0,1), 'Axis 1 reference pixel'),
                    'CRPIX2': (np.round(len(hdu.data)/2.0,1), 'Axis 2 reference pixel')
                }).get_header()
                hdu.header = header_template
            except montage.status.MontageError as e:
                self.print(f"Mosaicking Failed: {e}:",diagnostic_msg=traceback.format_exc(),show_caller=True)
                self.processing_status = processing_status.error
            except OSError as e:
                self.print(f"Badly formatted FITS file: {e}:",diagnostic_msg=traceback.format_exc(),show_caller=True)
                self.processing_status = processing_status.error
        elif len(hdul_tiles) == 1: # this shouldn't ever happen with webserver code
                hdu = hdul_tiles[0][0]
        self.mosaic_hdul_tiles_stack = hdul_tiles
        return hdu


    def trim_tile(self, hdu, position, size):
        if hdu is None:
            return None

        w = WCS(hdu.header)

        # trim to 2d from nd
        naxis = w.naxis
        while naxis > 2:
            w = w.dropaxis(2)
            naxis -= 1
        img_data = np.squeeze(hdu.data)

        stamp = Cutout2D(img_data, position, size, wcs=w, mode='trim', copy=True)
        hdu.header.update(stamp.wcs.to_header())
        trimmed = fits.PrimaryHDU(stamp.data, header=hdu.header)

        return trimmed


    def format_fits_hdu(self, hdu, position, size):
        if hdu is None:
            return None

        # hdu stuff...
        hdf = HeaderFilter(hdu.header,is_add_wcs=True).update({'DATE-OBS': ('na','')},is_overwrite_existing=False)
        data = hdu.data

        # standardize DATE-OBS to 'yyyy-mm-ddT00:00:00.000'
        hdf.update({
            'DATE-OBS': self.standardize_fits_header_DATE_and_DATE_OBS_fields(hdf.get_header()['DATE-OBS'])
        })

        # TODO (Issue #22): Need a cleaner way of centerng on the image, because NVSS truncates result to size of field,
        #       voiding the calculation below; however, it comes centered -- I think.
        # set (ra,dec) tile center ... rounded to 5dp -- more than good enough
        ra  = np.round(position.ra.to(u.deg).value,5)
        dec = np.round(position.dec.to(u.deg).value,5)
        hdf.update({
            'CRVAL1': (ra, 'RA at reference pixel'),
            'CRVAL2': (dec, 'Dec at reference pixel')
        })

        # set pixel reference position to center of tile
        x_pixels = len(data[0])
        y_pixels = len(data)
        hdf.update({
            'CRPIX1': (np.round(x_pixels/2.0, 1), 'Axis 1 reference pixel'),
            'CRPIX2': (np.round(y_pixels/2.0, 1), 'Axis 2 reference pixel')
        })


        # set survey name
        survey = type(self).__name__
        hdf.update({
            'SURVEY': (survey, 'Survey image obtained from')
        })

        # set default band
        hdf.update({
            'BAND': 'na',
        }, is_overwrite_existing=False)

        # TODO (Issue #6): Same thing as EQUINOX ... depricate.
        ## set epoch
        #hdf.update({
        #    'EPOCH': (2000.0, 'Julian epoch of observation')
        #}, is_overwrite_existing=False)

        # add survey dependent stuff...
        #hdf.update(self.get_fits_header_updates(hdf.get_header(),position,size))
        header_updates = self.get_fits_header_updates(hdf.get_header())
        if not header_updates is None and 'COMMENT' in header_updates:
            comment_updates = re.sub(r"\s*$"," ",header_updates['COMMENT'])
            del header_updates['COMMENT']
        else:
            comment_updates = ""
        hdf.update(header_updates)

        # set up new hdu
        ordered_keys   = hdf.saved_keys
        updated_header = hdf.get_header()
        new_hdu = fits.PrimaryHDU(data)
        for k in ordered_keys:
            new_hdu.header[k] = (updated_header[k], updated_header.comments[k])

        # add custom comments
        #new_hdu.header['COMMENT'] = ('This cutout was by the VLASS cross-ID working group within the CIRADA   project (www.cirada.ca)')
        comment_updates += 'This cutout was by the VLASS cross-ID working group within the CIRADA project (www.cirada.ca)'
        new_hdu.header['COMMENT'] = (comment_updates)

        # TODO (Issue #4): This is a kludge, for now, so as to make
        #       the claRAN machine learning code not bail.
        for field in ['LATPOLE','LONPOLE']:
            if field in new_hdu.header:
                del new_hdu.header[field]

        #self.print(f"  ===>  new_header:")
        #self.print(f"{get_header_pretty_string(new_hdu.header)}")

        return new_hdu


    def __pop_request_urls_stack(self):
        request_urls_stack = self.request_urls_stack
        self.request_urls_stack = list()
        return request_urls_stack


    def __pop_mosaic_hdul_tiles_stack(self):
        mosaic_hdul_tiles_stack = self.mosaic_hdul_tiles_stack
        self.mosaic_hdul_tiles_stack = list()
        return mosaic_hdul_tiles_stack


    def get_cutout(self, position, size):
        self.processing_status = processing_status.fetching
        try:
            tiles   = self.get_tiles(position,size)
            tile    = self.paste_tiles(tiles,position)
            trimmed = self.trim_tile(tile,position,size)
            cutout  = self.format_fits_hdu(trimmed,position,size)
            self.processing_status = processing_status.done
        except Exception as e:
            self.print(f"ERROR: {e}",diagnostic_msg=traceback.format_exc(),show_caller=True)
            self.processing_status = processing_status.error
            cutout = None

        self.print(f"J{self.get_sexadecimal_string(position)}[{self.get_ra_dec_string(position)} at {size}]: Processing Status = '{self.processing_status.name}'.")

        ## debug -- testing
        #for n,proc in enumerate(processing_status):
        #    if n == self.pid:
        #        self.print(f"DEBUGGING: setting status to {proc.name}...")
        #        self.processing_status = proc
        #        cutout = None
        return {
            'cutout':       cutout,
            'request_urls': self.__pop_request_urls_stack(),
            'raw_tiles':    self.__pop_mosaic_hdul_tiles_stack(),
            'message':      self.__pop_message_buffer(),
            'status':       self.__pop_processing_status()
        }


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


    @abstractmethod
    def get_fits_header_updates(self,header):
        pass
