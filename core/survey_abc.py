import io
import os
import sys
import traceback
import tempfile
import shutil
import datetime, timestring

# *** IO_WRAPPER ***
# TODO (Issue #11): ...
# Notes: https://stackoverflow.com/questions/1218933/can-i-redirect-the-stdout-in-python-into-some-sort-of-string-buffer/33979942#33979942
from io import TextIOWrapper, BytesIO
import io, os, shutil, tempfile, sys, base64, shutil, json
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
from .toolbox import *
from .FITS2DImageTools import *

from astropy import units as u

import montage_wrapper as montage
from astropy.nddata.utils import Cutout2D

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

# abstract class for a survey
from abc import ABC, abstractmethod
class SurveyABC(ABC):
    """
    The routine returns the following information:
               --------------------------------------
               | Dict Element | Contents            |
               |______________|_____________________|
               | cutout       | Processed Cutout    |
               | request_urls | Request URLS        |
               | raw_tiles    | Unprocessed Cutouts |
               | message      | Diagnotics Message  |
               | status       | Processing Status   |
               --------------------------------------
    """
    def __init__(self):# http_pool_manager = None, pid = None):
        ABC.__init__(self)
        self.processing_status = processing_status.idle
        self.pid = None
        self.http = None #http_pool_manager #None
        self.message_buffer = ""
        self.print_to_stdout = True
        self.request_urls_stack = list()
        self.mosaic_hdul_tiles_stack = list()
        # http request settings
        self.http_request_retries = 5
        self.http_wait_retry_s = 5
        self.tmp_dir = "/tmp"
        self.out_dir = None

    # save a list of dicts for HDU results into a folder with originals in nearby folder
    #save_orig_separately is for webserver
    @staticmethod
    def save_and_serialize(all_fits, originals_path_end="_ORIGINALS", save_orig_separately=False, save_dir="data_out/"):
        for f_dict in all_fits:
            if f_dict["out_dir"]:
                save_at = os.path.join(f_dict["out_dir"], f_dict['filename'])
            else:
                save_at = os.path.join(save_dir, f_dict['filename'])
            if f_dict['download']:
                try:
                    f_dict['download'].writeto(save_at, overwrite=True)
                except Exception as e:
                    print("EXCEPTION in TASKS", str(e))
                    if "Verif" in str(e): #skip verify warnings from fits standards
                        f_dict['download'].writeto(save_at, overwrite=True, output_verify='silentfix+warn')
                    else:
                        raise Exception(f"problem creating {f_dict['filename']} for {f_dict['survey']}")
                # add original raw tiles as extensions to mosaic
                if len(list(f_dict['originals']))>1 and save_orig_separately==False:
                    fits_img = fits.open(f_dict['filename'], mode='append')
                    for og in f_dict['originals']:
                        fits_img.append(f_dict['originals'][og]['tile'])
                        del f_dict['originals'][og]['tile']
                    fits_img.close(output_verify="silentfix")

                elif len(list(f_dict['originals']))>1 and save_orig_separately==True:
                    orig_dir = save_dir + originals_path_end
                    if os.path.exists(orig_dir):
                        shutil.rmtree(orig_dir)
                    os.makedirs(orig_dir)
                    # sort by date-obs
                    sorted_keys = sorted((list(f_dict['originals'])), key=lambda x: f_dict['originals'][x]['obs-date'])
                    for num, url in enumerate(sorted_keys, 1):
                        fname = str(num)+"-" + url.split('/')[-1]+'.fits'
                        f_dict['originals'][url]['tile'].writeto(orig_dir+'/'+fname, overwrite=True, output_verify='silentfix+warn')
                        del f_dict['originals'][url]['tile'] # remove fits image so json serializable
                else:  # still remove fits imagset_out_dire so json serializable
                    del f_dict['originals'][list(f_dict['originals'])[0]]['tile']
            else:
                print(f"Cutout at (RA, Dec) of ({f_dict['position'].ra.to(u.deg).value}, {f_dict['position'].dec.to(u.deg).value}) degrees /w size={f_dict['radius']} returned None for FITS data.",buffer=False)
            f_dict['download_path'] = save_at
            f_dict['thumbnail'] = base64.encodestring(get_thumbnail(f_dict['download'], f_dict['survey'])).decode('ascii')
            del f_dict['download'] # don't need FITZ image in memory because saved locally now

        return all_fits


    def __pop_processing_status(self):
        status = self.processing_status
        self.processing_status = processing_status.idle
        return status

    def set_tmp_dir(self, directory):
        self.tmp_dir= directory
        return self

    def set_out_dir(self, dir_path):
        self.out_dir = dir_path

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

    def __push_message_buffer(self,msg):
        self.message_buffer += msg+"\n"

    def __pop_message_buffer(self):
        msg = self.message_buffer
        self.message_buffer = ""
        return msg

    # std out printing functions
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
        message = self.sprint(**{key: value for key, value in locals().items() if key not in 'self'})
        if self.print_to_stdout:
            print(message)

    # currently only survey that packs data is NVSS
    def pack(self, url, payload=None):
        if payload:
            data = urllib.parse.urlencode(payload).encode('utf-8')
            # if self.http is None:
            #     request = urllib.request.Request(url, data)
            # else:
            request = url+"/post?"+data.decode('utf-8')
        else:
            request = url
        return request

    # get data over http post
    def send_request(self, url):
        potential_retries = self.http_request_retries
        print("sending request for fits")
        while potential_retries > 0:
            try:
                #webserver handles own process pool
                if self.http is None:
                    #response = urllib.request.urlopen(request)
                    response = requests.get(url, verify=False, timeout=40)
                else:
                    response = self.http.request('GET',url)
            except urllib.error.HTTPError as e:
                self.print(f"{e}",is_traceback=True)
            except ConnectionResetError as e:
                self.print(f"{e}",is_traceback=True)
            except urllib3.exceptions.MaxRetryError as e:
                self.print(f"{e}",is_traceback=True)
            except requests.exceptions.ReadTimeout as e:
                self.print(f"{e}", is_traceback=True)
            except requests.exceptions.ConnectionError as e:
                self.print(f"{e}", is_traceback=True)
                # raise Exception(str(e))
            except Exception as e:
                self.print("OTHER exception" + str(e))
            else:
                try:
                    #webserver way
                    if self.http is None:
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
            hdul = fits.open(fits_file, ignore_missing_end= True)
        except OSError as e:
            if "data is available" in str(data):
                raise Exception(f"{type(self).__name__}: error creating FITS "+ str(data.decode()))
            else:
                e_s = re.sub(r"\.$","",f"{e}")
                #self.print("Badly formatted FITS file: {0}\n\treturning None".format(str(e)), file=sys.stderr)
                self.print(f"OSError: {e_s}: Badly formatted FITS file: Cutout not found: Skipping...", is_traceback=True)
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
        # rotation_matrix_map = {
        #     'PC001001': 'PC1_1',
        #     'PC001002': 'PC1_2',
        #     'PC002001': 'PC2_1',
        #     'PC002002': 'PC2_2'
        # }
        # for key in rotation_matrix_map.keys():
        #     if key in hdul[0].header:
        #         hdul[0].header.insert(key,(rotation_matrix_map[key],hdul[0].header[key]),after=True)
        #         hdul[0].header.remove(key)
        return hdul

    def get_fits(self, url):
        self.print(f"Fetching: {url}")
        try:
            response = self.send_request(url)
        except Exception as e:
            print(self, "exception in surveyabc" + str(e))
            print(f"{type(self).__name__} EXCEPTION" + str(e))
            raise Exception(f"{type(self).__name__} EXCEPTION: " + str(e))
        if not response:
            print(self, "exception in surveyabc no response" + str(e))
            raise Exception(f"No FITS found at url {url} survey {type(self).__name__} !")
        hdul = self.create_fits(response)
        if not hdul:
            print(self, "NO HDUL" + str(e))
            raise Exception(f"{type(self).__name__}: error creating FITS")
        return (hdul[0], url)

    def get_tiles(self, position, size):
        self.print(f"getting tile urls for {str(position)}\n" )
        self.request_urls_stack = self.get_tile_urls(position,size)
        if not self.request_urls_stack:
            self.processing_status = processing_status.none
            raise Exception(f"no valid {type(self).__name__} urls found")
            return None
        if len(self.request_urls_stack) > 0:
            hdul_list = [hdul_tup for hdul_tup in [self.get_fits(url) for url in self.request_urls_stack] if hdul_tup[0]]
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
        tempfile.tempdir = self.tmp_dir
        td = tempfile.mkdtemp()
        input_dir = '{directory}/input'.format(directory=td)
        output_dir = '{directory}/output'.format(directory=td)
        self.__make_dir(input_dir)
        try:
            for i, c in enumerate(cutouts):
                with open('{directory}/{name}.fits'.format(directory=input_dir, name=i), 'wb') as tmp:
                    tmp.write(bytes(c))

            montage.mosaic(input_dir, output_dir, subset_fast =True)

            with open('{outdir}/mosaic.fits'.format(outdir=output_dir), 'rb') as f:
                merged = f.read()
        finally:
            shutil.rmtree(output_dir)
            shutil.rmtree(input_dir)
            shutil.rmtree(td)

        fits_file = io.BytesIO(merged)
        hdul = fits.open(fits_file)

        # TODO (Issue #8): Still require investigation...
        #      but much better -- I think.
        #wcs_header = HeaderFilter(hdul[0].header,is_add_wcs=True).get_header()
        wcs_header = WCS(hdul[0].header).to_header()
        hdul[0].header = wcs_header
        return hdul

    def get_image(self, hdu):
        # img_data = np.squeeze(hdu[0].data) # commenting this out doesn't change anything for Falon's code???
        img_data = hdu.data
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
        img = fits.PrimaryHDU(img_data, header=hdu.header)
        mem_file = io.BytesIO()
        img.writeto(mem_file)
        mem_file.seek(0)
        return mem_file.getvalue()

    def paste_tiles(self, hdul_tiles, position):
        if hdul_tiles is None:
            return None

        hdu = None
        if len(hdul_tiles) > 1:
            try:
                imgs = [img for img in [self.get_image(tile) for (tile,url) in hdul_tiles]]
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

    def format_fits_hdu(self, hdu, position, all_headers):
        if hdu is None:
            return None
        # hdu stuff...
        hdf = HeaderFilter(hdu.header,is_add_wcs=True).update({'DATE-OBS': ('na','')},is_overwrite_existing=False)
        data = hdu.data

        # standardize DATE-OBS to 'yyyy-mm-ddT00:00:00.000'
        date_obs = re.sub(r"^(\d\d\d\d)(\d\d)(\d\d)$",r"\1-\2-\3T00:00:00.000", sanitize_fits_date_fields(hdf.header['DATE-OBS']))
        hdf.update({
            'DATE-OBS': date_obs
        })
        # set (ra,dec) tile center ... rounded to 5dp -- more than good enough
        ra  = np.round(position.ra.to(u.deg).value,7)
        dec = np.round(position.dec.to(u.deg).value,7)
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

        # add survey dependent stuff...
        hdf.update(self.get_fits_header_updates(hdf.header, all_headers))

        # set up new hdu
        ordered_keys   = hdf.saved_keys
        updated_header = hdf.header
        new_hdu = fits.PrimaryHDU(data)
        for k in ordered_keys:
            if k.upper()=="COMMENT" or k.upper()=="HISTORY":
                new = str(updated_header[k]).replace("\n", " ")
                new_hdu.header[k] = (new, updated_header.comments[k])
            else:
                new_hdu.header[k] = (updated_header[k], updated_header.comments[k])
        # add custom comments
        # add more custom comments
        new_hdu.header['COMMENT'] = ('This cutout was provided by the CIRADA project (www.cirada.ca)')
        new_hdu.header['COMMENT'] = "The Canadian Initiative for Radio Astronomy Data Analysis (CIRADA) is funded " \
                        "by a grant from the Canada Foundation for Innovation 2017 Innovation Fund (Project 35999) " \
                        "and by the Provinces of Ontario, British Columbia, Alberta, Manitoba and Quebec, in " \
                        "collaboration with the National Research Council of Canada, the US National Radio Astronomy " \
                        "Observatory and Australia's Commonwealth Scientific and Industrial Research Organisation."
        new_hdu.header['COMMENT'] = "BMAJ, BMIN, BPA, MJD-OBS, and DATE-OBS only currently represent the values" \
                                    " from one of the input files."

        # TODO (Issue #4): Is this still needed?? This is a kludge, for now, so as to make
        #       the claRAN machine learning code not bail.
        # for field in ['LATPOLE','LONPOLE']:
        #     if field in new_hdu.header:
        #         del new_hdu.header[field]
        return new_hdu

    def __pop_request_urls_stack(self):
        request_urls_stack = self.request_urls_stack
        self.request_urls_stack = list()
        return request_urls_stack

    def __pop_mosaic_hdul_tiles_stack(self):
        mosaic_hdul_tiles_stack = self.mosaic_hdul_tiles_stack
        self.mosaic_hdul_tiles_stack = list()
        return mosaic_hdul_tiles_stack

    def group_tiles(self, tile_tups, rule):
        # rule has to be a valid fits HEADER value (such as DATE-OBS) OR of our defined groups "Mosaic" or "None"
        # tile_tups is a list of tuples with (hdul,url) for each
        groups = {}
        if rule=="MOSAIC" or rule=="None": # None is default
            return {rule:tile_tups}
        for (tile,tile_url) in tile_tups:
            if rule in tile.header:
                this_group = tile.header[rule]
                if rule=='DATE-OBS':
                    this_group = timestring.Date(tile.header[rule]).date.strftime("%Y-%m-%d")
                if this_group in list(groups):
                    groups[this_group].append((tile,tile_url))
                else:
                    groups[this_group]=[ (tile,tile_url) ]
            else: # None is default, don't group if category not found in header
                if "None" in list(groups):
                    groups["None"].append((tile,tile_url))
                else:
                    groups["None"]=[ ((tile,tile_url)) ]
        return groups

    def process_tile_group(self, tiles, position, size, group, index):
        fits_data = {}
        fits_data['epoch'] = None
        survey_name = type(self).__name__
        radius = size/2
        try:
            filter = self.get_filter_setting().name.lower()
        except Exception as e:
            print("no filter name" + str(e))
            filter = ""
        # try:
        if len(tiles)>1:
            all_headers = [t.header for (t, tile_url) in tiles]
            tile    = self.paste_tiles(tiles, position)
            if self.needs_trimming:
                tile = self.trim_tile(tile,position,size)
            cutout  = self.format_fits_hdu(tile,position,all_headers)
            fits_data['filename'] = get_mosaic_filename(position,radius,survey_name, filter=filter, group_title=group, extension=None)
        elif len(tiles)==0:
            print(f"no {survey_name} tiles found for {position}! ")
            return None
        else:
            if group == "MOSAIC":
                group="None"
            fits_data['filename'] = get_non_mosaic_filename(position, radius, survey_name, baseurl=tiles[0][1], index=index)
            cutout = tiles[0][0]
            if self.needs_trimming:
                cutout = self.trim_tile(cutout,position,size)
            if survey_name=="VLASS": # only label if not mosaicked for now in case multiple epochs
                fits_data['epoch'] = self.get_epoch(fits_data['filename'])

        # store original tiles
        # THIS SPECIFIC TO VLASS ONLY stokes being FILNAME09
        if survey_name == "VLASS":
            fits_data['originals'] = {tile_url:{'obs-date': get_header_value(tile, 'DATE-OBS'),
                                                'stokes': get_header_value(tile, 'FILNAM09'),
                                                'epoch': self.get_epoch(tile_url),
                                                'tile': tile} for (tile,tile_url) in tiles}
        else:
            fits_data['originals'] = {tile_url:{'obs-date': get_header_value(tile, 'DATE-OBS'),
                                                'tile': tile} for (tile,tile_url) in tiles}
        # except Exception as e:
        #     self.print(f"ERROR: {e}",diagnostic_msg=traceback.format_exc(),show_caller=True)
        #     self.processing_status = processing_status.error
            #return None

        self.processing_status = processing_status.done
        fits_data['download'] = cutout
        fits_data['out_dir'] = self.out_dir
        fits_data['group'] = group
        fits_data['survey'] = survey_name
        fits_data['filter'] = filter
        fits_data['position'] = f"{position.ra.degree}, {position.dec.degree}"
        fits_data['radius'] = radius.value
        # fits_data['status'] = self.__pop_processing_status()
        self.print(f"[Position: {position.ra.degree}, {position.dec.degree} at radius {size/2}]: Processing Status = '{self.processing_status.name}'.")
        return fits_data

        # return {
        #     'cutout':       cutout,
        #     'request_urls': self.__pop_request_urls_stack(),
        #     'raw_tiles':    self.__pop_mosaic_hdul_tiles_stack(),
        #     'message':      self.__pop_message_buffer(),
        #     'status':       self.__pop_processing_status()
        # }


    # main routine for CLI cutout processing
    def get_cutout(self, position, size, group_by="None"):
        if not group_by:
            group_by="None"
        self.processing_status = processing_status.fetching
        tiles   = self.get_tiles(position,size)
        if not tiles:
            print("NO {type(self).__name__} TILES FOUND for {position}")
        groups_dict = self.group_tiles(tiles, group_by) # to read header and separate tiles
        all_fits = []
        for group in list(groups_dict):
            if group=="None": # handle each individually if no grouping
                for single in groups_dict[group]:
                    # THIS COULD BE DANGEROUS? WILL THERE BE DOUBLES? THEY BE OVERWRITTEN WITH UPDATE
                    all_fits.append(self.process_tile_group([single], position, size, "None", groups_dict[group].index(single)))
            else:
                all_fits = all_fits+[self.process_tile_group(groups_dict[group], position, size, group, 0)]
        return all_fits

    # abstract base class functions required by survey/child classes
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
    def get_fits_header_updates(self,header, all_headers=None):
        pass
