import io
import os
import sys
import traceback
import tempfile
import shutil

# *** IO_WRAPPER ***
# Notes: https://stackoverflow.com/questions/1218933/can-i-redirect-the-stdout-in-python-into-some-sort-of-string-buffer/33979942#33979942
from io import TextIOWrapper, BytesIO

import urllib.request
import urllib.parse
import urllib.error
import urllib3
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

import montage_wrapper
from astropy.nddata.utils import Cutout2D

from enum import Enum
class processing_status(Enum):
    idle      = 0
    none      = 1
    fetching  = 2
    corrupted = 3
    error     = 4
    bailed    = 5
    done      = 6

# abstract class for a survey
from abc import ABC, abstractmethod
class SurveyABC(ABC):
    def __init__(self,
        http_pool_manager = None, # thread safe pool manager (i.e., urllib3)
        pid = None # 4 threads
    ):
        super().__init__()

        self.processing_status = processing_status.idle

        self.pid = None
        self.http = None

        self.message_buffer = ""


    def __pop_processing_status(self):
        status = self.processing_status
        self.processing_status = processing_status.idle
        return status


    def set_pid(self, pid):
        self.pid = pid
        return self


    def attach_http_pool_manager(self,http_pool_manager):
        self.http = http_pool_manager


    def get_sexadecimal_string(self, position):
        sexadecimal = "%02d%02d%02.0f" % position.ra.hms+re.sub(r"([+-])\d",r"\1","%+d%02d%02d%02.0f" % position.dec.signed_dms)
        return sexadecimal 


    def get_ra_dec_string(self, position):
        return "(%f,%+f) degrees" % (position.ra.to(u.deg).value,position.dec.to(u.deg).value)


    def __push_message_buffer(self,msg):
        self.message_buffer += msg+"\n"


    def __pop_message_buffer(self):
        msg = self.message_buffer
        self.message_buffer = ""
        return msg


    def sprint(self, msg, diagnostic_msg=None, show_caller=False, is_traceback=True):
        my_name    = type(self).__name__ + (f"[{sys._getframe(1).f_code.co_name}]" if show_caller else "")
        my_pid     = "" if self.pid is None else f"pid={self.pid}"
        my_filter  = (lambda f: "" if f is None else f"filter='{f.name}'")(self.get_filter_setting())
        prefix = f"{my_name}({my_pid}{'' if my_pid=='' or my_filter=='' else ','}{my_filter})"
        if not (diagnostic_msg is None):
            if isinstance(diagnostic_msg,fits.header.Header):
                msg_str = msg + ("\nHEADER:\n>%s" % "\n>".join(get_header_pretty_string(diagnostic_msg).splitlines()))
            else:
                msg_str = msg + ("\nTRACEBACK:" if is_traceback else "") + "\n>%s" % "\n> ".join(diagnostic_msg.splitlines())
        else:
            msg_str = msg
        prefixed_output = "\n".join([f"{prefix}: {s}" for s in msg_str.splitlines()])
        self.__push_message_buffer(prefixed_output)
        return prefixed_output


    def print(self, msg, diagnostic_msg=None, show_caller=False, is_traceback=True):
        print(self.sprint(msg,diagnostic_msg,show_caller,is_traceback))


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
    def __send_request(self, url, payload=None):

        potential_retries = 5

        request = self.pack(url,payload)

        while potential_retries > 0:

            try:
                if self.http is None:
                    response = urllib.request.urlopen(request)
                else:
                    response = self.http.request('GET',request)
            except urllib.error.HTTPError as e:
                pass
            except ConnectionResetError as e:
                pass
            except urllib3.exceptions.MaxRetryError as e:
                pass
            else:
                try:
                    if self.http is None:
                        response_data = bytearray(response.read())
                    else:
                        response_data = bytearray(response.data)
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

            # TODO: clean this up, as well as, retries, for configuration...
            duration_s = 60 if self.http else 25
            self.print(f"Taking a {duration_s}s nap...")
            sleep(duration_s)
            self.print("OK, lest trying fetching the cutout -- again!")

            potential_retries -= 1

        self.print(f"WARNING: Bailed on fetch '{url}'")
        self.processing_status = processing_status.bailed
        return None


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

        # TODO: So much for thread safe!
        ## *** IO_WRAPPER ***
        ## TODO: [1] Stream this
        ##       [2] Make a decorator
        ##       [3] Check if thread-safe
        ## setup the environment
        #old_stdout = sys.stdout
        #sys.stdout = TextIOWrapper(BytesIO(), sys.stdout.encoding)

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

        return merged


    def standardize_fits_header_DATE_and_DATE_OBS_fields(self, date_obs_value):
        # standardize formating to 'yyyy-mm-ddTHH:MM:SS[.sss]': cf., 'https://heasarc.gsfc.nasa.gov/docs/fcg/standard_dict.html'.
        return re.sub(r"^(\d\d\d\d)(\d\d)(\d\d)$",r"\1-\2-\3T00:00:00.000",sanitize_fits_date_fields(date_obs_value)) 

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
            self.print(f"WARINING: Ill-defined 'NAXIS/i': {'NAXIS=%d => no cutout found:' % header['NAXIS'] if 'NAXIS' in header else ''} skipping...")
            self.processing_status = processing_status.corrupted
            return None

        # get/check data field
        data = hdul[0].data
        if data.min() == 0 and data.max() == 0:
            self.print("WARNING: Fits file contains no data: skipping...")
            self.processing_status = processing_status.corrupted
            return None

        # sanitize the date-obs field
        if 'DATE-OBS' in hdul[0].header:
            hdul[0].header['DATE-OBS'] = sanitize_fits_date_fields(hdul[0].header['DATE-OBS'])

        # sanitize the rotation matrix
        # TODO: check for other antiquated stuff, i.e.,
        # http://tdc-www.harvard.edu/software/wcstools/cphead.wcs.html
        rotation_matrix_map = {
            'PC001001': 'PC1_1',
            'PC001002': 'PC1_2',
            'PC002001': 'PC2_1',
            'PC002002': 'PC2_2'
        }
        for key in rotation_matrix_map.keys():
            if key in hdul[0].header:
                hdul[0].header.insert(key,(rotation_matrix_map[key],hdul[0].header[key]),after=True)
                hdul[0].header.remove(key)
    
        return hdul


    def get_fits(self, url, position, payload=None):
        self.print(f"Fetching: {url}")
        response = self.__send_request(url, payload)
        # note that it returns None if the response isn't a valid fits
        return self.create_fits(response, position)


    def get_tiles(self, position, size):
        urls = self.get_tile_urls(position,size)
        if len(urls) > 0:
            hdul_list = [hdul for hdul in [self.get_fits(url,position) for url in urls] if hdul]
        else:
            self.processing_status = processing_status.none
            return None
        return hdul_list


    def get_image(self, hdu):
        img_data = np.squeeze(hdu[0].data)
        img = fits.PrimaryHDU(img_data, header=hdu[0].header)

        # writing to a pretend file in memory
        mem_file = io.BytesIO()
        img.writeto(mem_file)
        return mem_file.getvalue()


    def paste_tiles(self, hdul_tiles, position, size):
        if hdul_tiles is None:
            return None

        hdu = None
        if len(hdul_tiles) > 1:
            self.print(f"Pasting {len(hdul_tiles)} at J{self.get_sexadecimal_string(position)}")
            try:
                imgs = [img for img in [self.get_image(tile) for tile in hdul_tiles]]
                # TODO: Need to handle multiple COADDID's...
                header_template = hdul_tiles[0][0].header
                tiled_fits_file = io.BytesIO(self.mosaic(imgs))
                tiled_fits_file.seek(0)
                hdu = fits.open(tiled_fits_file)[0]
                for key in hdu.header:
                    if key != 'COMMENT':
                        header_template[key] = hdu.header[key]
                hdu.header = header_template
            except montage_wrapper.status.MontageError as e:
                self.print(f"Mosaicking Failed: {e}: file={sys.stderr}",show_caller=True)
                self.processing_status = processing_status.error
            except OSError as e:
                self.print("Badly formatted FITS file: {0}\n\treturning None".format(str(e)), file=sys.stderr)
                self.processing_status = processing_status.error
        elif len(hdul_tiles) == 1:
                hdu = hdul_tiles[0][0]
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
        img = fits.PrimaryHDU(stamp.data, header=hdu.header)
    
        # writing to a pretend file in memory
        #mem_file = io.BytesIO()
        #img.writeto(mem_file)
        return img


    def header_write(self, hdu, position):
        def yrmon_to_isot(yrmon):
            ###converts obs date of form 'yyyymm' to 'yyyy-mm-ddT00:00:00.000'
            ###allows atropy time to be used
            yr = yrmon[:4]
            mon = yrmon[4:6] ###set to assume str len == 6 accounts for some FIRST FITS adding dd at end
            ##can make more complex later to use dd info, at moment good enough as yyyymm will assume dd==15 T==00:00:00 for MJD
            isottime = yr+'-'+mon+'-15T00:00:00.000'
            return(isottime)

        if hdu is None:
            return None

        head = hdu.header
        data = hdu.data
        survey = type(self).__name__
        a = position.ra.to(u.deg).value
        d = position.dec.to(u.deg).value
        #head['DATE-OBS'] = self.standardize_fits_header_DATE_OBS_field(head['DATE-OBS']) if 'DATE-OBS' in head else 'na'

        ###data = cutout image data array
        ###head = original header from survey image - i.e. preserved old header
        ###survey =  survey image taken from
        ###a = target RA - not the image RA from the survey main image, this can be different
        ###d = target Dec - as RA, target RA/Dec should be obtained from the target input file
        ##round a, d to 5dp (more than good enough)
        a, d = np.round(a, 5), np.round(d, 5)
        
        ####remove superfluous nesting
        if data.ndim > 2:
            data = data[0]
        
        #set up pole_longitude
        if a>90 and a<=270:
            p_lon = 180.0
        else:
            p_lon = 0.0

        ##set up survey specific info
        if survey == 'PanSTARRS':
            sdict = {'BAND': ('i-band', 'Filter used in observation'),
                     'pos_units': 'deg',
                     'RADESYS': (head['RADESYS'], 'Coordinate system used'),
                     'DATE-OBS': (Time(head['MJD-OBS'], format='mjd').isot, 'Obs. date')}
        elif survey == 'WISE':
                sdict = {'BAND': ('W1', 'Filter used in observation'),
                         'pos_units': 'deg',
                         'RADESYS': ('FK5', 'Coordinate system used'),
                         'DATE-OBS': (head['MIDOBS'], 'Median observation date of stack')}
        elif survey == 'FIRST':
            sdict = {'BAND': ('1.4 GHz', 'Frequency of observation'),
                     'pos_units': 'deg',
                     'RADESYS': ('FK5', 'Coordinate system used'),
                     'DATE-OBS': (head['DATE-OBS'], 'Obs. date (yearmonth)'),
                     'FNAME': (head['FIELDNAM'], 'FIRST coadded image')}
        elif survey == 'VLASS': # michelle was here!
            ###complex file name - extract from header info
            fpartkeys = ['FILNAM01', 'FILNAM02', 'FILNAM03', 'FILNAM04', 'FILNAM05', 'FILNAM06',
                         'FILNAM07', 'FILNAM08', 'FILNAM09', 'FILNAM10', 'FILNAM11', 'FILNAM12']
            nameparts = [head[key] for key in fpartkeys]
            ###create single string - FILNAM12 goes after a constant
            vfile = nameparts[0]
            for i in range(len(nameparts)-2):
                vfile = vfile + '.' + nameparts[i+1]
            vfile = vfile + '.pbcor.' + nameparts[len(nameparts)-1] + '.subim.fits'
                             
            sdict = {'BAND': ('2-4 GHz', 'Frequency coverage of observation'),
                     'pos_units': head['CUNIT1'],
                     'RADESYS': (head['RADESYS'], 'Coordinate system used'),
                     'DATE-OBS': (head['DATE-OBS'], 'Obs. date'),
                     'STOKES': (head['BTYPE'], 'Stokes polarisation'),
                     'FNAME': (vfile, 'VLASS image file')}
        else:
            sdict = {'BAND': 'na', 'pos_units': 'deg', 'RADESYS': ('FK5', 'assumed'), 'DATE-OBS': 'na'}
        
        ###set up naxis len - vlass QL image headers bugger up just extracting this from header
        #####take from data instead - fixed by YG (12 Jun 2019 - 11:01)
        xlen = len(data[0])
        ylen = len(data)

        ### set up main header keys for all surveys
        ####add in FUNIT (z-axis unit, flux)?
        hkeys = {'WCSAXES': (2, 'Number of WCS axes'),
                 'CRVAL1': (a, 'RA at reference pixel'),
                 'CRVAL2': (d, 'Dec at reference pixel'),
                 'CRPIX1': (np.round(xlen/2, 1), 'Axis 1 reference pixel'),
                 'CRPIX2': (np.round(ylen/2, 1), 'Axis 2 reference pixel'),
                 'CDELT1': (-abs(head['CDELT1']), 'Axis 1 step size per pixel'),
                 'CDELT2': (head['CDELT2'], 'Axis 2 step size per pixel'),
                 'CUNIT1': (sdict['pos_units'], 'Unit for CDELT1'),
                 'CUNIT2': (sdict['pos_units'], 'Unit for CDELT2'),
                 'CTYPE1': (head['CTYPE1'], 'RA projection'),
                 'CTYPE2': (head['CTYPE2'], 'Dec projection'),
                 'LATPOLE': (d, 'Native latitude of celestial pole'),
                 'LONPOLE': (p_lon, 'Native longitude of celestial pole'),
                 'RADESYS': sdict['RADESYS'],
                 'SURVEY': (survey, 'Survey image obtained from'),
                 'BAND': sdict['BAND'],
                 'EPOCH': (2000.0, 'Julian epoch of observation'),
                 'DATE-OBS': sdict['DATE-OBS']}

        keylist = list(hkeys.keys())
        
        
        ##set up new header
        newhead = fits.PrimaryHDU(data).header
        
        ##append main keys+info to header
        for key in hkeys:
            newhead[key] = hkeys[key]

        ###add in MJD and radio specific header info
        if survey == 'FIRST':
            #newhead['MJD'] = (Time(yrmon_to_isot(newhead['DATE-OBS'])).mjd,
            #                  'Median MJD of obs month (00:00:00 on 15th)')
            newhead['MJD'] = (Time(newhead['DATE-OBS']).mjd,
                              'Median MJD of obs month (00:00:00 on 15th)')
            ###set up radio specific keys (beamsize etc)
            radkeys = {'BUNIT': ('Jy/beam', 'Pixel flux unit'),
                       'BMAJ': (head['BMAJ'], 'Beam major axis [deg]'),
                       'BMIN': (head['BMIN'], 'Beam minor axis [deg]'),
                       'BPA': (head['BPA'], 'Beam position angle')}
            rklist = list(radkeys.keys())
            ##append to header
            for key in rklist:
                newhead[key] = radkeys[key]
            
            ##add image filename (original image from survey not cutout)
            newhead['IMFILE'] = sdict['FNAME']

        elif survey == 'VLASS':
            newhead['MJD'] = (Time(newhead['DATE-OBS']).mjd, 'MJD of the observation date')
            
            ###set up radio specific keys (beamsize etc)
            radkeys = {'BUNIT': ('Jy/beam', 'Pixel flux unit'),
                       'BMAJ': (head['BMAJ'], 'Beam major axis [deg]'),
                       'BMIN': (head['BMIN'], 'Beam minor axis [deg]'),
                       'BPA': (head['BPA'], 'Beam position angle')}
            rklist = list(radkeys.keys())
            ##append to header
            for key in rklist:
                newhead[key] = radkeys[key]

            ##include polarisation info
            newhead['STOKES'] = sdict['STOKES']

            ### add original image filename
            newhead['IMFILE'] = sdict['FNAME']

        ### PanSTARRS/WISE specifics
        elif survey == 'PanSTARRS':
            newhead['STK_TYPE'] = (head['STK_TYPE'], 'PanSTARRS image stack type')
            newhead['STK_ID'] = (head['STK_ID'], 'PanSTARRS image stack ID')
            newhead['SKYCELL'] = (head['SKYCELL'], 'PanSTARRS image sky cell')
            newhead['TESS_ID'] = (head['TESS_ID'], 'PanSTARRS tesselation')
        
        ##wise specifics
        elif survey == 'WISE':
            newhead['IMFILE'] = (head['COADDID'], 'ATLAS image identifier')

        ## add comments to header
        newhead['Comment'] = ('This cutout was by the VLASS cross-ID ' +
                              'working group within the CIRADA project (www.cirada.ca)')

        ##DONE - return new header for use in output cutout fits file

        img = fits.PrimaryHDU(data,header=newhead)
        mem_file = io.BytesIO()
        img.writeto(mem_file)
        return img


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

        # TODO: Same thing as EQUINOX ... depricate.
        ## set epoch
        #hdf.update({
        #    'EPOCH': (2000.0, 'Julian epoch of observation')
        #}, is_overwrite_existing=False)

        #hdf.update({
        #    'COMMENT': ('This cutout was by the VLASS cross-ID working group within the CIRADA   project (www.cirada.ca)')
        #})

        # add survey dependent stuff...
        hdf.update(self.get_fits_header_updates(hdf.get_header(),position,size))

        # set up new hdu
        ordered_keys   = hdf.get_saved_keys()
        updated_header = hdf.get_header()
        new_hdu = fits.PrimaryHDU(data)
        for k in ordered_keys:
            new_hdu.header[k] = (updated_header[k], updated_header.comments[k])

        # add custom comments
        new_hdu.header['COMMENT'] = ('This cutout was by the VLASS cross-ID working group within the CIRADA   project (www.cirada.ca)')

        #self.print(f"  ===>  new_header:")
        #self.print(f"{get_header_pretty_string(new_hdu.header)}")

        return new_hdu


    def get_cutout(self, position, size):
        self.processing_status = processing_status.fetching
        try:
            tiles   = self.get_tiles(position,size)
            tile    = self.paste_tiles(tiles,position,size)
            trimmed = self.trim_tile(tile,position,size)
            # TODO: Remove this code, once vetted...
            # Yjan's code
            #cutout = self.header_write(trimmed,position)
            # Integrated code
            cutout = self.format_fits_hdu(trimmed,position,size)
            self.processing_status = processing_status.done
        except Exception as e:
            self.print(f"ERROR: {e}",diagnostic_msg=traceback.format_exc(),show_caller=True)
            self.processing_status = processing_status.error
            cutout = None

        #self.print(f"Finished processing J{self.get_sexadecimal_string(position)} (i.e., {self.get_ra_dec_string(position)}) cutout of size {size}.")
        self.print(f"J{self.get_sexadecimal_string(position)}[{self.get_ra_dec_string(position)} at {size}]: Procssing Status = '{self.processing_status.name}'.")

        ## debug -- testing
        #for proc in processing_status:
        #    if proc.value == self.pid:
        #        self.print(f"DEBUGGING: setting status to {proc.name}...")
        #        self.processing_status = proc
        #        cutout = None

        return {
            'cutout':  cutout,
            'message': self.__pop_message_buffer(),
            'status':  self.__pop_processing_status()
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
    def get_fits_header_updates(self,hdu,position,size):
        pass

