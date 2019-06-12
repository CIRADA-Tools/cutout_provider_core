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
        
        ## define the NAXISi=1,2 unit fields
        #if ('CUNIT1' in header) and ('CUNIT2' in header) and (header['CUNIT1'] != header['CUNIT2']):
        #    self.print(f"WARNING: ('CUNIT1'={header['CUNIT1']}) != ('CUNIT2'={header['CUNIT2']}) !")
        #elif not (('CUNIT1' in header) or ('CUNIT2' in header)):
        #     pos_units = 'deg'
        #     if 'CTYPE1' in header:
        #         header.insert('CTYPE1', ('CUNIT1', pos_units), after=True)
        #     else:
        #         self.print(f"WARNING: 'CTYPE1' undefined!")
        #         header['CUNIT1'] = pos_units
        #     if 'CTYPE2' in header:
        #         header.insert('CTYPE2', ('CUNIT2', pos_units), after=True)
        #     else:
        #         self.print(f"WARNING: 'CTYPE2' undefined!")
        #         header['CUNIT2'] = pos_units
        #
        ### set LATPOLE
        ##if not ('LATPOLE' in header):
        ##    header['LATPOLE'] = (position.dec.to(u.deg).value, 'Native latitude of celestial pole')
        #
        ### set LONPOLE
        ##if not ('LONPOLE' in header):
        ##    ra = position.ra.to(u.deg).value
        ##    pole_longitude = 180.0 if 90.0 < ra and ra <= 270.0 else 0.0
        ##    header['LONPOLE'] = (pole_longitude, 'Native longitude of celestial pole')
        #
        ### define RADESYS field
        ##if not ('RADESYS' in header):
        ##    header['RADESYS'] = ('FK5', 'assumed')
        header.update(WCS(header).to_header())

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


    def get_image(self,hdu):
        img_data = np.squeeze(hdu[0].data)
        img = fits.PrimaryHDU(img_data, header=hdu[0].header)

        # writing to a pretend file in memory
        mem_file = io.BytesIO()
        img.writeto(mem_file)
        return mem_file.getvalue()


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
        #img = fits.PrimaryHDU(stamp.data, header=stamp.wcs.to_header())
        # notes: Above gives...
        #
        #   In [157]: w.to_header()
        #   Out[157]: 
        #   WCSAXES =                    4 / Number of coordinate axes                      
        #   CRPIX1  =                137.0 / Pixel coordinate of reference point            
        #   CRPIX2  =                 34.0 / Pixel coordinate of reference point            
        #   CRPIX3  =                  1.0 / Pixel coordinate of reference point            
        #   CRPIX4  =                  1.0 / Pixel coordinate of reference point            
        #   CDELT1  =              -0.0005 / [deg] Coordinate increment at reference point  
        #   CDELT2  =               0.0005 / [deg] Coordinate increment at reference point  
        #   CDELT3  =           21875000.0 / [Hz] Coordinate increment at reference point   
        #   CDELT4  =                  1.0 / Coordinate increment at reference point        
        #   CUNIT1  = 'deg'                / Units of coordinate increment and value        
        #   CUNIT2  = 'deg'                / Units of coordinate increment and value        
        #   CUNIT3  = 'Hz'                 / Units of coordinate increment and value        
        #   CTYPE1  = 'RA---SIN'           / Right ascension, orthographic/synthesis project
        #   CTYPE2  = 'DEC--SIN'           / Declination, orthographic/synthesis projection 
        #   CTYPE3  = 'FREQ'               / Frequency (linear)                             
        #   CTYPE4  = 'STOKES'             / Coordinate type code                           
        #   CRVAL1  =          0.374722289 / [deg] Coordinate value at reference point      
        #   CRVAL2  =        1.08333333333 / [deg] Coordinate value at reference point      
        #   CRVAL3  =         1435100000.0 / [Hz] Coordinate value at reference point       
        #   CRVAL4  =                  1.0 / Coordinate value at reference point            
        #   LONPOLE =                  0.0 / [deg] Native longitude of celestial pole       
        #   LATPOLE =        1.08333333333 / [deg] Native latitude of celestial pole        
        #   RADESYS = 'FK5'                / Equatorial coordinate system                   
        #   EQUINOX =               2000.0 / [yr] Equinox of equatorial coordinates         
        #   DATE-OBS= '19951016'           / ISO-8601 observation date                      
        #   
        #   In [158]: header
        #   Out[158]: 
        #   SIMPLE  =                    T / conforms to FITS standard                      
        #   BITPIX  =                  -32 / array data type                                
        #   NAXIS   =                    2 / number of array dimensions                     
        #   NAXIS1  =                   99                                                  
        #   NAXIS2  =                   99                                                  
        #   WCSAXES =                    2 / Number of coordinate axes                      
        #   CRPIX1  =                136.0 / Pixel coordinate of reference point            
        #   CRPIX2  =                 33.0 / Pixel coordinate of reference point            
        #   CDELT1  =              -0.0005 / [deg] Coordinate increment at reference point  
        #   CDELT2  =               0.0005 / [deg] Coordinate increment at reference point  
        #   CUNIT1  = 'deg'                / Units of coordinate increment and value        
        #   CUNIT2  = 'deg'                / Units of coordinate increment and value        
        #   CTYPE1  = 'RA---SIN'           / Right ascension, orthographic/synthesis project
        #   CTYPE2  = 'DEC--SIN'           / Declination, orthographic/synthesis projection 
        #   CRVAL1  =          0.374722289 / [deg] Coordinate value at reference point      
        #   CRVAL2  =        1.08333333333 / [deg] Coordinate value at reference point      
        #   LONPOLE =                180.0 / [deg] Native longitude of celestial pole       
        #   LATPOLE =        1.08333333333 / [deg] Native latitude of celestial pole        
        #   RADESYS = 'FK5'                / Equatorial coordinate system                   
        #   EQUINOX =               2000.0 / [yr] Equinox of equatorial coordinates         
        #   DATE-OBS= '19951016'           / ISO-8601 observation date                      
        #   
        #   In [159]:
        #
        # Below appears preserve the header: i.e., cf., above with...
        #
        #   In [165]: w.to_header()
        #   Out[165]: 
        #   WCSAXES =                    4 / Number of coordinate axes                      
        #   CRPIX1  =                137.0 / Pixel coordinate of reference point            
        #   CRPIX2  =                 34.0 / Pixel coordinate of reference point            
        #   CRPIX3  =                  1.0 / Pixel coordinate of reference point            
        #   CRPIX4  =                  1.0 / Pixel coordinate of reference point            
        #   CDELT1  =              -0.0005 / [deg] Coordinate increment at reference point  
        #   CDELT2  =               0.0005 / [deg] Coordinate increment at reference point  
        #   CDELT3  =           21875000.0 / [Hz] Coordinate increment at reference point   
        #   CDELT4  =                  1.0 / Coordinate increment at reference point        
        #   CUNIT1  = 'deg'                / Units of coordinate increment and value        
        #   CUNIT2  = 'deg'                / Units of coordinate increment and value        
        #   CUNIT3  = 'Hz'                 / Units of coordinate increment and value        
        #   CTYPE1  = 'RA---SIN'           / Right ascension, orthographic/synthesis project
        #   CTYPE2  = 'DEC--SIN'           / Declination, orthographic/synthesis projection 
        #   CTYPE3  = 'FREQ'               / Frequency (linear)                             
        #   CTYPE4  = 'STOKES'             / Coordinate type code                           
        #   CRVAL1  =          0.374722289 / [deg] Coordinate value at reference point      
        #   CRVAL2  =        1.08333333333 / [deg] Coordinate value at reference point      
        #   CRVAL3  =         1435100000.0 / [Hz] Coordinate value at reference point       
        #   CRVAL4  =                  1.0 / Coordinate value at reference point            
        #   LONPOLE =                180.0 / [deg] Native longitude of celestial pole       
        #   LATPOLE =        1.08333333333 / [deg] Native latitude of celestial pole        
        #   RADESYS = 'FK5'                / Equatorial coordinate system                   
        #   EQUINOX =               2000.0 / [yr] Equinox of equatorial coordinates         
        #   DATE-OBS= '19951016'           / ISO-8601 observation date                      
        #   
        #   In [166]: header
        #   Out[166]: 
        #   SIMPLE  =                    T / conforms to FITS standard                      
        #   BITPIX  =                  -32 / array data type                                
        #   NAXIS   =                    2 / number of array dimensions                     
        #   NAXIS1  =                   99                                                  
        #   NAXIS2  =                   99                                                  
        #   DATATYPE= 'REAL*4  '           / Floating Point                                 
        #   DATE    = '2019-06-12'         / Creation UTC (CCCC-MM-DD) date of FITS header  
        #   COMMENT FITS (Flexible Image Transport System) format is defined in 'Astronomy  
        #   COMMENT and Astrophysics', volume 376, page 359; bibcode 2001A&A...376..359H    
        #   FIELDNAM= '00015+01050Z'       / FIRST coadded image name                       
        #   OBJECT  = 'J000140+010531'     / Object name (Ehhmmss+ddmmss, E=J or B)         
        #   OBJCTRA = '00 01 40.213'       / Object Right Ascension (J2000)                 
        #   OBJCTDEC= '+01 05 31.08'       / Object Declination                             
        #   OBJCTX  =              687.353 / Object X in image (pixels)                     
        #   OBJCTY  =              591.266 / Object Y in image                              
        #   XCORN   =                  637 / X corner in image (pixels)                     
        #   YCORN   =                  541 / Y corner in image                              
        #   TELESCOP= 'VLA     '           /                                                
        #   INSTRUME= 'VLA     '           /                                                
        #   OBSERVER= 'AB628   '           /                                                
        #   DATE-OBS= '19951016'           /                                                
        #   DATE-MAP= '20070828'           /                                                
        #   BUNIT   = 'JY/BEAM '           /                                                
        #   EPOCH   =              2000.00 /                                                
        #   EQUINOX =              2000.00 /                                                
        #   DATAMAX =            0.0306268 /                                                
        #   DATAMIN =         -0.000723907 /                                                
        #   CTYPE1  = 'RA---SIN'           /                                                
        #   CUNIT1  = 'deg     '                                                            
        #   CRVAL1  =       0.374722289000 /                                                
        #   CDELT1  =         -0.000500000 /                                                
        #   CRPIX1  =              137.000 /                                                
        #   CROTA1  =              0.00000 /                                                
        #   CTYPE2  = 'DEC--SIN'           /                                                
        #   CUNIT2  = 'deg     '                                                            
        #   CRVAL2  =        1.08333333333 /                                                
        #   CDELT2  =          0.000500000 /                                                
        #   CRPIX2  =              34.0000 /                                                
        #   CROTA2  =              0.00000 /                                                
        #   CTYPE3  = 'FREQ    '           /                                                
        #   CRVAL3  =        1435100000.00 /                                                
        #   CDELT3  =          2.18750E+07 /                                                
        #   CRPIX3  =              1.00000 /                                                
        #   CROTA3  =              0.00000 /                                                
        #   CTYPE4  = 'STOKES  '           /                                                
        #   CRVAL4  =        1.00000000000 /                                                
        #   CDELT4  =              1.00000 /                                                
        #   CRPIX4  =              1.00000 /                                                
        #   CROTA4  =              0.00000 /                                                
        #   BMAJ    =           1.7778E-03 / Beam major axis (deg)                          
        #   BMIN    =           1.5000E-03 / Beam minor axis (deg)                          
        #   BPA     =                 0.00 / Beam position angle (deg)                      
        #   SURVEY  = 'FIRST   '           / Survey image obtained from                     
        #   BAND    = 'na      '                                                            
        #   HISTORY AIPS   CLEAN BMAJ=  1.7778E-03 BMIN=  1.5000E-03 BPA=   0.00            
        #   HISTORY AIPS   CLEAN NITER=        3 PRODUCT=1   / NORMAL                       
        #   HISTORY AIPS   IMNAME='00015+01050Z' IMCLASS='COADD' IMSEQ=   1     /           
        #   HISTORY AIPS   USERNO=  111            /                                        
        #   HISTORY                                                                         
        #   HISTORY The FIRST VLA Survey: Copyright 1997 University of California           
        #   HISTORY Permission is granted for publication and reproduction of this          
        #   HISTORY material for scholarly, educational, and private non-commercial         
        #   HISTORY use.  Inquiries for potential commercial uses should be                 
        #   HISTORY addressed to:                                                           
        #   HISTORY                            Robert Becker                                
        #   HISTORY                            Physics Dept                                 
        #   HISTORY                            University of California                     
        #   HISTORY                            Davis, CA  95616                             
        #   HISTORY                                                                         
        #   HISTORY -----------------------------------------------------------             
        #   HISTORY                        The VLA FIRST Survey                             
        #   HISTORY                                                                         
        #   HISTORY          R.H. Becker            R.L. White         D.J. Helfand         
        #   HISTORY     UC Davis and IGPP, LLNL       STScI         Columbia University     
        #   HISTORY                                                                         
        #   HISTORY The FIRST survey (Faint Images of the Radio Sky at Twenty-cm) is a      
        #   HISTORY high angular resolution survey of the North Galactic Cap intended to    
        #   HISTORY match the 10,000 deg**2 survey area of the Sloan Digital Sky Survey     
        #   HISTORY (Gunn and Knapp 1992 PASP 43, 267).  Observations are made in the       
        #   HISTORY VLA's B configuration operating at frequencies of 1365 and 1435 MHz in  
        #   HISTORY the bandwidth synthesis mode. The synthesized beam has a FWHM of 5.4".  
        #   HISTORY Each field is observed for 3 minutes, resulting in a typical image rms  
        #   HISTORY of <0.15 mJy/beam.  The raw UV data from the observations are archived  
        #   HISTORY at the VLA and are available on the day they are taken without the      
        #   HISTORY standard proprietary period.  The fully self-calibrated UV datasets     
        #   HISTORY for each pointing are also archived at the NRAO.                        
        #   HISTORY                                                                         
        #   HISTORY Source positions derived from FIRST images have an error of 0.5" rms    
        #   HISTORY for the weakest discernible point sources (~0.75 mJy) and               
        #   HISTORY substantially smaller uncertainties for brighter sources. The source    
        #   HISTORY surface density is ~100 per square degree. The grid of images on the    
        #   HISTORY sky is sufficiently dense that considerable improvement in sensitivity  
        #   HISTORY can be obtained by coadding adjacent fields.  These coadded images      
        #   HISTORY are named HHMMM+-DDMMM (with implied decimal points before the final    
        #   HISTORY M's) corresponding to the J2000 coordinates of the field center to the  
        #   HISTORY nearest 6"; e.g., the field centered at 07 42 00 +30 45 36 is named     
        #   HISTORY 07420+30456. The images have 1.8" x 1.8" pixels, and vary in size       
        #   HISTORY depending on declination.                                               
        #   HISTORY                                                                         
        #   HISTORY Details of the survey design, pipeline processing algorithms, and       
        #   HISTORY catalog generation may be found in Becker et al. 1995 and White et al.  
        #   HISTORY 1997; these papers, along with further information on FIRST, are        
        #   HISTORY retrievable via the FIRST WWW home page: http://sundog.stsci.edu/ .     
        #   HISTORY Users interested in deriving quantitative information from the maps for 
        #   HISTORY scientific purposes are advised to read the FIRST papers first.         
        #   HISTORY Questions or comments on the FIRST survey should be addressed to        
        #   HISTORY Bob Becker at bob@igpp.ucllnl.org or 510-423-0664.                      
        #   
        #   In [167]: 
        #
        img = fits.PrimaryHDU(stamp.data, header=hdu[0].header)
    
        # writing to a pretend file in memory
        mem_file = io.BytesIO()
        img.writeto(mem_file)
        return img


    def get_cutout(self,position, size):
        try:
            tiles  = self.get_tiles(position,size)
            tile   = self.paste_tiles(tiles,position,size)
            cutout = self.trim_tile(tile,position,size)
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

