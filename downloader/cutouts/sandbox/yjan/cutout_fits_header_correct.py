###correct the fits headers of obtained fits images using michaels CIRADA cutout script

#################################################################################
#################################################################################
### import modules

import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy.time import Time

#################################################################################
#################################################################################
###define parameters/variables


fpath_in = '../../../vlass_source_finding/git/Continuum_common/downloader/cutouts/data/out/'

testimpath = '../../../vlass_source_finding/VLASS_source_test/testgal_J112029+040742/'


#################################################################################
#################################################################################

###def functions

def yrmon_to_isot(yrmon):
    ###converts obs date of form 'yyyymm' to 'yyyy-mm-ddT00:00:00.000'
    ###allows atropy time to be used
    yr = yrmon[:4]
    mon = yrmon[4:6] ###set to assume str len == 6 accounts for some FIRST FITS adding dd at end
    ##can make more complex later to use dd info, at moment good enough as yyyymm will assume dd==15 T==00:00:00 for MJD
    isottime = yr+'-'+mon+'-15T00:00:00.000'
    return(isottime)


def headwrite(data, head, survey, a, d):
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
    elif survey == 'VLASS':
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
        newhead['MJD'] = (Time(yrmon_to_isot(newhead['DATE-OBS'])).mjd,
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
    return(newhead)


#################################################################################
#################################################################################

####simple testing of function -- all below here is testing
#### michelle you just need the functions above (yr_mon_to_isot & headwrite)

####cutout images with headers as currently obtained by fetch_cutouts.py
hdu_ps = fits.open(fpath_in+'J112029+040742_s3arcmin_PanSTARRS.fits')
hdu_vlass = fits.open(fpath_in+'J112029+040742_s3arcmin_VLASS.fits')
hdu_wise = fits.open(fpath_in+'J112029+040742_s3arcmin_WISE.fits')
hdu_first = fits.open(fpath_in+'J112029+040742_s3arcmin_FIRST.fits')


####test cutouts manually downloaded from surveys - original survey headers
testf = fits.open(testimpath+'first_cutout.fits')
testw = fits.open(testimpath+'wise_w1_cutout.fits')
testp = fits.open(testimpath+'panstarrs_iband_cutout.fits')
testv = fits.open(testimpath+'vlass_ql.fits')

targra = 170.1218
targdec = 4.1284



####headers from original survey
hp, hv, hw, hf = hdu_ps[0].header, hdu_vlass[0].header, hdu_wise[0].header, hdu_first[0].header


####new headers for the cutouts using the header info from the original survey files
nhf = headwrite(data=hdu_first[0].data, head=testf[0].header, survey='FIRST',
                a=targra, d=targdec)
nhw = headwrite(data=hdu_wise[0].data, head=testw[0].header, survey='WISE',
                a=targra, d=targdec)
nhp = headwrite(data=hdu_ps[0].data, head=testp[0].header, survey='PanSTARRS',
                a=targra, d=targdec)
nhv = headwrite(data=hdu_vlass[0].data, head=testv[0].header, survey='VLASS',
                a=targra, d=targdec)








