### Produces 132x132px, 3x3' image cutouts for FIRST, WISE (W1), and VLASS imaging,
### products should be useable with ClaRAN AI source finding code

### A future version (v2) will grab the cutouts direct from the web, incorporating
### the astro image cutout tool developed by Michael Ramsey (ramseym@myumanitoba.ca)
#### Currently, and just for this test sample while the code is being developed,
#### the raw image files are stored locally (having been grabbed manually by the author)\

### This code has been principally developed by Yjan Gordon (yjan.gordon@umanitoba.ca)
### as part of the CIRADA project (www.cirada.org) in March 2019
#### WISE (W1) cutout images require no post-processing
#### FIRST cutouts need regridding from 100x100px to 132x132px
#### VLASS QL images need trimming 3x3' centred on the target
#### VLASS trimmed images then need regridding to 132x132px from 180x180px
##### potential edge sources in VLASS QL imaging undeline the need to incorporate
##### Michael Ramsey's cutout grabbing code in the next version

####do we need to mask the radio image at all? i.e. only use Signal above n-sigma?

#################################################################################
#################################################################################
### import modules

import numpy as np, pandas as pd, matplotlib.pyplot as plt, astropy.units as u, os
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import wcs
from FITS_tools.hcongrid import hcongrid


########################################
#time of code check - not important for code to run
import time
start_time = time.time()
########################################

#################################################################################
#################################################################################
###define parameters/variables

Testing = False ###is this a test for which a single target will suffice?
save_files = False ###run code saving output files

owrite = False ###overwrite files?

###file paths
survey_path = '../../../../../survey_data/'
images_in = '../../../../VLASS/vlass_source_finding/vlass_source_test/image_files/raw_files/'
images_out = 'data/' ####testing output path
tgal_fpath = 'data/' ##for test input targets


###other variables
###number of pixels to regrid to
npx = 132

##the maximum separation of components for a source in a 3x3' box will be
## sqrt(1.5^2 + 1.5^2), add 0.1 to dx/y for leeway
maxsep_am = np.sqrt(2*1.6**2)

###limit vlass to only components with s/n > some value (5?)?
VLASS_SigNoise = True
VLASS_SN_limit = 5 ##minimum s/n limit to use if only using reliable detections


###image plotting
plt.interactive(True)
plot_images = True
colmap = 'Blues_r' ### colour map to use
npx = 132 ##n pixels for regridded image (x=y=npx)


#################################################################################
#################################################################################
###define functions to use

def nrao_coords(cat):
    ### convert astropy skycoords into simple string for use with NRAO interface
    ra = cat.ra.hms
    dec = cat.dec.dms
    
    ###create hhmmss for ra
    hh = np.array(ra.h, dtype=int)
    rm = np.array(ra.m, dtype=int)
    rs = np.array(np.round(ra.s, 2), dtype=float)
    
    ##comnine to single string
    hms = [str(hh[i]) + ' ' + str(rm[i]) + ' ' + str(rs[i]) for i in range(len(hh))]
    
    ###repeat for ddmmss, dec
    dd = np.array(abs(dec.d), dtype=int)
    dm = np.array(abs(dec.m), dtype=int)
    ds = np.array(np.round(abs(dec.s), 2), dtype=float)
    
    ###mask for negative dec
    neg = (dec.d<0) | (dec.m<0) | (dec.s<0)
    negs = np.ones(len(neg), dtype=str)
    negs[neg==True] = '-'
    negs[neg==False] = '+'
    
    ###combine to single string
    dms = [negs[i]+str(dd[i]) + ' ' + str(dm[i]) + ' ' + str(ds[i]) for i in range(len(dd))]
    
    ###combine hms+/-dms
    sexco = [hms[i] + ' ' + dms[i] for i in range(len(hms))]
    return(sexco)


def targetname(ra, dec):
    ###get a target name from ra/dec, needs to generecise for operating on coordinate arrays
    coords = SkyCoord(ra=ra, dec=dec, unit='deg')
    hms = coords.ra.hms
    dms = coords.dec.dms
    
    ###define separator of ra/dec
    if coords.dec<0:
        pm = '-'
    else:
        pm = '+'
    
    ###list strings of integer h/d, m, s
    hms_s = [str(int(hms.h)), str(int(hms.m)), str(int(hms.s))]
    dms_s = [str(int(abs(dms.d))), str(int(abs(dms.m))), str(int(abs(dms.s)))]
    
    ##ensure that strings are of len == 2, i.e. 3->03
    for i in range(len(hms_s)):
        if len(hms_s[i]) < 2:
            hms_s[i] = '0'+hms_s[i]
        if len(dms_s[i]) < 2:
            dms_s[i] = '0'+dms_s[i]
    
    ###join elements into single 6digit ra/decs
    a, d = ''.join(hms_s), ''.join(dms_s)

    ## combine into single string target name
    tname = 'J' + a + pm + d
    return(tname)


def arcsinh_scale(data):
    ###scale image array using arcsinh scale of SAO DS9
    ###  y = arcsinh(10x)/3   for 0 < x < 1
    ###determine min/max values
    imin, imax = np.min(data), np.max(data)
    
    ###rescale so that max value == 10
    sfact = 10/imax
    data = sfact*data
    
    ###use arcsineh stretch and scale result by 1/3
    data = np.arcsinh(data)/3
    
    return(data)


def image_cutout(data, colourmap=colmap):
    ##create figure without a frame
    fig = plt.figure(figsize = (3,3), frameon=False)
    ax = fig.add_axes([0., 0., 1., 1.])
    ax.axis('off')
    
    ###add image
    ax.imshow(data, origin='lower', cmap=colourmap)
    
    ###set physical size of cutout
    fig.set_size_inches(2.56, 2.56)
    
    return(fig)


#################################################################################
#################################################################################
###import catalogue data

##test sample of known radio galaxies - taken from Gordon+19 (submitted)
tgals = pd.read_csv(tgal_fpath+'targs_for_michelle.csv')

###columns required (LERG: RA, DEC, merger status - use most_votes)
tcols = ['RA', 'DEC']
targets = tgals[tcols]


###vlass and first catalgues - just vlass catalogue required
###may need first in future version not using test input catalogue
###setup for my own personal fits file, but if using the NRAO catalogue use pandas
### i.e. vlass =  pd.read_csv(survey_path+'vlass_catalogue_filename.csv')
### if RA in decimal hours can either convert to degrees using astropy.coordinates.Angle
### or specify ra=ra*u.hour, dec=dec.deg when setting up vcat.
vlass = Table(fits.open(survey_path+'VLASS/test/'
                        +'vlass11_ql_catalogue_full_YGcorrectedRA_Mar19.fits')
              [1].data).to_pandas()

###SN limit on vlass
if VLASS_SigNoise == True:
    vlass = vlass[(vlass.Total_flux/vlass.E_Total_flux) > VLASS_SN_limit]

#################################################################################
#################################################################################
###cross match targets with VLASS to find nearest component and image info
#### this is likely not necessary when using Michaels code to obtain cutouts

#### MICHELLE: you should look into calling michaels code here to grab cutouts
#### using RA/DEC from tgals.
#### Having the cutout will also alter the


###create astropy position catalogue for targets
targcat = SkyCoord(ra=np.array(targets['RA']),
                   dec=np.array(targets['DEC']), unit='deg')


######find those with VLASS components
###set up astropy position catalogue for VLASS components:
vcat = SkyCoord(ra=np.array(vlass['RA']),
                dec=np.array(vlass['DEC']), unit='deg')


### x-match targetcat and vcat
tind, tsep = targcat.match_to_catalog_sky(vcat)[:2] ##gets nearest match index and sep


###assign new columns to targets
targets['closest_VLASS'] = tind
targets['sep_VLASS_arcmin'] = tsep.arcmin

###include Tile and Pointing info to get correct image
targets['Tile'] = np.array(vlass['Tile'].iloc[targets.closest_VLASS])
targets['Pointing'] = np.array(vlass['Pointing'].iloc[targets.closest_VLASS])

###create column for QL image file name
ims1 = 'VLASS1.1.ql.' ###first part of file string - doesn't vary with image
ims2 = '.10.2048.v1.I.iter1.image.pbcor.tt0.subim.fits' #end of fname - doesn't vary
imstring = [ims1 + list(targets.Tile)[i] + '.' + list(targets.Pointing)[i]
            + ims2 for i in range(len(targets))] ##creates fname from Tile and Pointing

targets['ql_image'] = imstring


#################################################################################
#################################################################################

#################################################################################
#################################################################################
####image processing

##test on first row only of tatgets, remove once comfortable code works
if Testing == True:
    targets = targets[:1]

tnames = []

###trim and regrid vlass images, test by plotting

for i in range(len(targets)):
    ####for grabbing image from QL file stored locally, below is good
    #### once grabbing cutout using micheals code the image array (ql_data here)
    #### needs to be the cutout instead
    #### does michael's cutout include wcs header info? if not that will need to be adressed
    
    targ = targets.iloc[i]
    ra, dec = targ['RA'], targ['DEC']
    vlass_ql_image_file = 'VLASS_QL/' + targ['ql_image']
    sourcename = targetname(ra=ra, dec=dec)
    tnames.append(sourcename)
    
    ##open QL image - again use array from michael's cutout here instead
    hdu = fits.open(images_in + vlass_ql_image_file)
    ql_head, ql_data = hdu[0].header, hdu[0].data[0][0]


    ##use wcs module
    wcsinfo = wcs.WCS(ql_head)
    x, y = wcsinfo.wcs_world2pix(ra, dec, ql_head['CRVAL3'], ql_head['CRVAL4'], 0)[:2]
    x, y = int(np.round(x)), int(np.round(y)) ###rounds to intergers

    dpix = 90 ##number of px either side of target to use, VLASS QL has 1"/px
    
    ###set up min/max x,y accounting for potential edge of tile (i.e. in range [0,len(image)-1])
    minvals, maxvals = [x-dpix, y-dpix], [x+dpix, y+dpix]
    
    edgecase = False ###assume not an edgecase, will be altered if untrue
    
    for j in range(len(minvals)):
        if minvals[j] < 0:
            minvasl[j] = 0
            edgecase = True
        if maxvals[j] > len(ql_data)-1:
            maxvals[j] = len(ql_data)-1
            edgecase = True
    
    xmin, ymin = minvals
    xmax, ymax = maxvals

    imtrim = ql_data[ymin:ymax, xmin:xmax]


    ### from here if trimhead == michaels cutout array then this should work

    ##regrid to 132x132px
    ##set up header info for trimmed image required for regrid
    trimhead = fits.PrimaryHDU(imtrim).header

    ##new header info
    trimhead['CTYPE1'] = ql_head['CTYPE1']
    trimhead['CRVAL1'] = np.round(ra, 5)
    trimhead['CDELT1'] = ql_head['CDELT1']
    trimhead['CRPIX1'] = dpix
    trimhead['CUNIT1'] = ql_head['CUNIT1']

    trimhead['CTYPE2'] = ql_head['CTYPE2']
    trimhead['CRVAL2'] = np.round(dec, 5)
    trimhead['CDELT2'] = ql_head['CDELT2']
    trimhead['CRPIX2'] = dpix
    trimhead['CUNIT2'] = ql_head['CUNIT2']

    ##set up dummy header info for new grid shape
    nghead = fits.PrimaryHDU(np.zeros(shape=(npx, npx))).header

    ##list of keys to add to new header, reduce to only those needed to add
    keylist = list(trimhead.keys())
    keylist = keylist[keylist.index('CTYPE1'):keylist.index('CUNIT2')]

    ##add to new header
    for key in keylist:
        if 'TYPE' in key or 'VAL' in key or 'UNIT' in key:
            nghead[key] = trimhead[key]
        elif 'DELT' in key:
            nghead[key] = (trimhead['NAXIS1']/npx)*trimhead[key]
        elif 'PIX' in key:
            nghead[key] = int(npx/2)

    ###regrid image array to new shape
    regrid = hcongrid(imtrim, trimhead, nghead)

    ##set up fits format output, header info and hdu list to write to file
    new_hdu = fits.PrimaryHDU(regrid)
    new_hdu.header = nghead
    new_hdu.header['QL_file'] = targ['ql_image']
    new_hdu.header['comment'] = ('3x3 arcmin cutout taken from the VLASS '
                                 + 'quick look image named in QL_file. '
                                 + 'This image has been regrided to 132x132px '
                                 + 'by Yjan Gordon (yjan.gordon@umanitoba.ca) '
                                 + 'as part of the CIRADA project (www.cirada.org)')

    new_hdu = fits.HDUList(new_hdu)

    ###below deals with local FIRST and WISE cutouts
    ###these should be substituted for michael's cutouts
    #############################################################################
    ### regrid FIRST image
    firstfile = 'FIRST/' + sourcename + '.fits'
    fhdu = fits.open(images_in + firstfile)
    fhead, fdata = fhdu[0].header, fhdu[0].data

    ##stupid first header is header of main image rather than cutout
    ### set up new header for cutout
    fchead = fits.PrimaryHDU(fdata).header

    ##new columns
    fchead['CTYPE1'] = fhead['CTYPE1']
    fchead['CRVAL1'] = np.round(ra, 5)
    fchead['CDELT1'] = fhead['CDELT1']
    fchead['CRPIX1'] = int(len(fdata)/2)
    fchead['CUNIT1'] = ql_head['CUNIT1'] #grab unit from vlass header...

    fchead['CTYPE2'] = fhead['CTYPE2']
    fchead['CRVAL2'] = np.round(dec, 5)
    fchead['CDELT2'] = fhead['CDELT2']
    fchead['CRPIX2'] = int(len(fdata)/2)
    fchead['CUNIT2'] = ql_head['CUNIT2'] #grab unit from vlass header...


    ###set up new header for regrid
    ngheadf = fits.PrimaryHDU(np.zeros(shape=(npx, npx))).header
    ##add to new header
    for key in keylist:
        if 'TYPE' in key or 'VAL' in key or 'UNIT' in key:
            ngheadf[key] = fchead[key]
        elif 'DELT' in key:
            ngheadf[key] = (fchead['NAXIS1']/npx)*fchead[key]
        elif 'PIX' in key:
            ngheadf[key] = int(npx/2)

    ###regrid image
    rgfirst = hcongrid(fdata, fchead, ngheadf)

    ##set up fits output and add to header
    nf_hdu = fits.PrimaryHDU(rgfirst)
    nf_hdu.header = ngheadf
    nf_hdu.header['FIELD'] = fhead['FIELDNAM']
    nf_hdu.header['comment'] = ('3x3 arcmin FIRST cutout image, from the FIRST '
                                + 'coadded image named in FIELD. '
                                + 'This image has been regrided to 132x132px '
                                + 'by Yjan Gordon (yjan.gordon@umanitoba.ca) '
                                + 'as part of the CIRADA project (www.cirada.org)')

    nf_hdu = fits.HDUList(nf_hdu)

    #############################################################################
    ### adjust WISE header info and output png for adrian (as well as FITS)
    wisefile = 'WISE/' + sourcename + '.fits'
    whdu = fits.open(images_in + wisefile)
    whead, wdata = whdu[0].header, whdu[0].data

    ##stupid wise header is header of main image rather than cutout
    ### set up new header for cutout
    wchead = fits.PrimaryHDU(wdata).header

    ##new columns
    wchead['CTYPE1'] = whead['CTYPE1']
    wchead['CRVAL1'] = np.round(ra, 5)
    wchead['CDELT1'] = whead['CDELT1']
    wchead['CRPIX1'] = int(len(wdata)/2)
    wchead['CUNIT1'] = ql_head['CUNIT1'] #grab unit from vlass header...

    wchead['CTYPE2'] = whead['CTYPE2']
    wchead['CRVAL2'] = np.round(dec, 5)
    wchead['CDELT2'] = whead['CDELT2']
    wchead['CRPIX2'] = int(len(wdata)/2)
    wchead['CUNIT2'] = ql_head['CUNIT2'] #grab unit from vlass header...

    ###no regridding necessary just header tidying
    nw_hdu = fits.PrimaryHDU(wdata)
    nw_hdu.header = wchead
    nw_hdu.header['AT_IM'] = whead['COADDID'] #Atlas image identifier
    nw_hdu.header['comment'] = ('3x3 arcmin WISE W1 cutout image, from the '
                                + 'coadded image named in AT_IM (Atlas identifier).')

    nw_hdu = fits.HDUList(nw_hdu)



    #############################################################################
    ###want to plot images?
    if plot_images == True and save_files==False:
        image_cutout(arcsinh_scale(regrid), colourmap='Blues_r')
        image_cutout(arcsinh_scale(rgfirst), colourmap='Greys_r')
        image_cutout(arcsinh_scale(wdata), colourmap='Reds_r')

    #############################################################################


    ##need to write to file
    if save_files == True:
        
        ###setup wisepng
        wisepng = image_cutout(arcsinh_scale(nw_hdu[0].data), colourmap='YlOrRd_r')
        
        ###set up directory list in images_out
        dlist = os.listdir(images_out)

        ##create a folder for each source to write WISE/FIRST/VLASS files to
        if sourcename not in dlist:
            os.mkdir(images_out+sourcename)

        ##write regridded fits images to file
        ##use format SURVEY+sourcename,
        ##e.g. the vlass image of source J104921-004005 would be VLASS_J104921-004005.fits


        ##VLASS output
        new_hdu.writeto(images_out+sourcename+'/VLASS_'+sourcename+'.fits', overwrite=owrite)

        ##FIRST output
        nf_hdu.writeto(images_out+sourcename+'/FIRST_'+sourcename+'.fits', overwrite=owrite)

        ##WISE output - FITS file and png
        nw_hdu.writeto(images_out+sourcename+'/WISE_'+sourcename+'.fits', overwrite=owrite)
        wisepng.savefig(images_out+sourcename+'/WISE_W1_'+sourcename+'.png', dpi=100)



#############################################################################
targets['name'] = tnames

#####need to incorporate Michaels code - seems to be slight issue with finding x/y pos of coord
#####is this also true for cutouts from michael? if not problem solved, if so need to investigate

#################################################################################
#################################################################################
#################################################################################
#################################################################################
#time taken for code to run
print('END: ', ("--- %s seconds ---" % (time.time() - start_time)))

