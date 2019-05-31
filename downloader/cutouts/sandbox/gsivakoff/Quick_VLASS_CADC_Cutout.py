#!/usr/bin/env python
# coding: utf-8

# Quick_VLASS_CADC_Cutout.py
# Authors: Gregory Sivakoff
# Version: 0.1
# Version Date: 2019-05-30

#
# This is a quick python script to do an image cutout from VLASS
# I have not had time to comment this or do a bunch of test cases
# 
# 
# You will need astropy and astroquery installed.
# 
#     pip install --pre astroquery


from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import (imshow_norm,
                                   MinMaxInterval,
                                   AsinhStretch)
from astroquery.cadc import Cadc

import shutil
import tempfile
import urllib.request
import urllib.parse

import matplotlib.pyplot as plt
import numpy as np


# Enter position and radius of cutout below.
# If you would like to save the fits file, uncomment the second savefile
# If you like grids on your image, change showgrid to True

#position = SkyCoord('08h45m07.5s', '+54d18m00s', frame='icrs', unit='deg')
position = SkyCoord(131.28125, 54.3, frame='icrs', unit='deg')
radius = 0.05*u.deg

savefile = ''
#savefile = '/Users/gsivakoff/Desktop/temp.fits'

showgrid = False


def open_VLASS_I_mJy(file):
    hdu = fits.open(file)
    extension = 0
    wcs = WCS(hdu[extension].header)
    image_data = hdu[extension].data[0,0,:,:]*1000
    # places image into mJy
    #I'm not sure why the image axes seem to be reversed from the FITS definition
    
    return wcs, image_data


def asinh_plot_VLASS_mJy(wcs, image_data, showgrid=False):
    cmap = plt.get_cmap('hot')

    plt.rcParams.update({'font.size': 18})
    fig = plt.figure(figsize=(10.5, 8), dpi=75)
    ax = plt.subplot(projection=wcs.celestial)
    im, norm = imshow_norm(image_data, ax, origin='lower',
                       interval=MinMaxInterval(), 
                       stretch=AsinhStretch(),
                       cmap=cmap)
    cbar = fig.colorbar(im, cmap=cmap)
    ax.set_xlabel('RA J2000')
    ax.set_ylabel('Dec J2000')
    cbar.set_label('mJy')
    if showgrid:
        ax.grid(color='white', ls='solid')
    plt.show()


def plot_http_fits(url,showgrid=False,savefile=''):
    with urllib.request.urlopen(url) as response:
        with tempfile.NamedTemporaryFile(delete=True) as tmp_file:
            shutil.copyfileobj(response, tmp_file)

            wcs, image_data = open_VLASS_I_mJy(tmp_file.name)
            asinh_plot_VLASS_mJy(wcs, image_data, showgrid)
            if savefile != '':
                shutil.copy(tmp_file.name, savefile)

            tmp_file.close()


def construct_cadc_url(baseurl, position, radius):
    ICRS_position = position.transform_to('icrs')
    basefile = baseurl.split('pub/')[1].split('?')[0]
    if (basefile[-10:] == 'subim.fits' and basefile[:6] == 'VLASS/'):
        url = ( 'https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/cutout?uri=ad:' 
               + urllib.parse.quote(basefile) 
               + ('&cutout=Circle+ICRS+{}+{}+{}').format(ICRS_position.ra.degree,
                                                         ICRS_position.dec.degree,
                                                         radius.to(u.degree).value))
        return url
    else:
        print('CADC URL appears to be incorrect: {}'.format(basefile))
        return None


def get_VLASS_images(position, radius, showgrid=False, savefile=''):

    cadc = Cadc()
    result = cadc.query_region(position, collection='VLASS')
    if len(result) == 0:
        print('No Data at CADC at position: {}'.format(position.to_string('hmsdms')))
    else:
        urls = cadc.get_data_urls(result)
        for url in urls:
            cutout_url = construct_cadc_url(url, position, radius)
            plot_http_fits(cutout_url, showgrid, savefile)


def process_cutout_input(position, radius, showgrid=False, savefile=''):
    radius_deg = radius.to(u.degree).value
    if position.dec.degree < -30:
        print('VLASS only covers Declinations > -30 deg')
    elif (radius_deg >1 or radius_deg < 0.0103):
        print('Radius must be between 0.0103 and 1 deg')
    else:
        get_VLASS_images(position, radius, showgrid, savefile)


process_cutout_input(position, radius, showgrid, savefile)