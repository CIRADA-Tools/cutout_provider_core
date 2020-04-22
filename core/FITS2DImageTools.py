import io
import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import mad_std
from astropy.visualization import imshow_norm, ZScaleInterval, MinMaxInterval, AsinhStretch, ManualInterval, ImageNormalize


def asinh_plot(wcs, image_data, survey, showgrid=False):
    cmap = plt.get_cmap('gray')
    if survey.upper() in ["PANSTARRS", "WISE", "SDSS"]:
        interval = ZScaleInterval()
        cbar_label = ''
        stretch = AsinhStretch()
    else:
        image_data *= 1000 # convert to mJy
        cbar_label = 'mJy'
        interval = MinMaxInterval()
        stretch  = AsinhStretch(asinh_soften_for_noise_RMS(image_data,
                                                           3,
                                                           interval))
    return CIRADA_image_plot(wcs, image_data,
                      cbar_label, cmap,
                      showgrid,
                      interval,
                      stretch)

def CIRADA_image_plot(wcs, image_data,cbar_label, cmap,showgrid=False,
                      interval = MinMaxInterval(),stretch = AsinhStretch()):
    plt.rcParams.update({'font.size': 24})
    fig = plt.figure(figsize=(10.5, 8), dpi=75)
    ax = plt.subplot(projection=wcs.celestial)
    im, norm = imshow_norm(image_data, ax, origin='lower',
                       interval=interval,
                       stretch=stretch,
                       cmap=cmap)
    cbar = fig.colorbar(im, cmap=cmap)
    ax.set_xlabel('RA J2000')
    ax.set_ylabel('Dec J2000')
    cbar.set_label(cbar_label)
    if showgrid:
        ax.grid(color='white', ls='solid')
    output = io.BytesIO()
    plt.savefig(output, bbox_inches="tight")
    plt.close(fig)
    return output

def rms_mad(data):
  # Calculates the standard deviation of data using
  # the median absolute deviation, but forcing the median to be 0
  # i.e., Calculates an RMS
  return 1.482602218505602*np.median(np.absolute(data))

def asinh_soften_for_noise_RMS(image_data_2D, factor=3, interval=MinMaxInterval()):
    (image_MIN, image_MAX, image_MEDIAN,
     image_RMS, image_STD) = robust_stats_radio(image_data_2D)
    (image_MIN, image_MAX) = interval.get_limits(image_data_2D)
    if factor*image_RMS >= image_MIN:
        parameter = (factor*image_RMS-image_MIN)/(image_MAX-image_MIN)
    else:
        parameter = 0.1 # Default for asinh
    return parameter

def robust_stats_radio(image_data_2D):
    # Calculates robust statistics useful for images
    # Flatten array
    #potentially unneccessary
    #image_data_1D = image_data_2D.flatten()
    # remove non finite values
    image_data_2D = image_data_2D[np.where(np.isfinite(image_data_2D))]
    image_RMS    = rms_mad(image_data_2D)
    image_MEDIAN = np.median(image_data_2D)
    image_STD    = mad_std(image_data_2D)
    image_MIN    = np.min(image_data_2D)
    image_MAX    = np.max(image_data_2D)
    # Alternative for looking at RMS of radio data
    # Do not use at this point
    #neg_data = image_data_1D[np.where(image_data_1D <0)]
    #neg_data =np.concatenate((neg_data,-neg_data), axis=None)
    #image_RMS  = biweight_scale(neg_data,M=0)
    return(image_MIN, image_MAX, image_MEDIAN, image_RMS, image_STD)

# def trim_axes(hdu, wcs):
#     # trim to 2d from nd to make thumbnail
#     # naxis = wcs.naxis
#     # # TODO: this needs to be checked for working with fits cubes (as per greg)
#     # this doesn't do anything because wcs isn't returned?......
#     # while naxis > 2:
#     #     w = wcs.dropaxis(2)
#     #     naxis -= 1
#     image_data = np.squeeze(hdu.data)
#     return image_data

# main, calls other methods to create and return thumbnail in image buffer
def get_thumbnail(hdu, survey):
    wcs = WCS(hdu.header)
    # # squeeze to 2d from nd to make thumbnail
    image_data = np.squeeze(hdu.data) #trim_axes(hdu, wcs)
    img_buffer = asinh_plot(wcs, image_data, survey)
    return img_buffer.getvalue()

# def get_image_stats(hdu):
#     wcs = WCS(hdu.header)
#     #TODO check if these stats for only 2D image or okay to put in FITS header
#     # and update FITS header to include stats before saving image
#     image_data_2D = np.squeeze(hdu.data) #trim_axes(hdu, wcs)
#     # what interval to do stats for??
#     interval = MinMaxInterval()
#     (image_MIN, image_MAX) = interval.get_limits(image_data_2D)
#
#     (image_MIN, image_MAX, image_MEDIAN,
#      image_RMS, image_STD) = robust_stats_radio(image_data_2D)
#      return {"image_MIN":image_MIN, "image_MAX":image_MAX, "image_MEDIAN":image_MEDIAN,
#       "image_RMS":image_RMS, "image_STD":image_STD}
