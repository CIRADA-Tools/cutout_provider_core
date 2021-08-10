import io
import numpy as np
from astropy import units as u
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import mad_std
from astropy.visualization import imshow_norm, ZScaleInterval, MinMaxInterval, AsinhStretch, ManualInterval, ImageNormalize
from .statistics import rms_mad, error_median, robust_stats_radio


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
    # output.seek(0)
    return output


def asinh_soften_for_noise_RMS(image_data_2D, factor=3, interval=MinMaxInterval()):
    (image_MIN, image_MAX, image_MEDIAN,
     image_RMS, image_STD, image_MEAN) = robust_stats_radio(image_data_2D)
    (image_MIN, image_MAX) = interval.get_limits(image_data_2D)
    if factor*image_RMS >= image_MIN:
        parameter = (factor*image_RMS-image_MIN)/(image_MAX-image_MIN)
    else:
        parameter = 0.1 # Default for asinh
    return parameter

# main, calls other methods to create and return thumbnail in image buffer
def get_thumbnail(hdu, survey):
    wcs = WCS(hdu.header)
    # # squeeze to 2d from nd to make thumbnail
    image_data = np.squeeze(hdu.data) #WILL NEED A BETTER TRIMMER FOR CUBES LIKELY trim_axes(hdu, wcs)
    img_buffer = asinh_plot(wcs, image_data, survey)
    return img_buffer.getvalue()


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
