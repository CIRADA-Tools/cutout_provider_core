import io

import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.visualization import imshow_norm, MinMaxInterval, AsinhStretch


def asinh_plot_VLASS_mJy(wcs, image_data, showgrid=False):
    image_data *= 1000 # convert to mJy
    cmap = plt.get_cmap('gray')
    plt.rcParams.update({'font.size': 24})
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
    output = io.BytesIO()
    plt.savefig(output, bbox_inches="tight")
    return output


def trim_axes(hdu, wcs):
    # trim to 2d from nd to make thumbnail
    naxis = wcs.naxis
    # TODO: this needs to be checked for working with fits cubes (as per greg)
    while naxis > 2:
        w = wcs.dropaxis(2)
        naxis -= 1
    image_data = np.squeeze(hdu.data)
    return image_data



def get_thumbnail(hdu):
    wcs = WCS(hdu.header)
    # # trim to 2d from nd to make thumbnail
    image_data = trim_axes(hdu, wcs)
    img_buffer = asinh_plot_VLASS_mJy(wcs, image_data)
    return img_buffer.getvalue()
