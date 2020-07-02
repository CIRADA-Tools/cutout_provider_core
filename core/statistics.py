import io
import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import mad_std
from astropy.visualization import imshow_norm, ZScaleInterval, MinMaxInterval, AsinhStretch, ManualInterval, ImageNormalize




def mse(errs):
    # mean_square_error
    return sum([e**2 for e in errs])/len(errs)

#using list of Fluxs
# Flux Excess normalized variance
# NormalizedExcessVariance ref: https://arxiv.org/pdf/astro-ph/0307420.pdf before eq 10
def flux_nxs(Fluxs, errs):
    mean_flux = np.mean(Fluxs)
    sample_variance = sum([(x-mean_flux)**2 for x in Fluxs])/(len(Fluxs)-1)
    xs = sample_variance - mse(errs)
    return xs/(mean_flux**2)

# ref: https://arxiv.org/pdf/astro-ph/0307420.pdf before eq 11
def error_flux_nxs(flux_nxs, Fluxs, errs):
    f_rms_var = np.sqrt(flux_nxs) # fractional root mean square (rms) variability amplitude
    mean_flux = np.mean(Fluxs)
    f_mse = mse(errs)
    return np.sqrt((np.sqrt(2/len(Fluxs))*(f_mse/(mean_flux**2)))**2 + (np.sqrt(f_mse/len(Fluxs))*(2*f_rms_var/mean_flux))**2)

#calculates the weigted mean using inverse square variances as weights
def variance_weighted_mean(items, item_variances):
    #if variance is 0 make it 0.000001 no divide by zero!
    # test using
    # items= [1.2, 1.3, 1.4, 1.2, 1.3, 1.4]
    # item_variances = [0.0025, 0.0225, 0.01, 0.0625, 0.5625, 0.0006]
    # variance_weighted_mean = 1.389068165073782
    return sum([a/(max(e,0.000001)**2) for a,e in zip(items, item_variances)])/sum([1/(max(e,0.000001)**2) for e in item_variances])

def error_variance_weighted_mean(variances):
    # test:
    # variances = [0.0025, 0.0225, 0.01, 0.0625, 0.5625, 0.0006]
    # error variance weighted mean = 0.0005822213011960972
    #if variance is 0 make it 0.000001 no divide by zero!
    return np.sqrt(1/sum([1/(max(e,0.000001)**2) for e in variances]))

def rms_mad(data):
  # Calculates the standard deviation of data using
  # the median absolute deviation, but forcing the median to be 0
  # i.e., Calculates an RMS
  return 1.482602218505602*np.median(np.absolute(data))

def error_median(all_data):
    # calculates the error of the median of the data
    return rms_mad(all_data)/np.sqrt(len(all_data)-1)

def robust_stats_radio(image_data_2D):
    # Calculates robust statistics useful for images
    # Flatten array
    #image_data_1D = image_data_2D.flatten() #potentially unneccessary
    # remove non finite values
    image_data_2D = image_data_2D[np.where(np.isfinite(image_data_2D))]
    image_RMS     = rms_mad(image_data_2D)
    image_MEDIAN  = np.median(image_data_2D)
    image_STD     = mad_std(image_data_2D)
    image_MIN     = np.min(image_data_2D)
    image_MAX     = np.max(image_data_2D)
    image_MEAN    = np.mean(image_data_2D)
    # Alternative for looking at RMS of radio data
    # Do not use at this point
    #neg_data = image_data_1D[np.where(image_data_1D <0)]
    #neg_data =np.concatenate((neg_data,-neg_data), axis=None)
    #image_RMS  = biweight_scale(neg_data,M=0)
    return(image_MIN, image_MAX, image_MEDIAN, image_RMS, image_STD, image_MEAN)
