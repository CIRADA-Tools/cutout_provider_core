import io
import numpy as np
from astropy import units as u
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import mad_std
from astropy.visualization import imshow_norm, ZScaleInterval, MinMaxInterval, AsinhStretch, ManualInterval, ImageNormalize


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

def overall_variability_t_stat(Fluxs, errs):
    # still deciding how to do this STOPGAP SOLUTION
    wm_flux = variance_weighted_mean(Fluxs, errs)
    chi_sq = np.sum(((np.array(Fluxs)-wm_flux)/np.array(errs))**2)
    return chi_sq

def overall_modulation_index(Fluxs, errs):
    fnsx = flux_nxs(Fluxs, errs)
    if fnsx<0:
        print("stats error!! fnsx<0 ", str(e))
        return -1.0
    f_var = np.sqrt(fnsx)
    return np.abs(2*(f_var-1)/(f_var+1))

# called with 2 or more fluxs
# ref eq.(1) from https://arxiv.org/pdf/1601.01693.pdf
def variability_t_stat(Fluxs, errs):
    if len(Fluxs)>2:
        return overall_variability_t_stat(Fluxs, errs)
    try:
        return np.abs((Fluxs[1]-Fluxs[0])/np.sqrt(errs[0]*2 + errs[1]*2))
    except Exception as e:
        print("stats error!! variability_t_stat ", str(e))
        return -1

# called with 2 or more fluxs
# ref eq.(2) from https://arxiv.org/pdf/1601.01693.pdf
def modulation_index(Fluxs, errs):
    if len(Fluxs)>2:
        return overall_modulation_index(Fluxs, errs)
    try:
        return np.abs((Fluxs[1]-Fluxs[0])/variance_weighted_mean(Fluxs, errs))
    except Exception as e:
        print("stats error!! modulation index ", str(e))
        return -1

def mse(errs):
    # mean_square_error
    return sum([e**2 for e in errs])/len(errs)

#using list of Fluxs
# Flux Excess normalized variance
# NormalizedExcessVariance ref: https://arxiv.org/pdf/astro-ph/0307420.pdf before eq 10
def flux_nxs(Fluxs, errs):
    try:
        mean_flux = np.mean(Fluxs)
        sample_variance = sum([(x-mean_flux)**2 for x in Fluxs])/(len(Fluxs)-1)
        xs = sample_variance - mse(errs)
        return xs/(mean_flux**2)
    except Exception as e:
        print(str(e))
        return -1 # return dummy variable

# ref: https://arxiv.org/pdf/astro-ph/0307420.pdf before eq 11 but use eq. (B1) for cases! and
# use the first case for S^2 < 3 * mean(sigma_err^2) and the second for S^2 >= 3 * mean(sigma_err^2)
def error_flux_nxs(nxs, Fluxs, errs):
    if nxs == -1: # if nxs raised error then store -1 here too
        return -1
    f_mse = mse(errs)
    mean_flux = np.mean(Fluxs)
    sample_variance = nxs*(mean_flux**2) + f_mse
    if (sample_variance < (3*f_mse)):
        return np.sqrt(2/len(Fluxs)) * (f_mse/(mean_flux**2)) # (eq B1 first case)
    else:
        f_rms_var = np.sqrt(nxs) # fractional root mean square (rms) variability amplitude
        return np.sqrt(f_mse/len(Fluxs))*(2*f_rms_var/mean_flux) # (eq B1 2nd case)

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
    # returns in mJy/beam
    # Flatten array
    #image_data_1D = image_data_2D.flatten() #potentially unneccessary
    # remove non finite values
    image_data_2D = image_data_2D[np.where(np.isfinite(image_data_2D))]
    image_RMS     = rms_mad(image_data_2D)*1000 # mJy/beam
    image_MEDIAN  = np.median(image_data_2D)*1000
    image_STD     = mad_std(image_data_2D)*1000
    image_MIN     = np.min(image_data_2D)*1000
    image_MAX     = np.max(image_data_2D)*1000
    image_MEAN    = np.mean(image_data_2D)*1000
    # Alternative for looking at RMS of radio data
    # Do not use at this point
    #neg_data = image_data_1D[np.where(image_data_1D <0)]
    #neg_data =np.concatenate((neg_data,-neg_data), axis=None)
    #image_RMS  = biweight_scale(neg_data,M=0)
    return(image_MIN, image_MAX, image_MEDIAN, image_RMS, image_STD, image_MEAN)
