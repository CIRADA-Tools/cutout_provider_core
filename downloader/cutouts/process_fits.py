import os
import errno

import numpy as np

from astropy.io import fits
from astropy.visualization import make_lupton_rgb

import matplotlib.pyplot as plt

'''
Grab FITS files from /out, process into JPGs and save to /processed

TODO:
Should probably be refactored into the Survey objects
'''


def vlass_stretch(data):

    vmin = -0.000235248
    vmax = 0.00798479

    return arcsinh_stretch(data, vmin, vmax)


def nvss_stretch(data):

    vmin = -0.00070356
    vmax = 0.162288

    return arcsinh_stretch(data, vmin, vmax)


def first_stretch(data):

    vmin = -0.00551651
    vmax = 0.0634952

    return arcsinh_stretch(data, vmin, vmax)


def sdss_stretch(data):

    vmin = -0.05538
    vmax = 20.67399

    return lupton_stretch(data, vmin, vmax)

def panstarrs_stretch(data):

    vmin = -250.226
    vmax = 80822.4

    fit = np.interp(data, (vmin, vmax), (0, 10))

    fit = make_lupton_rgb(fit[2], fit[1], fit[0], stretch=0.7, Q=11)

    return fit


# make the directory structure if it doesn't exist
def make_dir(dirname):

    try:
        os.makedirs(dirname)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def lupton_stretch(data, vmin, vmax):

    fit = np.interp(data, (vmin, vmax), (0, 10))

    fit = make_lupton_rgb(fit[2], fit[1], fit[0], stretch=0.5, Q=9)

    return fit


def arcsinh_stretch(data, vmin, vmax):

    fit = np.interp(data, (vmin, vmax), (0, 1))

    fit = np.arcsinh(10*fit)/3

    return fit


def get_fits_data(in_file):
    f = fits.getdata(in_file)

    # trim any 1-long axes
    f = np.squeeze(f)
    
    return f


def save_image(data, out_file):
    
    # a color image
    if data.shape == 3:
        plt.imsave(fname=out_file, arr=data)
    else:
        plt.imsave(fname=out_file, arr=data, cmap='gray')
    
    
def process_fits(in_file, process_func, out_file):

    d = get_fits_data(in_file)

    processed = process_func(d)

    save_image(processed, out_file)


def format_infile(directory, filename):
    return '{0}/{1}'.format(directory, filename)


def format_outfile(directory, filename):
    return '{0}/{1}.jpg'.format(directory, filename[:-4])


def survey_filter(survey):

    def f(filename):

        return survey.lower() in filename.lower()

    return f


def survey_list(files, in_dir, out_dir, process):

    files = list(files)

    infiles = (format_infile(in_dir, f) for f in files)
    outfiles = (format_outfile(out_dir, f) for f in files)

    return ((i, process, o) for (i, o) in zip(*(infiles, outfiles)))


def main():

    in_dir = 'data/out'
    out_dir = 'data/processed'

    make_dir(out_dir)

    (_, _, files) = next(os.walk(in_dir))
    fits_files = [f for f in files if 'fits' in f[-4:].lower()]

    vlass_files = filter(survey_filter('vlass'), fits_files)
    nvss_files = filter(survey_filter('nvss'), fits_files)
    first_files = filter(survey_filter('first'), fits_files)
    sdss_files = filter(survey_filter('sdss'), fits_files)
    panstarrs_files = filter(survey_filter('panstarrs'), fits_files)

    survey_processes = [
                        (vlass_files, vlass_stretch),
                        (nvss_files, nvss_stretch),
                        (first_files, first_stretch),
                        (sdss_files, sdss_stretch),
                        (panstarrs_files, panstarrs_stretch)
                        ]

    processes = (survey_list(s, in_dir, out_dir, p) for (s, p) in survey_processes)

    flatten = lambda l: [item for sublist in l for item in sublist]

    for image in flatten(processes):

        print(image[0])

        process_fits(*image)


if __name__ == '__main__':
    main()
