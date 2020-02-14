## This Cutout Core repository based on common utilities and functionality forked from http://orbit.dunlap.utoronto.ca/michelle.boyce/Continuum_common

### Applications currently using this common Core Cutout code include:    
- a Command Line Interface with instructions included below      
- a public web service found at http://206.12.91.186/  with project code hosted at http://orbit.dunlap.utoronto.ca/falon3/cirada_cutouts

       
           
# Cutout Gathering

Grabbing cutouts from various surveys

**The code currently is in a testing phase**, and at the moment grabs ~140 Best & Heckman galaxies as FITS files from all included surveys.

### Currently included surveys:

| Survey | Band |
|---|---|
| FIRST | Radio |
| NVSS | Radio|
| VLASS| Radio|
|WISE|Infrared|
| SDSS| Optical|
| PanSTARRS| Optical|

### Surveys to include:

| Survey | Band|
|---|---|
|LOTSS|Radio|
|Akari|Infrared|
|BOSS|Optical|
|GALEX|Ultraviolet|
|3XMMRLQSO|X-Ray|
|RASS|X-Ray|
|FERMI (LAT)|Gamma-Ray|
|HESS|Gamma-Ray|

---

### How to run
from the command line:    
`$python3 fetch_cutouts.py`    
with Commands:    
  `fetch        Single cutout fetching command.   `     
  `fetch_batch  Batch cutout fetching command.   `     
  
Options:    
      -c, --coords TEXT     [required one of either -c or -n]    
      -n, --name TEXT
      -f, --file TEXT       batch file(s) name(s)  [required if fetch_batch]   
      -r, --radius INTEGER  [required]     
      -s, --surveys TEXT   
      -o, --output TEXT   
      -g, --groupby TEXT   
      -cf, --config TEXT   [optional]    
      --overwrite           overwrite existing duplicate target files (default
                            True)   
      --flush               flush existing target files (supersedes --overwrite)   
      --help                Show this message and exit.   
      
Argument Descriptions:    
`-c 'coords' for Source coordinates OR`    
`-n 'name' for Source name`    
     
      example accepted coordinate formats:    
      > RA,DEC or 'RA, DEC' in degrees    
      > '00h42m30s', '+41d12m00s' or 00h42m30s,+41d12m00s    
      > '00 42 30 +41 12 00'    
      > '00:42.5 +41:12'    
      if name:    
      > The name of the object to get coordinates for, e.g. 'M42'    
      
`-r 'radius' is the Integer search radius around the specified source location in arcmin.`    
      The cutouts will be of maximum width and height of 2*radius    

`-s 'surveys' is one or several surveys comma separated without spaces between.`    
      Implemented surveys include: FIRST,VLASS,WISE,SDSS,PANSTARRS,NVSS       
           
      Filters for each survey may be specified in the following formats:    
      > "WISE(w2),SDSS[g,r]"    
      > "WISE[w1],VLASS"    
      > WISE,VLASS    
      
      If no filters are specified then the default filter is used for each.    
      If surveys argument is not specified then will fetch from ALL implemented    
      surveys with default filters for each survey.    
    
`-o 'output' is the directory location to save output FITS images to.`    
      Output will be furthered separated into subfolders for the corresponding survey.    
      Default location is a folder named 'data_out/' in this current directory.    

`-g 'groupby' is an option to separate FITS results by "MOSAIC", "DATE-OBS", or "NONE" (default).`     
    
      > "MOSAIC": if the requested position and radius straddle boundaries in multiple      
                  FITS images for a given survey a mosaicked FITS file will be generated    
                  from all of these input images with each input image as an extension of    
                  the corresponding mosaicked FITS. Mosaics are largely provided for visual    
                  use only.    
      > "DATE-OBS": For surveys VLASS, FIRST, NVSS, or PanSTARRS a Mosaicked FITS is made    
                  (when needed) for every unique DATE-OBS.     
      > "NONE" (default): All resulting FITS images in the requested survey are returned    
                  without doing any mosaicking    
    
`-cf 'config' is to specify a YAML config file for settings, ex."config.yml".`    
      *Note: Specified command line args will overwrite these settings.`          
   
`-f "file" FOR FETCH_BATCH ONLY. The CSV file(s) name. `      

       CSV must at least have separate columns named "RA" and "Dec"    
       (or any of the variants below, but there can only be one variant of    
       RA and one of Dec per file). A column labelled "Name" or "NAME" may also be used.   
       For a given source, coordinates will be evaluated via "RA" and "Dec" if   
       they are non-empty. If a line does not have a valid coordinate position,   
       but does have a "Name" column value, the service will attempt to resolve   
       the source name.   
           
       Accepted variants of RA and Dec Column header names are:    
       R.A.   
       Right Ascension   
       RA (J2000)   
       R.A. (J2000)   
       Right Ascension (J2000)   
       RAJ2000   
       DEC   
       DEC.   
       Declination   
       DEC (J2000)   
       DEC. (J2000)   
       Declination (J2000)   
       DecJ2000   
         
       Source names will be resolved via the Sesame Name Resolver:    
       http://vizier.u-strasbg.fr/viz-bin/Sesame    
        
Sample command looks like:    
`python3 fetch_cutouts.py fetch -n M87 -s VLASS,WISE -r 3 -g MOSAIC`    
   
This will fill `data_out` with the FITS files separated by Survey name directory.    

**FITS Processing to JPEG**

`python3 process_fits.py`

This takes whatever supported FITS files are in `out` and processes them into `.jpg`'s.
_All included surveys (see above) are supported_

### Installation Notes

You'll need to install <a target=_blank href="https://montage-wrapper.readthedocs.io/en/v0.9.5/#installation">Montage</a>, among other standard fair, which can be a littly tricky. Here's the Coles Notes:

 * Download `http://montage.ipac.caltech.edu/download/Montage_v5.0.tar.gz` from <a target=_blank href="http://montage.ipac.caltech.edu/docs/download2.html">`http://montage.ipac.caltech.edu/docs/download2.html`</a>.
 * `tar xvzf Montage_v5.0.tar.gz`
 * Got into the `Montage_v5.0` and type `make` to build.
 * Move the `Montage_v5.0` to `~/.montage/Montage_v5.0`, say, and add `~/.montage/Montage_v5.0/bin` to `$PATH`.
 * To test, run `mAdd` and you should see something like,<br>```[struct stat="ERROR", msg="Usage: mAdd [-d level] [-p imgdir] [-n(o-areas)] [-a mean|median|count] [-e(xact-size)] [-s statusfile] images.tbl template.hdr out.fits"]```<br>indicating it is installed correctly.


## NOTE: if you get pyopenssl ssl.SSLError('bad handshake....
Check which version of requests you are using.    
`import requests`    
`print(requests.__version__)`    
Try downgrading to version 2.11.1. It worked as seen from https://stackoverflow.com/questions/40741361/python-requests-gives-me-bad-handshake-error/40741362    

`pip3 uninstall requests`    
`pip3 install requests==2.11.1`     
