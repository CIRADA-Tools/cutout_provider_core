## This Cutout Core repository based on common utilities and functionality forked from http://orbit.dunlap.utoronto.ca/michelle.boyce/Continuum_common

### Applications currently using this common Core Cutout code include:    
- a Command Line Interface with instructions included below      
- a public web service found at http://206.12.91.186/  with project code hosted at http://orbit.dunlap.utoronto.ca/falon3/cirada_cutouts

### Contents:
+ [Included surveys](#included-surveys)
+ [Installation](#installation)
+ [How to run](#how-to-run)
+ [Output](#output)

# Cutout Gathering

Grabbing cutouts from various surveys

\***Note:** The code currently is in a developing and testing phase

### Included surveys

| Survey | Band |
|---|---|
| VLASS| Radio|
| GLEAM | Radio |
| FIRST | Radio |
| NVSS | Radio|
| WISE |Infrared|
| PanSTARRS| Optical|
| SDSS| Optical|
---

### Installation     
(this was developed on an Ubuntu 18.04 OS with Python 3.6.8)

First clone this repo and cd into it    
then create a virtualenvironment with python virtualenv    
```bash
$ cd Cutout_Core
$ virtualenv -p python3 venv    
```

activate it      
```bash
$ . venv/bin/activate    
```

install all requirements from requirements.txt    
```bash
$ pip3 install -r requirements.txt    
```

if you get an astroquery version error you must install astroquery by           
(from within your virtualenv still) but cd out of th repo directory:    
```bash
$ cd ..     
$ git clone https://github.com/astropy/astroquery.git        
$ cd astroquery    
$ python setup.py install
```

Then remove the astroquery line from requirements.txt and run this again:  
```bash
$ pip3 install -r requirements.txt
```


You'll need to install <a target=_blank href="https://montage-wrapper.readthedocs.io/en/v0.9.5/#installation">Montage</a>, among other standard fair, which can be a littly tricky. Here's the Coles Notes:

 * Download `http://montage.ipac.caltech.edu/download/Montage_v5.0.tar.gz` from <a target=_blank href="http://montage.ipac.caltech.edu/docs/download2.html">`http://montage.ipac.caltech.edu/docs/download2.html`</a>.
 * `tar xvzf Montage_v5.0.tar.gz`
 * Got into the `Montage_v5.0` and type `make` to build.
 * Move the `Montage_v5.0` to `~/.montage/Montage_v5.0`, say, and add `~/.montage/Montage_v5.0/bin` to `$PATH`.
 * To test, run `mAdd` and you should see something like,<br>```[struct stat="ERROR", msg="Usage: mAdd [-d level] [-p imgdir] [-n(o-areas)] [-a mean|median|count] [-e(xact-size)] [-s statusfile] images.tbl template.hdr out.fits"]```<br>indicating it is installed correctly.

### How to run
from the command line:  
`$python3 fetch_cutouts.py   `

with Commands:    
  `fetch        Single cutout fetching command.   `     
  `fetch_batch  Batch cutout fetching command.   `     

Options:   
```text
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
```

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
```text
      Implemented surveys include:    
         - VLASS   
         - GLEAM    
            frequencies: f1 (072-103 MHz), f2 (103-034 MHz), f3 (139-170 MHz), f4 (170-231 MHz default)    
         - FIRST    
         - NVSS    
         - WISE    
            wavelengths: W1 (3.4μm default),  W2 (4.6μm),  W3 (12μm),  W4 (22μm)    
         - PANSTARRS    
            filters: g, r, i (default), z, y    
         - SDSS-I/II    
            filters: g (default), r, i    
            
        Filters/Frequencies/Wavelengths for each survey may be specified in the following formats:        
         > "WISE(w2),SDSS[g,r]"    
         > "WISE[w1],VLASS"    
         > "GLEAM(f1,f3)"    
         > WISE,VLASS    
 ```
 
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

### Output
This will fill `data_out` with the FITS files separated by Survey name directory.    
Success or failure results will be written to `OUTlog.txt`


