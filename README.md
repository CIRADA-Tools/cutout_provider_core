## Continuum Common

This repo contains common components for development of the VLASS Continuum Database Components.

### Cutout Software Installation Notes

You'll need to install <a target=_blank href="https://montage-wrapper.readthedocs.io/en/v0.9.5/#installation">Montage</a>, among other standard fair, which can be a littly tricky. Here's the Coles Notes:

 * Download `http://montage.ipac.caltech.edu/download/Montage_v5.0.tar.gz` from <a target=_blank href="http://montage.ipac.caltech.edu/docs/download2.html">`http://montage.ipac.caltech.edu/docs/download2.html`</a>.
 * `tar xvzf Montage_v5.0.tar.gz`
 * Got into the `Montage_v5.0` and type `make` to build.
 * Move the `Montage_v5.0` to `~/.montage/Montage_v5.0`, say, and add `~/.montage/Montage_v5.0/bin` to `$PATH`.
 * To test, run `mAdd` and you should see something like,<br>```[struct stat="ERROR", msg="Usage: mAdd [-d level] [-p imgdir] [-n(o-areas)] [-a mean|median|count] [-e(xact-size)] [-s statusfile] images.tbl template.hdr out.fits"]```<br>indicating it is installed correctly.

### VOSpace Installation Notes:

RE: http://www.canfar.net/en/docs/storage/

Once anaconda is installed do,

    pip install --user -U vos

and then add ~/.local/bin to $PATH.

Then to initialize, run

    getCert

and enter your username and password.

### Montage Installation Notes:

This installation is trick, so here's the Coles Notes:

 * Download `http://montage.ipac.caltech.edu/download/Montage_v5.0.tar.gz` from <a target=_blank href="http://montage.ipac.caltech.edu/docs/download2.html">`http://montage.ipac.caltech.edu/docs/download2.html`</a>.
 * `tar xvzf Montage_v5.0.tar.gz`
 * Got into the `Montage_v5.0` and type `make` to build.
 * Move the `Montage_v5.0` to `~/.montage/Montage_v5.0`, say, and add `~/.montage/Montage_v5.0/bin` to `$PATH`.
 * To test, run `mAdd` and you should see something like,<br>```[struct stat="ERROR", msg="Usage: mAdd [-d level] [-p imgdir] [-n(o-areas)] [-a mean|median|count] [-e(xact-size)] [-s statusfile] images.tbl template.hdr out.fits"]```<br>indicating it is installed correctly.

### Pyhton Packages:

Also lot of the sandard astronomy pakacges are used, but here's the hightlights:

* astropy
* atroquery
* urllib
* pickle
* yaml

