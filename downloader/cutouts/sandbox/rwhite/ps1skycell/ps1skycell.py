"""ps1skycell.py    Locate PanSTARRS skycell for a given RA/Dec

vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4 ai :

R. White, 2019 July 12
"""

from __future__ import print_function
import sys, os
import numpy as np
from astropy.io import fits

# pixel scale is 0.25 arcsec
pixscale = 0.25

# read table of rings info from table in first extension of FITS file
gridfits = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'ps1grid.fits')
rings = fits.open(gridfits)[1].data
dec_limit = rings.field('dec_min').min()

def findskycell(ra, dec, all=False):

    """Given input arrays RA, DEC (degrees), returns a dictionary with the skycell info

    RA and Dec can be scalars or 1-D arrays.

    all is a boolean indicating that all overlapping subcells should be returned.  By default
    this returns just the best sky cell for each ra/dec (the one where the position is closest to the center).
    If all=True, this returns neighboring skycells **from the same projection cell** that also
    cover the RA/Dec position.  Note that this does not attempt to identify overlapping projection
    cells in regions near the edges of the large 4x4 degree cells.
    If all=True, multiple matches for a ra,dec position are grouped together with the best image
    listed first.

    This uses the ps1grid.fits file for info on the tessellation.

    The return dictionary is indexed by column name, with each column being an array of values:
    Column      Value
    index       Index into original input list
    ra          Input position
    dec         
    projcell    Projection cell (0 if outside coverage)
    subcell     Subcell (0..99)
    crval1      Sky position of reference pixel for projection cell
    crval2      
    crpix1      Reference pixel in this skycell
    crpix2      
    x           Source pixel position in this skycell
    y           
    """

    # force ra and dec to be arrays
    if np.isscalar(ra):
        ra = [ra]
    if np.isscalar(dec):
        dec = [dec]
    ra = np.asarray(ra)
    dec = np.asarray(dec)
    if ra.shape != dec.shape or ra.ndim != 1:
        raise ValueError("ra and dec must both be scalars or be matching 1-D arrays")

    # find dec zone where rings.dec_min <= dec < rings.dec_max
    idec = np.searchsorted(rings.field('dec_max'), dec)

    # special handling at pole where overlap is complicated
    # do extra checks for top 2 rings
    # always start with the ring just below the pole
    nearpole = np.where(idec >= len(rings)-2)[0]
    idec[nearpole] = len(rings)-2

    getfield = rings[idec].field
    nband = getfield('nband')
    # get normalized RA in range 0..360
    nra = ra % 360.0
    ira = np.rint(nra*nband/360.0).astype(int) % nband

    projcell = getfield('projcell') + ira
    dec_cen = getfield('dec')
    ra_cen = ira*360.0/nband

    # locate subcell(s) within the projection cell

    # use tangent project to get pixel offsets
    x, y = sky2xy_tan(nra, dec, ra_cen, dec_cen)

    pad = 480

    if len(nearpole) > 0:
        # handle the points near the pole (if any)
        # we know that this "ring" has only a single field
        getfield2 = rings[-1].field
        projcell2 = getfield2('projcell')
        dec_cen2 = getfield2('dec')
        ra_cen2 = 0.0
        x2, y2 = sky2xy_tan(nra[nearpole], dec[nearpole], ra_cen2, dec_cen2)
        # compare x,y and x2,y2 to image sizes to select best image
        # returns a Boolean array with true for values where 2nd image is better
        use2 = poleselect(x[nearpole], y[nearpole], x2, y2, rings[-2], rings[-1], pad)
        if use2.any():
            wuse2 = np.where(use2)[0]
            w2 = nearpole[wuse2]
            idec[w2] = len(rings)-1
            getfield = rings[idec].field
            nband[w2] = 1
            ira[w2] = 0
            projcell[w2] = projcell2
            dec_cen[w2] = dec_cen2
            ra_cen[w2] = ra_cen2
            x[w2] = x2[wuse2]
            y[w2] = y2[wuse2]

    # compute the subcell from the pixel location
    px = getfield('xcell')-pad
    py = getfield('ycell')-pad
    k = np.rint(4.5+x/px).astype(int).clip(0,9)
    j = np.rint(4.5+y/py).astype(int).clip(0,9)
    subcell = 10*j + k

    # get pixel coordinates within the skycell image
    crpix1 = getfield('crpix1') + px*(5-k)
    crpix2 = getfield('crpix2') + py*(5-j)
    ximage = x + crpix1
    yimage = y + crpix2

    # insert zeros where we are below lowest dec_min
    w = np.where(dec < dec_limit)
    projcell[w] = 0
    subcell[w] = 0
    crpix1[w] = 0
    crpix2[w] = 0
    ximage[w] = 0
    yimage[w] = 0
    index = np.arange(len(dec))

    if all:
        # check neighboring skycells in same projection cell to see if there is overlap
        dec0 = dec.copy()
        for xoff in [-1,0,1]:
            for yoff in [-1,0,1]:
                if (xoff == 0 and yoff == 0):
                    continue
                k2 = k+xoff
                j2 = j+yoff
                w = np.where((k2 >= 0) & (k2 <= 9) & (j2 >= 0) & (j2 <= 9) & (dec0 >= dec_limit))[0]
                if len(w) == 0:
                    continue
                # indexing below works because i,j are not modified and so the w array indexes into the 
                # original values for projcell, subcell, ra_cen, dec_cen, etc.
                k2 = k2[w]
                j2 = j2[w]
                px2 = px[w]
                py2 = py[w]
                projcell2 = projcell[w]
                subcell2 = 10*j2 + k2
                ra2 = ra[w]
                dec2 = dec[w]
                ra_cen2 = ra_cen[w]
                dec_cen2 = dec_cen[w]
                crpix12 = crpix1[w] - px2*xoff
                crpix22 = crpix2[w] - py2*yoff
                ximage2 = x[w] + crpix12
                yimage2 = y[w] + crpix22
                wk = np.where((ximage2 >= 0) & (ximage2 < px2+pad) & (yimage2 >= 0) & (yimage2 < py2+pad))[0]
                if len(wk) > 0:
                    ra = np.append(ra, ra2[wk])
                    dec = np.append(dec, dec2[wk])
                    projcell = np.append(projcell, projcell2[wk])
                    subcell = np.append(subcell, subcell2[wk])
                    ra_cen = np.append(ra_cen, ra_cen2[wk])
                    dec_cen = np.append(dec_cen, dec_cen2[wk])
                    crpix1 = np.append(crpix1, crpix12[wk])
                    crpix2 = np.append(crpix2, crpix22[wk])
                    ximage = np.append(ximage, ximage2[wk])
                    yimage = np.append(yimage, yimage2[wk])
                    index = np.append(index, w[wk])
        # re-sort to get entries for same object grouped together
        # use stable sort so best image comes first
        isort = np.argsort(index, kind='mergesort')
        ra = ra[isort]
        dec = dec[isort]
        projcell = projcell[isort]
        subcell = subcell[isort]
        ra_cen = ra_cen[isort]
        dec_cen = dec_cen[isort]
        crpix1 = crpix1[isort]
        crpix2 = crpix2[isort]
        ximage = ximage[isort]
        yimage = yimage[isort]
        index = index[isort]

    # return a dictionary
    return {'index': index, 'ra': ra, 'dec': dec, 'projcell': projcell, 'subcell': subcell,
        'crval1': ra_cen, 'crval2': dec_cen, 'crpix1': crpix1, 'crpix2': crpix2,
        'x': ximage, 'y': yimage}


def poleselect(x1, y1, x2, y2, rings1, rings2, pad):

    """Compares x,y values from 2 images to determine which is best

    Returns boolean array with True where x2,y2 is best
    """

    nx1 = 10*(rings1['xcell']-pad)+pad
    ny1 = 10*(rings1['ycell']-pad)+pad
    nx2 = 10*(rings2['xcell']-pad)+pad
    ny2 = 10*(rings2['ycell']-pad)+pad
    # compute minimum distances to image edges
    # note negative values are off edge
    d1 = np.minimum(np.minimum(x1+nx1//2,nx1//2-1-x1), np.minimum(y1+ny1//2,ny1//2-1-y1))
    d2 = np.minimum(np.minimum(x2+nx2//2,nx2//2-1-x2), np.minimum(y2+ny2//2,ny2//2-1-y2))
    return (d1 < d2)


def getskycell_center(projcell, subcell):

    """Return position of the center (RA, DEC) of a skycell or list of skycells

    projcell and subcell can be scalars or arrays.
    This uses the ps1grid.fits file for info on the tessellation.
    This just clips to the nearest cell in the grid for values outside the bounds

    The return dictionary is indexed by column name, with each column being an array of values:
    Column      Value
    projcell    Input parameters
    subcell     
    ra          Central position (degrees)
    dec         
    crval1      Sky position of reference pixel for projection cell
    crval2      
    crpix1      Reference pixel in this skycell
    crpix2      
    x           Center pixel position in this skycell
    y           
    """

    # force projcell and subcell to be arrays
    if np.isscalar(projcell):
        projcell = [projcell]
    if np.isscalar(subcell):
        subcell = [subcell]
    projcell = np.asarray(projcell)
    subcell = np.asarray(subcell)
    if projcell.shape != subcell.shape:
        raise ValueError("projcell and subcell must both be scalars or be matching shape arrays")

    # find dec zone where rings.projcell <= projcell
    idec = (np.searchsorted(rings.field('projcell'), projcell+1)-1).clip(0)
    getfield = rings[idec].field
    nband = getfield('nband')
    ira = (projcell - getfield('projcell')).clip(0,nband)
    dec_cen = getfield('dec')
    ra_cen = ira*360.0/nband

    # locate subcell within the projection cell
    k = subcell % 10
    j = (subcell//10) % 10

    pad = 480
    px = getfield('xcell')-pad
    py = getfield('ycell')-pad
    ximage = 0.5*(px+pad-1)
    yimage = 0.5*(py+pad-1)

    # get pixel coordinates within the skycell image
    crpix1 = getfield('crpix1') + px*(5-k)
    crpix2 = getfield('crpix2') + py*(5-j)
    # position in projection cell
    x = ximage - crpix1
    y = yimage - crpix2

    ra, dec = xy2sky_tan(x, y, ra_cen, dec_cen)

    # return a dictionary
    return {'ra': ra, 'dec': dec, 'projcell': projcell, 'subcell': subcell,
        'crval1': ra_cen, 'crval2': dec_cen, 'crpix1': crpix1, 'crpix2': crpix2,
        'x': ximage, 'y': yimage}


def sky2xy_tan(ra, dec, ra_cen, dec_cen, crpix=(0.0,0.0)):

    """Convert RA,Dec sky position (degrees) to X,Y pixel position

    ra[n], dec[n] are input arrays in degrees
    ra_cen[n], dec_cen[n] are image centers in degrees
    crpix is the reference pixel position (x,y)
    Returns tuple (x,y) where x and y are arrays with pixel position for each RA,Dec
    """

    dtor = np.pi/180
    cd00 = -pixscale*dtor/3600
    cd01 = 0.0
    cd10 = 0.0
    cd11 = -cd00
    determ = cd00*cd11-cd01*cd10
    cdinv00 =  cd11/determ
    cdinv01 = -cd01/determ
    cdinv10 = -cd10/determ
    cdinv11 =  cd00/determ

    cos_crval1 = np.cos(dtor*dec_cen)
    sin_crval1 = np.sin(dtor*dec_cen)

    radif = (ra - ra_cen)*dtor
    w = np.where(radif > np.pi)
    radif[w] -= 2*np.pi
    w = np.where(radif < -np.pi)
    radif[w] += 2*np.pi

    decrad = dec*dtor
    cos_dec = np.cos(decrad)
    sin_dec = np.sin(decrad)
    cos_radif = np.cos(radif)
    sin_radif = np.sin(radif)
    h = sin_dec*sin_crval1 + cos_dec*cos_crval1*cos_radif
    xsi = cos_dec*sin_radif/h
    eta = (sin_dec*cos_crval1 - cos_dec*sin_crval1*cos_radif)/h
    xdif = cdinv00*xsi + cdinv01*eta
    ydif = cdinv10*xsi + cdinv11*eta
    return (xdif+crpix[0], ydif+crpix[1])


def xy2sky_tan(x, y, ra_cen, dec_cen, crpix=(0.0,0.0)):

    """Convert X,Y pixel position to RA,Dec sky position (degrees)

    x[n], y[n] are input arrays in pixels
    ra_cen[n], dec_cen[n] are image centers in degrees
    crpix is the reference pixel position (x,y)
    Returns tuple (ra,dec) where ra and dec are arrays with pixel position for each x,y
    """

    dtor = np.pi/180
    cd00 = -pixscale*dtor/3600
    cd01 = 0.0
    cd10 = 0.0
    cd11 = -cd00

    cos_crval1 = np.cos(dtor*dec_cen)
    sin_crval1 = np.sin(dtor*dec_cen)

    xdif = x - crpix[0]
    ydif = y - crpix[1]
    xsi = cd00*xdif + cd01*ydif
    eta = cd10*xdif + cd11*ydif
    beta = cos_crval1 - eta*sin_crval1
    ra = np.arctan2(xsi, beta) + dtor*ra_cen
    gamma = np.sqrt(xsi**2 + beta**2)
    dec = np.arctan2(eta*cos_crval1+sin_crval1, gamma)
    return (ra/dtor, dec/dtor)


if __name__ == "__main__":
    # test findskycell at a bunch of points, returning all matching skycells for each point
    tdec = np.arange(31)*3.95 - 29.1
    tra = np.arange(31)*12.
    table = findskycell(tra,tdec,all=True)
    projcell = table['projcell']
    subcell = table['subcell']
    ra = table['ra']
    dec = table['dec']
    index = table['index']
    for i in range(len(dec)):
        print("%2d %11.6f %10.6f skycell.%4.4d.%3.3d" % (index[i], ra[i], dec[i], projcell[i], subcell[i]))
