## list of some plotting routines used in UCCON1/5


import numpy as np, matplotlib.pyplot as plt, aplpy
from astropy import units as u
from astropy.coordinates import SkyCoord


###################################################################
###def params
testing = True ##is code being tested (i.e. don't apply to whole sample)
make_images = False
plt.interactive(testing)

fpath_in = '../data/rgz_qc/'
fpath_imout = '../data/rgz_qc/images/test_samp_n=200/'
fits_in = ('../VLASS/vlass_source_finding/cutouts/data/' +
           'Sample1_N=2465_27Jun2019/fits_4surveys/')

rms = 0.00013
rms_v = 0.00012
rms_f = 0.00013

rmarker = 'x'
rcolour = 'k'

imarker = 'o'
icolour = 'k'



######specific files
tpath = '../VLASS/vlass_source_finding/cutouts/cutout_testing/'
#tpath = '../VLASS/vlass_source_finding/cutouts/cutout_testing/'



vlass_file = tpath + 'data/Adrian/ClaRAN_testims_v1_26Jun19/J001108+123531_s3arcmin_VLASS.fits'
ps_file = tpath + 'data/Adrian/ClaRAN_testims_v1_26Jun19/J001108+123531_s3arcmin_PanSTARRS-i.fits'
wise_file = tpath + 'data/Adrian/ClaRAN_testims_v1_26Jun19/J001108+123531_s3arcmin_WISE-w1.fits'
first_file = tpath + 'data/out_3/J114522+152943_s3arcmin_FIRST.fits'

###################################################################
###################################################################
###def functions

###define contour levels to use
def rclevs(sigma, minlev=5, nlevs=10, levstep=2):
    ##define radio contour levels to use based on survey noise level
    level = minlev*sigma
    contour_levels = [level]
    for i in range(nlevs):
        if len(contour_levels) < nlevs:
            level = levstep*level
            contour_levels.append(level)
    
    return(contour_levels)

###creates an unframed png (if save==True) of optical/ir image with radio contour overlays
###contour levels defines above
def aplover_noframe(image, contours, fsize=4, cmap='Reds_r', ccol='darkblue', lw=0.8,
                    save=False, fname='testfig.png', vmx=10):
    ###create an aplpy contour overlay without a frame (suitable for Zoo)
    fig = plt.figure(figsize=(fsize, fsize), frameon=False)
    #    fig.set_size_inches(2.56, 2.56)
    #    fig.set_size_inches(5.12, 5.12)
    
    ###set up main image
    im = aplpy.FITSFigure(image, figure=fig, subplot=[0.0, 0.0, 1.0, 1.0])
    
    ##hide ticks, tick_labels, and axis_labels, set frame lw==0
    im.axis_labels.hide()
    im.tick_labels.hide()
    im.ticks.hide()
    im.frame.set_linewidth(0)
    
    ##plot colourmap of image
    im.show_colorscale(cmap=cmap, vmax=vmx, stretch='asinh')
    
    ##overlay contours
    im.show_contour(contours, levels=rclevs(sigma=rms), colors=ccol, linewidths=lw)
    
    if save == True:
        fig.savefig(fname, dpi=150)
        plt.close(fig=fig)
    return

###creates an unframed png (if save==True) of optical/ir image with radio contours and rgz components overlaid
###contour levels defines above
def aplover_noframe_rgz(image, contours, ir_ra, ir_dec, c_ra, c_dec, rpcol=rcolour,
                        rpmark=rmarker, ipcol=icolour, ipmark=imarker, fsize=5,
                        cmap='Reds_r', ccol='darkblue', lw=1, save=False,
                        fname='testfig.png', vmx=10):
    ###create an aplpy contour overlay without a frame (suitable for Zoo)
    fig = plt.figure(figsize=(fsize, fsize), frameon=False)
    
    ###set up main image
    im = aplpy.FITSFigure(image, figure=fig, subplot=[0.0, 0.0, 1.0, 1.0])
    
    ##hide ticks, tick_labels, and axis_labels, set frame lw==0
    im.axis_labels.hide()
    im.tick_labels.hide()
    im.ticks.hide()
    im.frame.set_linewidth(0)
    
    ##plot colourmap of image
    im.show_colorscale(cmap=cmap, vmax=vmx, stretch='asinh')
    
    ##overlay contours
    im.show_contour(contours, levels=rclevs(sigma=rms), colors=ccol, linewidths=lw, zorder=1)
    
    ###mark ir host and radio components
    im.show_markers(xw=c_ra, yw=c_dec, edgecolor=rpcol, facecolor=rpcol,
                    marker=rpmark, s=80, zorder=2)
    im.show_markers(xw=ir_ra, yw=ir_dec, edgecolor=ipcol, facecolor=ipcol,
                    marker=ipmark, s=100, zorder=3)
                    
    if save == True:
        fig.savefig(fname, dpi=150)
        plt.close(fig=fig)
    return

###comparison W+F and P+V images and versions with no overlaid contours for clarity
###contour levels defines above
def fw_pv_RGZcomp_incNO(wiseim, firstim, panim, vim, ir_ra, ir_dec, c_ra, c_dec,
                        rpcol=rcolour, rpmark=rmarker, ipcol=icolour, ipmark=imarker,
                        cmap_w='Reds_r', cmap_p='Blues_r', concol_f='darkblue', concol_v='r',
                        lw=0.8, rms_f=0.00015, rms_v=0.00013, fwidth=8, fheight=7,
                        maxw=10, maxp=10000, save=False, fname='test_comp.png'):
    
    fig = plt.figure(figsize=(fwidth, fheight))
    
    irim = aplpy.FITSFigure(wiseim, figure=fig, subplot=(2,2,1))
    #hide x-tick and axes labels
    irim.tick_labels.hide_x()
    irim.axis_labels.hide_x()
    
    opim = aplpy.FITSFigure(panim, figure=fig, subplot=(2,2,2))
    #hide x- and y-tick and axes labels
    opim.tick_labels.hide()
    opim.axis_labels.hide()
    
    irimo = aplpy.FITSFigure(wiseim, figure=fig, subplot=(2,2,3))
    opimo = aplpy.FITSFigure(panim, figure=fig, subplot=(2,2,4))
    #hide y-tick and axes labels
    opimo.tick_labels.hide_y()
    opimo.axis_labels.hide_y()
    
    ###show image
    irim.show_colorscale(cmap=cmap_w, stretch='asinh', vmax=maxw)
    opim.show_colorscale(cmap=cmap_p, stretch='asinh', vmax=maxp)
    irimo.show_colorscale(cmap=cmap_w, stretch='asinh', vmax=maxw)
    opimo.show_colorscale(cmap=cmap_p, stretch='asinh', vmax=maxp)
    
    ###add radio_contours
    irimo.show_contour(firstim, levels=rclevs(sigma=rms_f), colors=concol_f,
                       linewidths=lw, zorder=3)
    opimo.show_contour(vim, levels=rclevs(sigma=rms_v), colors=concol_v,
                       linewidths=lw, zorder=3)
                       
    ##add rgz info
    ###mark ir host and radio components
    irimo.show_markers(xw=c_ra, yw=c_dec, edgecolor=rpcol, facecolor=rpcol,
                       marker=rpmark, s=50, zorder=1)
    irimo.show_markers(xw=ir_ra, yw=ir_dec, edgecolor=ipcol, facecolor='none',
                       marker=ipmark, s=100, zorder=2)
                       
    opimo.show_markers(xw=c_ra, yw=c_dec, edgecolor=rpcol, facecolor=rpcol,
                       marker=rpmark, s=50, zorder=1)
    opimo.show_markers(xw=ir_ra, yw=ir_dec, edgecolor=ipcol, facecolor='none',
                       marker=ipmark, s=100, zorder=2)
                       
    if save == True:
        fig.savefig(fname, dpi=250)
        plt.close(fig=fig)

    return


###def vmin/max for ps -500, 5000

def aplover(image, contour, levs, colmap='Greys', concols='r'):
    fig = aplpy.FITSFigure(image, figsize=(6, 5))
    ##PS image
#    fig.show_colorscale(cmap=colmap, stretch='asinh', vmin=-500, vmax=12000)
    ##WISE image
    fig.show_colorscale(cmap=colmap, stretch='asinh', vmin=3.7, vmax=13.6)
    fig.show_contour(contour, levels=levs, colors=concols, linewidths=1.5)
    
    return


###extract coordinates from target name (of form Jhhmmss±ddmmss)
###may need adjusting to be more generic, works with target names we've used
def coords_from_name(names):
    ###returns astropy catalogue of coords from list of names of form Jhhmmss±ddmmss
    names = np.array(names)
    
    ra = [name[1:3]+'h'+name[3:5]+'m'+name[5:7]+'s' for name in names]
    dec = [name[7:10]+'d'+name[10:12]+'m'+name[12:]+'s' for name in names]
    
    cat = SkyCoord(ra=ra, dec=dec)
    
    return cat

###################################################################
###################################################################

