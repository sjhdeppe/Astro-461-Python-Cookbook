"""
Modified and simplified from Ian Cross's (UCSC) code available here:
https://people.ucsc.edu/~ianc/python/_modules/

2019-05-17 Stephanie Hamilton
"""

import numpy as np

try:
    from astropy.io import fits as pyfits
except:
    import pyfits

def calc_snr(data, axis=None):
    """
    Compute the quantity:
      data.mean(axis=axis)/data.std(axis=axis)
    for the specified data array/vector along the specified axis.

    Output will be a scalar (axis is None) or numpy array, as
    appropriate.
    """
    # 2010-09-02 08:10 IJC: Created
    # 2011-12-16 14:51 IJMC: Added optional nsigma flag
    # 2014-07-27 14:09 IJMC: Briefly documented nsigma flag.

    data = np.array(data)
    ret = data.mean(axis=axis)/data.std(axis=axis)

    return  ret

def makexflat(subreg, xord, nsigma=3, minsnr=10, minfrac=0.5, niter=1):
    """Helper function for XXXX.

    :INPUTS:
      subreg : 2D NumPy array
        spectral subregion, containing spectral background, sky,
        and/or target flux measurements.

      xord : scalar or sequence
        Order of polynomial by which each ROW will be normalized. If
        niter>0, xord can be a sequence of length (niter+1).  A good
        approach for, e.g., spectral dome flats is to set niter=1 and
        xord=[15,2].

      nsigma : scalar
        Sigma-clipping level for calculation of column-by-column S/N

      minsnr : scalar
        Minimum S/N value to use when selecting 'good' columns for
        normalization.

      minfrac : scalar, 0 < minfrac < 1
        Fraction of columns to use, selected by highest S/N, when
        selecting 'good' columns for normalization.

      niter : int
        Number of iterations.  If set to zero, do not iterate (i.e.,
        run precisely once through.)

    :NOTES:
      Helper function for functions XXXX
    """
    # 2013-01-20 14:20 IJMC: Created

    ny,  nx = subreg.shape
    xall = np.arange(nx)
    subreg_new = subreg.copy()

    iter = 0
    if not hasattr(xord, '__iter__'): xord = [xord]*(niter+1)
    while iter <= niter:
        snr = calc_snr(subreg_new, axis=0)
        ind = (snr > np.sort(snr)[-int(minfrac*snr.size)]) * (snr > minsnr)
        xxx = ind.nonzero()[0]
        norm_subreg = subreg[:,ind] / np.median(subreg[:,ind], 0)
        coefs = np.array([polyfitr(xxx, row, xord[iter], 3) for row in norm_subreg])
        xflat = np.array([np.polyval(coef0, xall) for coef0 in coefs])
        iter += 1
        subreg_new = subreg / xflat
    return xflat

def extractSubregion(fitsfilename, corners=None, dx=None, dy=None, kw=None, retall=False):
    """Extract a specified rectangular subregion from a FITS file.

    :INPUTS:
      fitsfilename : str
        Name of the (2D) FITS file.

      corners : str, 4-sequence
        if sequence: [x0, x1, y0, y1], corners of subregion.

        if str: header keyword containing this sequence.

        In either case, the extracted subregion (when dx=dy=0) will be:
          data[corners[2]:corners[3], corners[0]:corners[1]]

      dx : None, 2-sequence
        If sequence: [x0, x1] will become [x0-dx[0], x1+dx[1]]

      dy : None, 2-sequence
        If sequence: [y0, y1] will become [y0-dy[0], y1+dy[1]]

      kw : None, str
        If str: this header keyword will be updated with the new
        corners (possibly modified by dx, dy)

    :OUTPUTS:
      (subregion_data, [fitsfile_header, corners_used])

      If the specified header keyword is not found, or the specified
      corners return an error, then this function will crash inelegantly.

    :NOTE:
      WCS headers will not be updated, so be careful when using this
      routine for imaging data!
    """
    # 2012-08-28 15:25 IJMC: Created
    # 2013-01-20 14:50 IJMC: Clarified documentation of 'corners' input.

    try:
        from astropy.io import fits as pyfits
    except:
        import pyfits


    if dx is None:
        dx = 0
    if dy is None:
        dy = 0
    if not hasattr(dx, '__iter__'):
        dx = [dx]
    if not hasattr(dy, '__iter__'):
        dy = [dy]
    if len(dx) < 2:
        dx = [dx[0], dx[0]]
    if len(dy) < 2:
        dy = [dy[0], dy[0]]



    if not isinstance(fitsfilename, np.ndarray):
        data = pyfits.getdata(fitsfilename)
        header = pyfits.getheader(fitsfilename)
    else:
        data = fitsfilename
        header = pyfits.Header()

    ny, nx = data.shape


    newcorners = [max(0, corners[0]-dx[0]), min(nx, corners[1]+dx[1]), max(0, corners[2]-dy[0]), min(ny, corners[3]+dy[1])]
    subreg = data[newcorners[2]:newcorners[3], newcorners[0]:newcorners[1]]

    if kw is not None:
        header.update(kw, str(newcorners))

    if retall:
        ret = subreg, header, newcorners
    else:
        ret = subreg

    return ret


def polyfitr(x, y, N, s, fev=100, w=None, diag=False, clip='both',
             verbose=False, plotfit=False, plotall=False, eps=1e-13, catchLinAlgError=False):
    """Matplotlib's polyfit with weights and sigma-clipping rejection.

    :DESCRIPTION:
      Do a best fit polynomial of order N of y to x.  Points whose fit
      residuals exeed s standard deviations are rejected and the fit is
      recalculated.  Return value is a vector of polynomial
      coefficients [pk ... p1 p0].

    :OPTIONS:
        w:   a set of weights for the data; uses CARSMath's weighted polynomial
             fitting routine instead of numpy's standard polyfit.

        fev:  number of function evaluations to call before stopping

        'diag'nostic flag:  Return the tuple (p, chisq, n_iter)

        clip: 'both' -- remove outliers +/- 's' sigma from fit
              'above' -- remove outliers 's' sigma above fit
              'below' -- remove outliers 's' sigma below fit

        catchLinAlgError : bool
          If True, don't bomb on LinAlgError; instead, return [0, 0, ... 0].

    :REQUIREMENTS:
       :doc:`CARSMath`

    :NOTES:
       Iterates so long as n_newrejections>0 AND n_iter<fev.


     """
    # 2008-10-01 13:01 IJC: Created & completed
    # 2009-10-01 10:23 IJC: 1 year later! Moved "import" statements within func.
    # 2009-10-22 14:01 IJC: Added 'clip' options for continuum fitting
    # 2009-12-08 15:35 IJC: Automatically clip all non-finite points
    # 2010-10-29 09:09 IJC: Moved pylab imports inside this function
    # 2012-08-20 16:47 IJMC: Major change: now only reject one point per iteration!
    # 2012-08-27 10:44 IJMC: Verbose < 0 now resets to 0
    # 2013-05-21 23:15 IJMC: Added catchLinAlgError

    #from CARSMath import polyfitw
    from numpy import polyfit, polyval, isfinite, ones, array, std
    from numpy.linalg import LinAlgError
    from pylab import plot, legend, title, figure

    if verbose < 0:
        verbose = 0

    xx = array(x, copy=False)
    yy = array(y, copy=False)
    noweights = (w==None)
    if noweights:
        ww = ones(xx.shape, float)
    else:
        ww = array(w, copy=False)

    ii = 0
    nrej = 1

    if noweights:
        goodind = isfinite(xx)*isfinite(yy)
    else:
        goodind = isfinite(xx)*isfinite(yy)*isfinite(ww)

    xx2 = xx[goodind]
    yy2 = yy[goodind]
    ww2 = ww[goodind]

    while (ii<fev and (nrej>0)):
        p = polyfit(xx2,yy2,N)
        residual = yy2 - polyval(p,xx2)
        stdResidual = std(residual)
        clipmetric = s * stdResidual

        if clip=='both':
            worstOffender = abs(residual).max()
            #pdb.set_trace()
            if worstOffender <= clipmetric or worstOffender < eps:
                ind = ones(residual.shape, dtype=bool)
            else:
                ind = abs(residual) < worstOffender
        elif clip=='above':
            worstOffender = residual.max()
            if worstOffender <= clipmetric:
                ind = ones(residual.shape, dtype=bool)
            else:
                ind = residual < worstOffender
        elif clip=='below':
            worstOffender = residual.min()
            if worstOffender >= -clipmetric:
                ind = ones(residual.shape, dtype=bool)
            else:
                ind = residual > worstOffender
        else:
            ind = ones(residual.shape, dtype=bool)

        xx2 = xx2[ind]
        yy2 = yy2[ind]
        if (not noweights):
            ww2 = ww2[ind]
        ii = ii + 1
        nrej = len(residual) - len(xx2)
        if plotall:
            figure()
            plot(x,y, '.', xx2,yy2, 'x', x, polyval(p, x), '--')
            legend(['data', 'fit data', 'fit'])
            title('Iter. #' + str(ii) + ' -- Close all windows to continue....')

        if verbose:
            print( str(len(x)-len(xx2)) + ' points rejected on iteration #' + str(ii))

    if (plotfit or plotall):
        figure()
        plot(x,y, '.', xx2,yy2, 'x', x, polyval(p, x), '--')
        legend(['data', 'fit data', 'fit'])
        title('Close window to continue....')

    if diag:
        chisq = ( (residual)**2 / yy2 ).sum()
        p = (p, chisq, ii)

    return p


def make_spectral_flats(domeflat, subreg_corners, badpixelmask=None, xord_pix=[15,2], 
                        xord_sky=[2,1], yord=2, minsnr=5, minfrac_pix=0.7, minfrac_sky=0.5, 
                        locs=None):
    """
    Construct appropriate corrective frames for multi-object
    spectrograph data.  Specifically: corrections for irregular slit
    widths, and pixel-by-pixel detector sensitivity variations.

    sky : 2D NumPy array
       Master spectral sky frame, e.g. from median-stacking many sky
       frames or masking-and-stacking dithered science spectra frames.
       This frame is used to construct a map to correct science frames
       (taken with the identical slit mask!) for irregular sky
       backgrounds resulting from non-uniform slit widths (e.g.,
       Keck/MOSFIRE).

       Note that the dispersion direction should be 'horizontal'
       (i.e., parallel to rows) in this frames.

    domeflat : 2D NumPy array
       Master dome spectral flat, e.g. from median-stacking many dome
       spectra.  This need not be normalized in the dispersion or
       spatial directions. This frame is used to construct a flat map
       of the pixel-to-pixel variations in detector sensitivity.

       Note that the dispersion direction should be 'horizontal'
       (i.e., parallel to rows) in this frames.

    subreg_corners : sequence of 2- or 4-sequences
        Indices (or merely starting- and ending-rows) for each
        subregion of interest in 'sky' and 'domeflat' frames.  Syntax
        should be that of :func:`tools.extractSubregion`, or such that
        subreg=sky[subreg_corners[0]:subreg_corners[1]]

    badpixelmask : 2D NumPy array, or str
        Nonzero for any bad pixels; these will be interpolated over
        using :func:`nsdata.bfixpix`.

    xord_pix : sequence
        Polynomial orders for normalization in dispersion direction of
        pixel-based flat (dome flat) on successive iterations; see
        :func:`makexflat`.

    xord_sky : sequence
        Polynomial orders for normalization in dispersion direction of
        slit-based correction (sky frame) on successive iterations;
        see :func:`makexflat`.

    yord : scalar
        Polynomial order for normalization of pixel-based (dome) flat
        in spatial direction.

    locs : None, or sequence
        Row-index in each subregion of the location of the
        spectral-trace-of-interest.  Only used for rectifying of sky
        frame; leaving this at None will have at most a mild
        effect. If not None, should be the same length as
        subreg_corners.  If subreg_corners[0] = [800, 950] then
        locs[0] might be set to, e.g., 75 if the trace lies in the
        middle of the subregion.


    :RETURNS:
        wideslit_skyflat, narrowslit_domeflat

    :EXAMPLE:
     ::


        try:
            from astropy.io import fits as pyfits
        except:
            import pyfits

        import spec
        import ir

        obs = ir.initobs('20121010')
        sky = pyfits.getdata(obs._raw + 'masktersky.fits')
        domeflat = pyfits.getdata(obs._raw + 'mudflat.fits')
        allinds = [[7, 294], [310, 518], [532, 694], [710, 960], [976, 1360], [1378, 1673], [1689, 2022]]
        locs = [221, 77, 53, 88, 201, 96, 194]

        skycorrect, pixcorrect = spec.make_spectral_flats(sky, domeflat, allinds, obs.badpixelmask, locs=locs)
    """
    # 2013-01-20 14:40 IJMC: Created

    # Check inputs:
    if not isinstance(domeflat, np.ndarray):
        domeflat = pyfits.getdata(domeflat)

    ny0, nx0 = domeflat.shape

    narrowslit_domeflat = np.ones((ny0, nx0), dtype=float)

    # Loop through all subregions specified:
    for jj in range(len(subreg_corners)):
        # Define extraction & alignment indices:
        specinds = np.array(subreg_corners[jj]).ravel().copy()
        if len(specinds)==2:
            specinds = np.concatenate(([0, nx0], specinds))
        if locs is None:
            loc = np.mean(specinds[2:4])
        else:
            loc = locs[jj]

        domeflatsub = extractSubregion(domeflat, specinds, retall=False)
        ny, nx = domeflatsub.shape
        yall = np.arange(ny)

        # Normalize Dome Spectral Flat in X-direction (rows):
        xflat_dome = makexflat(domeflatsub, xord_pix, minsnr=minsnr,
                               minfrac=minfrac_pix, niter=len(xord_pix)-1)
        ymap = domeflatsub*0.0
        xsubflat = domeflatsub / np.median(domeflatsub, 0) / xflat_dome

        # Normalize Dome Spectral Flat in Y-direction (columns):
        ycoefs = [polyfitr(yall, xsubflat[:,ii], yord, 3) for ii in range(nx)]
        yflat = np.array([np.polyval(ycoef0, yall) for ycoef0 in ycoefs]).transpose()
        pixflat = domeflatsub / (xflat_dome * yflat * np.median(domeflatsub, 0))
        bogus = (pixflat< 0.02) + (pixflat > 50)
        pixflat[bogus] = 1.

        narrowslit_domeflat[specinds[2]:specinds[3], specinds[0]:specinds[1]] = pixflat

        print("Done with subregion %i" % (jj+1))

    return narrowslit_domeflat
