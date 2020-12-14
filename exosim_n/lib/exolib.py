import numpy as np
from scipy import signal
from scipy import interpolate
from scipy.integrate import cumtrapz
import scipy.special
from astropy import units as u
from astropy import constants as const
from astropy.io import fits
import sys, os
import matplotlib.pyplot as plt
import pandas as pd


def exosim_error(error_msg):
    sys.stderr.write("Error code: {:s}\n".format(error_msg))
    sys.exit(0)
    
def exosim_msg(msg, prefix = None):
  msg = msg if prefix==None else "[%s]: %s\n"%(prefix, msg)
  print(msg)
  
def exosim_plot(figure_name, show, image=False, xlabel=None, 
                 ylabel=None, image_data=None, xdata=None, ydata=None, marker='-',
                 xlim=None, ylim=None, aspect=None, interpolation=None, linewidth=1, alpha=1,
                 grid=False, label=None):
  if hasattr(image_data, "unit"): # fid for astropy image issue with units
      image_data = image_data.value
  if show == 1:
      if image == True:
          plt.figure(figure_name)
          if aspect or interpolation:
              plt.imshow(image_data, aspect=aspect, interpolation=None)
          else:    
              plt.imshow(image_data)
          if xlabel:
              plt.xlabel(xlabel)
          if ylabel:
              plt.ylabel(ylabel)         
              
      else:
          plt.figure(figure_name)
 
          if xdata is not None:
                plt.plot(xdata, ydata, marker, linewidth=linewidth, alpha=alpha, label=label)
          else:
                plt.plot(ydata, marker,linewidth=linewidth, alpha=alpha, label=label)
                        
          if xlabel:
              plt.xlabel(xlabel)
          if ylabel:
              plt.ylabel(ylabel)
          if ylim:
              plt.ylim(ylim[0],ylim[1])
          if xlim:
              plt.xlim(xlim[0],xlim[1])
          if grid==True:
              plt.grid()
          if label is not None:
              plt.legend(loc='best')

  
  
def calc_logg(m,r):
    m = m.to(u.kg)
    r = r.to(u.m)
    # G = spc.G*u.m**3/u.kg**1/u.s**2
    G = const.G
    g = G*m /r**2
    g0 = g.to(u.cm/u.s**2)
    return g, np.log10(g0.value)
    
def calc_EqT(T_s,R_s,a,A,tidal):   
    if tidal == 1: 
        factor = 2.
    else:
        factor = 4.
    T_p = T_s*(R_s.to(u.m)/a.to(u.m))**0.5*((1-A)/factor)**0.25
     
    return T_p   
  
def add_nonlin(data, FW):
  count = data
  x= count
  y1 = x
  a = (0.05*FW-1) / FW**2
  b=1;c=1
  y2 = a*x**2 + b*x + c
  y = y1-(y2-y1)
  mod_count = y
  
  #make a non-linearity coefficient matrix for data reduction
  coeff= [a,b,c]
  lincorr_matrix = np.zeros((data.shape[0], data.shape[1], 3))
  for i in range(lincorr_matrix.shape[2]):
      lincorr_matrix[...,i] = lincorr_matrix[...,i] + coeff[i]
  
  return mod_count, lincorr_matrix
  
  
def add_nonlin_2(data, FW):
    count = data
    aa = np.array(pd.read_csv('/Users/user1/Downloads/nonlin.csv', header=None))
    x = aa[:,0]*1e4
    y = aa[:,1]*1e4
    r = 3 # degree 7 works well      
    z0 = np.polyfit(x, y, r)
    y0 =0
    for i in range (0,r+1):
        y0 = y0 + z0[i]*count**(r-i) 
        
    count = y0
    
    return count 
  
   
  
def logbin(x, a,  R, xmin=None, xmax=None):
  n = a.size
  imin = 0
  imax = n-1
  
  if xmin == None or xmin < x.min(): xmin = x.min()
  if xmax == None or xmax > x.max(): xmax = x.max()
  
  idx = np.argsort(x)
  xp = x[idx]
  yp = a[idx]
  
  delta_x = xmax/R
  N = 20.0 * (xmax-xmin)/delta_x
  _x = np.linspace(xmin,xmax, N)
  _y = np.interp(_x, xp, yp)
  
  nbins = 1+np.round( np.log(xmax/xmin)/np.log(1.0 + 1.0/R) ).astype(np.int)
  bins = xmin*np.power( (1.0+1.0/R), np.arange(nbins))
  
  slices  = np.searchsorted(_x, bins)
  counts = np.ediff1d(slices)
  
  mean = np.add.reduceat(_y, slices[:-1])/(counts)
  bins = 0.5*(bins[:-1] + bins[1:])
  return bins[:-1], mean[:-1]
 
def rebin(x, xp, fp):
  ''' Resample a function fp(xp) over the new grid x, rebinning if necessary, 
    otherwise interpolates
    Parameters
    ----------
    x	: 	array like
	New coordinates
    fp 	:	array like
	y-coordinates to be resampled
    xp 	:	array like
	x-coordinates at which fp are sampled
	
    Returns
    -------
    out	: 	array like
	new samples  
  '''
  
  if (x.unit != xp.unit):
    print (x.unit, xp.unit)
    exosim_error('Units mismatch')
  
  idx = np.where(np.logical_and(xp > 0.9*x.min(), xp < 1.1*x.max()))[0]

  xp = xp[idx]
  fp = fp[idx]
  

  if np.diff(xp).min() < np.diff(x).min():
   
    # Binning!
    c = cumtrapz(fp, x=xp)
    xpc = xp[1:]
        
    delta = np.gradient(x)
    new_c_1 = np.interp(x-0.5*delta, xpc, c, 
                        left=0.0, right=0.0)
    new_c_2 = np.interp(x+0.5*delta, xpc, c, 
                        left=0.0, right=0.0)
    new_f = (new_c_2 - new_c_1)/delta
    

 
  else:
    # Interpolate !
    new_f = np.interp(x, xp, fp, left=0.0, right=0.0)
    
  new_f = (new_f.value)*fp.unit
  
#    func = interpolate.interp1d(xp, fp, kind='quadratic', bounds_error=None, fill_value=0.0)
#    new_f  = func(x)*fp.unit
  '''
  import matplotlib.pyplot as plt
  plt.plot(xp, fp, '-')
  plt.plot(x, new_f, '.-')
  plt.show()
  # check
  print np.trapz(new_f, x)
  idx = np.where(np.logical_and(xp>= x.min(), xp <= x.max()))
  print np.trapz(fp[idx], xp[idx])
  '''
  return x, new_f
  
def fast_convolution(im, delta_im, ker, delta_ker):
  """ fast_convolution.
    Convolve an image with a kernel. Image and kernel can be sampled on different
      grids defined.
    
    Parameters
    __________
      im : 			array like
				the image to be convolved
      delta_im :		scalar
				image sampling interval
      ker : 			array like
				the convolution kernel
      delta_ker :		scalar
				Kernel sampling interval
    Returns
    -------
      spectrum:			array like
				the image convolved with the kernel.
  """
  fc_debug = False
  # Fourier transform the kernel
  kerf = (np.fft.rfft2(ker))
  ker_k = [ np.fft.fftfreq(ker.shape[0], d=delta_ker),
	   np.fft.rfftfreq(ker.shape[1], d=delta_ker) ]
  ker_k[0] = np.fft.fftshift(ker_k[0])
  kerf     = np.fft.fftshift(kerf, axes=0)
  
  # Fourier transform the image
  imf  = np.fft.rfft2(im)
  im_k = [ np.fft.fftfreq(im.shape[0], d=delta_im),
	   np.fft.rfftfreq(im.shape[1], d=delta_im) ]
  im_k[0] = np.fft.fftshift(im_k[0])
  imf     = np.fft.fftshift(imf, axes=0)
  
  # Interpolate kernel 
  kerf_r = interpolate.RectBivariateSpline(ker_k[0], ker_k[1],
					   kerf.real)
  kerf_i = interpolate.RectBivariateSpline(ker_k[0], ker_k[1],
					   kerf.imag)
  if (fc_debug):
    pl.plot(ker_k[0], kerf[:, 0].real,'.r')
    pl.plot(ker_k[0], kerf[:, 0].imag,'.g')
    pl.plot(im_k[0], kerf_r(im_k[0], im_k[1])[:, 0],'-r')
    pl.plot(im_k[0], np.abs(imf[:, 0]),'-b')

  # Convolve
  imf = imf * (kerf_r(im_k[0], im_k[1]) + 1j*kerf_i(im_k[0], im_k[1])) 
  
  if (fc_debug):
    pl.plot(im_k[0], np.abs(imf[:, 0]),'-y')

  imf = np.fft.ifftshift(imf, axes=0)
  
  return np.fft.irfft2(imf)*(delta_ker/delta_im)**2

   
def planck(wl, T):
  """ Planck function. 
    
    Parameters
    __________
      wl : 			array
				wavelength [micron]
      T : 			scalar
				Temperature [K]
				Spot temperature [K]
    Returns
    -------
      spectrum:			array
				The Planck spectrum  [W m^-2 sr^-2 micron^-1]
  """
    
  a = np.float64(1.191042768e8)*u.um**5 *u.W/ u.m**2 /u.sr/u.um
  b = np.float64(14387.7516)*1*u.um * 1*u.K
  try:
    x = b/(wl*T)
    bb = a/wl**5 / (np.exp(x) - 1.0)
  except ArithmeticError:
    bb = np.zeros(np.size(wl))
  return bb
 
 
def sed_propagation(sed, transmission, emissivity=None, temperature = None):
  sed.sed = sed.sed*transmission.sed
  if emissivity and temperature:
    sed.sed = sed.sed + emissivity.sed*planck(sed.wl, temperature)

  return sed
  
def Psf_Interp(zfile, delta_pix, WavRange):
    ''' 
    PSF Interpolation
    Parametes
    ---------
        zfile : string
            input PSF fits file
        Delta : scalar
            Sampling interval in micron
        WavRange : ndarray
            array of wavelengths in micron
    
    Returns
    -------
        PSF interpolated data cube. Area normalised to unity.
        
    '''
    delta_pix =  delta_pix.value

    hdulist = pyfits.open(zfile)    
    NAXIS1, NAXIS2 = hdulist[0].header['NAXIS1'], hdulist[0].header['NAXIS2']
    in_ph_size_x, in_ph_size_y = hdulist[0].header['CDELT1']*NAXIS1, hdulist[0].header['CDELT2']*NAXIS2
    num_pix_x, num_pix_y = np.trunc(in_ph_size_x/delta_pix).astype(np.int), np.trunc(in_ph_size_y/delta_pix).astype(np.int)
               
    inwl   = np.zeros(len(hdulist))

    redata = np.zeros((num_pix_y, num_pix_x, len(hdulist)))

    xin = np.linspace(-1.0, 1.0, NAXIS1)
    yin = np.linspace(-1.0, 1.0, NAXIS2)

    xout = np.linspace(-1.0, 1.0, num_pix_x)
    yout = np.linspace(-1.0, 1.0, num_pix_y)

    for i, hdu in enumerate(hdulist):
        inwl[i]   = np.float64(hdu.header['WAV'])        
        f = interpolate.RectBivariateSpline(xin, yin, hdu.data)
        redata[..., i] = f(xout,yout)
        redata[..., i] /= redata[..., i].sum()
    return interpolate.interp1d(inwl, redata, axis=2, bounds_error=False, fill_value=0.0, kind='quadratic')(WavRange)


def Psf_photometer(zfile, delta_pix, osf, fpn, WavRange):
    ''' 
    PSF Interpolation
    Parametes
    ---------
        zfile : string
            input PSF fits file
        Delta : scalar
            Sampling interval in micron
        WavRange : ndarray
            array of wavelengths in micron
    
    Returns
    -------
        PSF interpolated data cube. Area normalised to unity.
        
    '''
    
    delta_pix =  delta_pix.value

    hdulist = pyfits.open(zfile)    
    NAXIS1, NAXIS2 = hdulist[0].header['NAXIS1'], hdulist[0].header['NAXIS2']
    in_ph_size_x, in_ph_size_y = hdulist[0].header['CDELT1']*NAXIS1, hdulist[0].header['CDELT2']*NAXIS2
    num_pix_x, num_pix_y = np.trunc(in_ph_size_x/delta_pix).astype(np.int), np.trunc(in_ph_size_y/delta_pix).astype(np.int)
            
    inwl   = np.zeros(len(hdulist))
    redata = np.zeros((num_pix_y, num_pix_x, len(hdulist)))
        
    xin = np.linspace(-1.0, 1.0, NAXIS1)
    yin = np.linspace(-1.0, 1.0, NAXIS2)

    xout = np.linspace(-1.0, 1.0, num_pix_x)
    yout = np.linspace(-1.0, 1.0, num_pix_y)
      
    for i, hdu in enumerate(hdulist):
        inwl[i]   = np.float64(hdu.header['WAV'])        
        f = interpolate.RectBivariateSpline(xin, yin, hdu.data)
        redata[..., i] = f(yout, xout)
        
    x = np.zeros((redata[...,i].shape[0],redata[...,i].shape[1], len(WavRange)))
    for i in range (x.shape[2]):
        x[...,i] = redata[...,0]
        
#==============================================================================
#Trim to square shape    
#==============================================================================        
    xmax,ymax = np.unravel_index(np.argmax(x[...,0], axis=None), x[...,0].shape)
    if x[...,0].shape[1]> x[...,0].shape[0]:
        sub = -5+x[...,0].shape[0]/2
        x = x[xmax-sub: xmax+sub, ymax-sub: ymax+sub]
    elif x[...,0].shape[1]< x[...,0].shape[0]:
        sub = -5+x[...,0].shape[1]/2
        x = x[xmax-sub: xmax+sub, ymax-sub: ymax+sub]
        
    for i in range (x.shape[2]):
        x[...,i] = x[...,0]/  x[...,0].sum()        
    
#    import matplotlib.pyplot as plt
#    plt.figure('psf 1')
#    plt.imshow(x[...,0])
#    plt.figure('psf 2')      
#    plt.imshow(x[...,-1])
 

#==============================================================================
#     PSF array must fit into the fpn
#==============================================================================

    if fpn[0] <= x[...,0].shape[0]/ osf  or fpn[0] <= x[...,0].shape[0]/ osf:
         print ('Trimming PSF array to fit FPA')
         sub =  int( ( int(( (x.shape[0]/ osf) -fpn[0] ) /2 ) + 5  ) * osf)
         x = x[sub:-sub, sub: -sub, :]
         

#    plt.figure('psf 1a')
#    plt.imshow(x[...,0])
#    plt.figure('psf 2a')
#    
#    plt.imshow(x[...,-1])
  
#==============================================================================
# normalise to 1   
#==============================================================================
    for i in range (x.shape[2]):
        x[...,i] = x[...,i]/  x[...,i].sum() 
        
#    plt.figure('psf 1b')
#    plt.imshow(x[...,0])
  
    
    return x
    

def Psf_spectrometer(zfile, delta_pix, osf, fpn,  WavRange):
    delta_pix = delta_pix.to(u.m)
    delta_pix =  delta_pix.value    
    hdul = fits.open(zfile)
    img_list = []
    shape_x =[]
    shape_y = []
    xmax_list = []; ymax_list = []
    wl_list = []
    img_list2=[]
    
#==============================================================================
#   Convert from pixels in original image to subpixels in exosim, 
#==============================================================================
    
    for i in range(1,len(hdul)):
#        print i
        Wav =  hdul[i].header['WAV'] 
       
        NAXIS1, NAXIS2 = hdul[i].header['NAXIS1'], hdul[i].header['NAXIS2']
        in_ph_size_x, in_ph_size_y = hdul[i].header['CDELT1']*NAXIS1, hdul[i].header['CDELT2']*NAXIS2
        num_pix_x, num_pix_y = np.trunc(in_ph_size_x/delta_pix).astype(np.int), np.trunc(in_ph_size_y/delta_pix).astype(np.int)

        redata = np.zeros((num_pix_y, num_pix_x))

        xin = np.linspace(-1.0, 1.0, NAXIS1)
        yin = np.linspace(-1.0, 1.0, NAXIS2)

        xout = np.linspace(-1.0, 1.0, num_pix_x)
        yout = np.linspace(-1.0, 1.0, num_pix_y)
    
        f = interpolate.RectBivariateSpline(xin, yin, hdul[i].data)
        
        redata = f(yout, xout)
        redata /= redata.sum()
        ymax, xmax = np.unravel_index(np.argmax(redata, axis=None), redata.shape)
        img_list.append(redata)
        img_list2.append(hdul[i].data)
        shape_x.append(redata.shape[1]) ; shape_y.append(redata.shape[0])  
        xmax_list.append(xmax); ymax_list.append(ymax)
        wl_list.append(Wav)
        
#        print redata.shape,xmax,ymax

#    plt.figure ('NS 1')
#    plt.imshow(img_list[2])
      
    idx1 = np.argmin(shape_x)
    idx2 = np.argmin(shape_y)
    shape_min_x = shape_x[idx1]; shape_min_y = shape_y[idx2]
#    print idx1, idx2
    
    
#==============================================================================
# Must centre maximum at same point on all images and cut to same array size prior to interpolation
#==============================================================================
    img_array = np.zeros((shape_min_y, shape_min_x, len(wl_list)))
    for i in range(len(img_list)):        
        x = img_list[i]
        xmax,ymax = np.unravel_index(np.argmax(x, axis=None), x.shape)
#        print x.shape,xmax, ymax, xmax_list[i]
        x = x[ymax_list[i]-ymax_list[idx2]: ymax_list[i]-ymax_list[idx2]+shape_y[idx2], xmax_list[i]-xmax_list[idx1]: xmax_list[i]-xmax_list[idx1]+shape_x[idx1]]
#        xmax,ymax = np.unravel_index(np.argmax(x, axis=None), x.shape)    
#        print x.shape, xmax, ymax
#        print ""
        img_array[...,i] = x
 
     
#==============================================================================
#Trim to square shape    
#==============================================================================        
    x = img_array   
    xmax,ymax = np.unravel_index(np.argmax(x[...,0], axis=None), x[...,0].shape)
    if x[...,0].shape[1]> x[...,0].shape[0]:
        sub = -5+x[...,0].shape[0]/2
        x = x[xmax-sub: xmax+sub, ymax-sub: ymax+sub]
    elif x[...,0].shape[1]< x[...,0].shape[0]:
        sub = -5+x[...,0].shape[1]/2
        x = x[xmax-sub: xmax+sub, ymax-sub: ymax+sub]
        
    for i in range (x.shape[2]):
        x[...,i] = x[...,i]/  x[...,i].sum()        
    
#    plt.figure('psf 1')
#    plt.imshow(x[...,0])
#    plt.figure('psf 2')      
#    plt.imshow(x[...,-1])
 
#    print "PSF 1 sum", np.sum(x[...,0])
#    print "Psf shape",  x.shape
    img_array = x
#    
#    maxx =[]
#    for i in range(img_array.shape[2]):   
#       maxx.append((img_array[...,i]/img_array[...,i].sum()).max())
#   
#   
#    plt.figure(888)
#    plt.plot(wl_list, maxx)       
#    xxxx
#
#   
  
#==============================================================================
# interpolate to x_wav_osr
#==============================================================================
#    psf_stack = interpolate.interp1d(wl_list, img_array, axis=2, bounds_error=False, fill_value=0.0, kind='quadratic')(WavRange)
    psf_stack = interpolate.interp1d(wl_list, img_array, axis=2, bounds_error=False, fill_value=0.0, kind='linear')(WavRange)

#==============================================================================
#     PSF array must fit into the fpn
#==============================================================================
    print ("PSF array width in whole pixel", x.shape[0]/ osf)
    print ("FPA size", fpn)
#    print psf_stack[...,0].shape
    if fpn[0] <= psf_stack[...,0].shape[0]/ osf:
         print ('Trimming PSF array to fit FPA')
         sub =  int( ( int(( (psf_stack[...,0].shape[0]/ osf) -fpn[0] ) /2 ) + 5  ) * osf)
         psf_stack = psf_stack[sub:-sub, sub: -sub, :]
    print ("PSF Stack shape", psf_stack.shape)
#    plt.figure('psf 1b')
#    plt.imshow(psf_stack[...,200]) 
    
    print ("pSF stack 200 sum",  psf_stack[...,200].sum())

#==============================================================================
# smoothing function    moving average of 9 - needed since rect spline intorduces lack of smooth PSFs
#==============================================================================
    psf_stack0 =psf_stack*1
    box = 9 #must be odd - could increase to 29
    a = int(box/2)
    b = a+1
    if a+b != box:
        print ("error")
        xxx

    for i in range (a, psf_stack.shape[2]-b):
      
        psf_stack[...,i]=  psf_stack0[...,i-a:i+b].sum(axis=2) /(a+b)
#        print psf_stack[...,i].sum()

    for i in range (psf_stack.shape[2]):
        if np.sum(psf_stack[...,i])!=0:
            psf_stack[...,i] = psf_stack[...,i]/np.sum(psf_stack[...,i])

    print ("pSF stack std",  psf_stack.std())
 
    return psf_stack
    

    
 
def Psf(wl, fnum_x, fnum_y, delta, nzero = 4, shape='airy'):
  '''
  Calculates an Airy Point Spread Function arranged as a data-cube. The spatial axies are 
  0 and 1. The wavelength axis is 2. Each PSF area is normalised to unity.
  
  Parameters
  ----------
  wl	: ndarray [physical dimension of length]
    array of wavelengths at which to calculate the PSF
  fnum : scalar
    Instrument f/number
  delta : scalar
    the increment to use [physical units of length]
  nzero : scalar
    number of Airy zeros. The PSF kernel will be this big. Calculated at wl.max()
  shape : string
    Set to 'airy' for a Airy function,to 'gauss' for a Gaussian
  
  Returns
  ------
  Psf : ndarray
    three dimensional array. Each PSF normalised to unity
  '''
#  fnum_y =  fnum_x
  delta = delta.rescale(wl.units)
  Nx = int(np.round(scipy.special.jn_zeros(1, nzero)[-1]/(2.0*np.pi) * fnum_x*wl.max()/delta).astype(np.int))
  
  Ny = Nx = np.int(Nx)
#  Ny = int(np.round(scipy.special.jn_zeros(1, nzero)[-1]/(2.0*np.pi) * fnum_y*wl.max()/delta).astype(np.int))

  if shape=='airy':
    d = 1.0/(fnum_x*(1.0e-30*delta.units+wl))
  elif shape=='gauss':
    sigma = 1.029*fnum_x*(1.0e-30*delta.units+wl)/np.sqrt(8.0*np.log(2.0))
    d     = 0.5/sigma**2
    
  x = np.linspace(-Nx*delta.item(), Nx*delta.item(), 2*Nx+1)*delta.units
  y = np.linspace(-Ny*delta.item(), Ny*delta.item(), 2*Ny+1)*delta.units
  
  yy, xx = np.meshgrid(y, x)
 
  if shape=='airy':
    arg = 1.0e-20+np.pi*np.multiply.outer(np.sqrt(yy**2 + xx**2), d)
    img   = (scipy.special.j1(arg)/arg)**2
  elif shape=='gauss':
    arg = np.multiply.outer(yy**2 + xx**2, d)
    img = np.exp(-arg)
  
  if fnum_y !=  fnum_x:
      x_pix = img.shape[0]  
      stretch = fnum_y/fnum_x
      y_pix = int(np.round(x_pix*stretch))
      
      x_pos = np.linspace(0,1,x_pix)
      y_pos = np.linspace(0,1,y_pix)
      new_img = np.zeros((y_pix, x_pix, img.shape[2]))
#      print x_pix, img[...,0].shape
      for i in range(img.shape[2]):
          img_ = interpolate.interp2d(x_pos,x_pos,img[...,i], kind='linear')(x_pos,y_pos)
          new_img[...,i]=img_
#      plt.figure('img')
#      plt.imshow(img[...,1000])      
#      plt.figure('img2')
#      plt.imshow(new_img[...,1000])
      img = new_img

  norm = img.sum(axis=0).sum(axis=0)
  img /= norm
  
  idx = np.where(wl <= 0.0)
  if idx:
    img[..., idx] *= 0.0
    
  print (img.shape)

  
  return img
   
   
   
   
def Psf_Twinkle(WavRange):
    
    from scipy.io import loadmat
    import matplotlib.pyplot as plt

    a = ['/Users/user1/Desktop/Twinkle_PSF_SWIR_STR.mat' ,'/Users/user1/Desktop/Twinkle_PSF_MWIR_STR.mat']
    x1 = loadmat(a[0])
    wav1 = x1['lamb_save'][0]
    x2 = loadmat(a[1])
    wav2 = x2['lamb_save'][0]
    p2 = x2['PSF_nat2'][0:-1,0:-1,...] # to make same shape - it's okay since max is on the same pixel
    
    wav_list =[]
     
    img_stack = np.zeros((p2.shape[0]*3, p2.shape[1]*3, len(wav1)+len(wav2)))
#    img_stack2 = np.zeros((p2.shape[0], p2.shape[1], len(wav1)+len(wav2)))

    #A,B = np.unravel_index(p2[...,0].argmax(), p2[...,0].shape)
    #
    #print A,B
    
    #print img_stack.shape
    
    ct=0 
    for ii in range (len(a)):
        x = loadmat(a[ii])
        ## one-liner to read a single variable
    #    p1 = x['PSF_fgs2']
        p2 = x['PSF_nat2']
        if ii == 1:
            p2 = x['PSF_nat2'][0:-1,0:-1,...]
        
        wav = x['lamb_save'][0]
        for j in range(len(wav)):
            wav_list.append(wav[j])
            xin = np.linspace(-1.0, 1.0, p2[...,j].shape[1])
            yin = np.linspace(-1.0, 1.0, p2[...,j].shape[0])
    
            xout = np.linspace(-1.0, 1.0, p2[...,j].shape[1]*3)
            yout = np.linspace(-1.0, 1.0, p2[...,j].shape[1]*3)   
            
            f = interpolate.interp2d(xin, yin, p2[...,j], kind='linear', copy=True, bounds_error=False, fill_value=None)
            redata = f(yout, xout)
            
            img_stack[...,ct] = redata/redata.sum() #normalise
    #        img_stack2[...,ct] = p2[...,j]/p2[...,j].sum()
    
            ct+=1
            
    #for i in range (len(wav_list)):
    #    plt.figure('%s'%(wav_list[i]))
    #    plt.imshow(img_stack[...,i])          
            
    psf_stack = interpolate.interp1d(wav_list, img_stack, axis=2, bounds_error=False, fill_value=0.0, kind='linear')(WavRange)
    #psf_stack2 = interpolate.interp1d(wav_list, img_stack2, axis=2, bounds_error=False, fill_value=0.0, kind='linear')(WavRange)
    
    
    #for i in range (len(WavRange)):
    #    plt.figure('%s'%(WavRange[i]))
    #    plt.imshow(psf_stack[...,i])
    #    print np.sum(psf_stack[...,i])
        
    #for i in range (len(WavRange)):
    #    plt.figure('%s2'%(WavRange[i]))
    #    plt.imshow(psf_stack2[...,i])    
    #    
    #    
        
   
    return psf_stack
   
   
   
  
def PixelResponseFunction(opt, psf_shape, osf, delta, lx = 1.7*u.um, ipd = 0.0*u.um):
  '''
  Estimate the detector pixel response function with the prescription of 
  Barron et al., PASP, 119, 466-475 (2007).
  
  Parameters
  ----------
  psf_shape	: touple of scalars 
		  (ny, nx) defining the PSF size	
  osf		: scalar
		  number of samples in each resolving element. The 
		  final shape of the response function would be shape*osf
  delta 	: scalar
		  Phisical size of the detector pixel in microns
  lx		: scalar
		  diffusion length in microns
  ipd           : scalar
		  distance between two adjacent detector pixels 
		  in microns
		 
  Returns
  -------
  kernel	: 2D array
		  the kernel image
  kernel_delta  : scalar
                  the kernel sampling interval in microns
  '''
  if type(osf) != int: osf = np.int(osf)
    
#  lx = 0*u.um # top hat 
  lx += 1e-8*u.um # to avoid problems if user pass lx=0
#==============================================================================
#  lx = 3.7*u.um # approximates Hardy et al
#==============================================================================
  lx = lx.to(delta.unit)
 
  exosim_msg ("diffusion length in IPRF %s"%(lx), opt.diagnostics)
 
#==========FOR IMAGES ONLY====================================================================
#  osf*=33 # for image demo only 
#==============================================================================
  kernel = np.zeros( (psf_shape[0]*osf, psf_shape[1]*osf) )
  
#===========FOR IMAGES ONL===================================================================
#  kernel = np.zeros( (3*osf, 3*osf) )   # for image demo only
#==============================================================================
    
  kernel_delta = delta/osf
  yc, xc = np.array(kernel.shape) // 2
  yy = (np.arange(kernel.shape[0]) - yc) * kernel_delta 
  xx = (np.arange(kernel.shape[1]) - xc) * kernel_delta 
  mask_xx = np.where(np.abs(xx) > 0.5*(delta-ipd))
  mask_yy = np.where(np.abs(yy) > 0.5*(delta-ipd))
  xx, yy = np.meshgrid(xx, yy)

  kernel = np.arctan(np.tanh( 0.5*( 0.5*delta.value - xx.value)/lx.value )) - \
	   np.arctan(np.tanh( 0.5*(-0.5*delta.value - xx.value)/lx.value ))
	 	 
  kernel*= np.arctan(np.tanh( 0.5*( 0.5*delta.value - yy.value)/lx.value )) - \
  	   np.arctan(np.tanh( 0.5*(-0.5*delta.value - yy.value)/lx.value )) 

# for cross-talk (e.g if bell shaped function like Hardy) comment out these lines below.
  kernel[mask_yy, ...] = 0.0  
  kernel[..., mask_xx] = 0.0
#
#  # Normalise the kernel such that the pixel has QE=1
  kernel *= osf**2/kernel.sum()
  
 
  
#==============================================================================
# Below is for images only  
# 
#   
#  from mpl_toolkits.mplot3d import Axes3D
#  import matplotlib.pyplot as plt
#  from matplotlib import cm
#  from matplotlib.ticker import LinearLocator, FormatStrFormatter
#
#
#  fig = plt.figure('IPRF')
#  ax = fig.gca(projection='3d')
# 
#  surf = ax.plot_surface(xx/18e-6, yy/18e-6, kernel, cmap=cm.coolwarm,
#                       linewidth=0, antialiased=False)
#                       
#  ax.set_xlim(-1.5, 1.5)
#  ax.set_ylim(-1.5, 1.5)
#  ax.set_zlim(0, 1.0)
#  ax.set_xlabel('X axis (pixels)')
#  ax.set_ylabel('Y axis (pixels)')
#  ax.set_zlabel('Relative count per pixel')                     
# 
#  fig.colorbar(surf, shrink=0.5, aspect=5)
#
#  plt.show()  
#  
#  plt.figure('X-section through IPRF')
#  X,Y = np.unravel_index(kernel.argmax(), kernel.shape)
#  
#  import pandas as pd
#  qq = np.array(pd.read_csv('/Users/user1/Desktop/JWST_IRPRF_Hardy.csv'))
#  qq = qq[qq[:,0].argsort()]
# 
#  xx0 = np.linspace(3,8,1000)
# 
#  yy0  = np.interp(xx0,qq[:,0],qq[:,1])
#  idx = np.argwhere(yy0>=0)
#  yy0= yy0[idx]
#  xx0=xx0[idx]
#  idx = np.argmax(yy0)
#  xx0 =xx0-xx0[idx]
#  yy0 = yy0/yy0.max()
#  plt.plot(xx0+0.012,yy0, 'r-', linewidth=2, label='Hardy et al. 2014')
#  plt.plot(xx[0]/18e-6, kernel[X]/kernel[X].max(), 'b--',linewidth=2, label = 'JexoSim')
#  plt.grid(True)
#  plt.xlabel('Distance (pixels)')
#  plt.ylabel('Relative response')
#  
#  ax = plt.gca() 
#
#
#  legend = ax.legend(loc='upper right', shadow=True)
#  frame = legend.get_frame()
#  frame.set_facecolor('0.90')
#  for label in legend.get_texts():
##    label.set_fontsize('medium')
#    label.set_fontsize(22)
#  for label in legend.get_lines():
#    label.set_linewidth(1.5)  # the legend line width
#  for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
#             ax.get_xticklabels() + ax.get_yticklabels()):
#    item.set_fontsize(22)
#    
#  xxx
 
  
#==============================================================================
 
  kernel = np.roll(kernel, -xc, axis=1)
  kernel = np.roll(kernel, -yc, axis=0)
  
  return kernel, kernel_delta


def pointing_add_scan(pointing_timeline, scan_throw_arcsec, frame_time, frame_osf, exposure_time):  
  ''' Superimpose saw-tooth scan mode to a pointing jitter timeline.
  The period of a scan is equal exposure time
  
  Parameters
  ----------
  pointing_timeline: Quantities Array 
         Poitning timeline (yaw/pitch) in deg
  scan_throw_arcsec: scalar
         The scan throw in units of arcseconds. 
  frame_time: scalar
      detector frame time in units of time
  frame_osf: scalar
      Frame oversampling factor
  exposure_time: scalar
      time for one exposure containing set of NDRs
  Returns
  -------
  pointing_timeline: Quantities Array
      Pointing timeline  updated with saw-tooth scan pointing superimposed. 
  '''

  tArr = np.arange(len(pointing_timeline) +1)*frame_time/frame_osf   
  tArr = tArr[1:]

  saw_tooth = (tArr.value %exposure_time.value /exposure_time.value)    
  saw_tooth = np.where(saw_tooth > 0.0, saw_tooth, 1.0)
  saw_tooth = saw_tooth*  scan_throw_arcsec.rescale(u.deg)
  pointing_timeline = saw_tooth + pointing_timeline 
        
  return pointing_timeline   
  


def pointing_add_drift(pointing_timeline, drift_rate, frame_time, frame_osf, exposure_time):  
  ''' Superimpose saw-tooth scan mode to a pointing jitter timeline.
  The period of a scan is equal exposure time
  
  Parameters
  ----------
  pointing_timeline: Quantities Array 
         Poitning timeline (yaw/pitch) in deg
  drift_rate: scalar
         The drift rate in units of arcseconds/sec
  frame_time: scalar
      detector frame time in units of time
  frame_osf: scalar
      Frame oversampling factor
  exposure_time: scalar
      time for one exposure containing set of NDRs
  Returns
  -------
  pointing_timeline: Quantities Array
      Pointing timeline  updated with saw-tooth scan pointing superimposed. 
  '''

  print ("adding drift to jitter......")

  tArr = np.arange(len(pointing_timeline) +1)*frame_time/frame_osf   
  tArr = tArr[1:]
  
  drift_rate = drift_rate*u.arc_second/u.s
  drift_rate = drift_rate.rescale(u.deg/u.s)
  tArr = tArr
  
  saw_tooth = tArr*drift_rate

  plt.figure('drift')
  plt.plot(tArr, saw_tooth)
  
  plt.figure('pointing timeline')
  plt.plot(tArr, pointing_timeline)  
  
  pointing_timeline = saw_tooth + pointing_timeline 
    
  plt.figure('pointing timeline + drift')
  plt.plot(tArr, pointing_timeline)  
  
  return pointing_timeline   


def pointing_jitter(opt):
    
  jitter_file = opt.jitter_psd_file 
  total_observing_time = opt.total_observing_time
  frame_time = opt.frame_time
  rms = opt.simulation.sim_pointing_rms.val
    
  ''' Estimate pointing jitter timeline
  
  Parameters
  ----------
  jitter_file: string
	       filename containing CSV columns with 
	       frequency [Hz], Yaw PSD [deg**2/Hz], Pitch [deg**2/Hz]
	       If only two columns given, then it is assumed that 
	       the second column is the PSD of radial displacements
  totoal_observing_time: scalar
      total observing time in units of time
  frame_time: scalar
      detector frame time in units of time
  rms: scalar
      renormalisation rms in units of angle
      
  Returns
  -------
  yaw_jit: jitter timeline in units of degrees
  pitch_jit: jitter rimeline in units of degrees
  osf: number of additional samples in each frame_time needed to capture
       the jitter spectral information
  '''
  
  data = np.genfromtxt(jitter_file, delimiter=',')
  psd_freq = data[..., 0]

  if data.shape[1] > 2:
    psd_yaw = data[..., 1]
    psd_pitch = data[..., 2]
  else:
    psd_yaw = data[..., 1]/2
    psd_pitch = psd_yaw

  # each frame needs to be split such that jitter is Nyquis sampled
  jitter_sps = 2.0*psd_freq.max()/u.s
  
  osf = (np.ceil((frame_time).to(u.s) * jitter_sps).take(0).astype(np.int)).value

  if osf < 1: osf = 1
      
  exosim_msg ("Frame OSF (number of jitter frames per NDR %s"%(osf) , opt.diagnostics)

  number_of_samples_ = np.int(osf*np.ceil(total_observing_time/frame_time))  + 100
  N0 = number_of_samples_ 

  number_of_samples = 2**(np.ceil(np.log2(number_of_samples_)))
  number_of_samples = int( (number_of_samples /2)+1  )

  freq = np.linspace(0.0, 0.5* osf/ (frame_time).to(u.s).value, number_of_samples)

  npsd_yaw   = 1.0e-30+np.interp(freq, psd_freq, psd_yaw, left=0.0, right=0.0)
  npsd_pitch = 1.0e-30+np.interp(freq, psd_freq, psd_pitch, left=0.0, right=0.0)

  exosim_msg ("interpolation of PSD done", opt.diagnostics)
  # import matplotlib.pyplot as plt
  # plt.figure(33)
  # plt.plot(psd_freq,psd_yaw, 'rx-')
  # plt.plot(freq,npsd_yaw,'bx-')
  
  # xxxx
  npsd_yaw    = np.sqrt(npsd_yaw   * np.gradient(freq))
  npsd_pitch  = np.sqrt(npsd_pitch * np.gradient(freq))

  yaw_jit_re   = np.random.normal(scale=npsd_yaw/2.0)
  yaw_jit_im   = np.random.normal(scale=npsd_yaw/2.0)
  pitch_jit_re = np.random.normal(scale=npsd_pitch/2.0)
  pitch_jit_im = np.random.normal(scale=npsd_pitch/2.0)
  
  pitch_jit_im[0] = pitch_jit_im[-1] = 0.0
  yaw_jit_im[0]   = yaw_jit_im[-1]   = 0.0

  norm = 2*(number_of_samples-1)
 
  exosim_msg ("starting irfft" , opt.diagnostics)
  yaw_jit = norm*np.fft.irfft(yaw_jit_re + 1j * yaw_jit_im)*u.deg
  pitch_jit = norm*np.fft.irfft(pitch_jit_re + 1j * pitch_jit_im)*u.deg
  exosim_msg ("completed.....", opt.diagnostics)
  
  if opt.simulation.sim_adjust_rms.val ==1:
      
    norm = ((rms/3600.0/1000.0)**2/((yaw_jit.value).var()+ (pitch_jit.value).var()))
    yaw_jit *= np.sqrt(norm)
    pitch_jit *= np.sqrt(norm)
 
  if len(yaw_jit) > N0:
      yaw_jit = yaw_jit[0:N0]
      pitch_jit = pitch_jit[0:N0] 

  exosim_msg ("jitter RMS in mas %s %s"%(np.std (yaw_jit.value)*3600*1000, np.std(pitch_jit.value)*3600*1000) , opt.diagnostics)
    
  return yaw_jit, pitch_jit, osf
 
def oversample(fp, ad_osf):
    
    xin = np.linspace(0,fp.shape[1]-1,fp.shape[1])
    yin = np.linspace(0,fp.shape[0]-1,fp.shape[0])
    x_step =  abs(xin[1]) - abs(xin[0])
    y_step =  abs(yin[1]) - abs(yin[0])
    
    # calculates the new step sizes for new grid
    x_step_new = np.float(x_step/ad_osf)
    y_step_new = np.float(y_step/ad_osf)
    
    # new grid must start with an exact offset to produce correct number of new points
    x_start = -x_step_new * np.float((ad_osf-1)/2)
    y_start = -y_step_new * np.float((ad_osf-1)/2)
    
    # new grid points- with correct start, end and spacing
    xout = np.arange(x_start, x_start + x_step_new*fp.shape[1]*ad_osf, x_step_new)
    yout = np.arange(y_start, y_start + y_step_new*fp.shape[0]*ad_osf, y_step_new)
    
    # interpolate fp onto new grid
    fn = interpolate.RectBivariateSpline(yin,xin, fp)
    new_fp = fn(yout,xout)
    
    return new_fp


def animate(Data):
    
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

#    Data = data['channel']['SWIR'].timeline

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    wframe = None
    
    for j in range(0,Data.shape[2]):
    
            oldcol = wframe

            X = np.arange(0, Data.shape[1])
            Y = np.arange(0, Data.shape[0])
            X, Y = np.meshgrid(X, Y)
            
           
        
            Z = Data[...,j]
#            print Z.sum()
            
            wframe = ax.plot_wireframe(X, Y, Z, rstride=2, cstride=2)
            
#            ax.set_zlim(0,20000)
#            ax.set_title(j)
#            
     
        

    
        # Remove old line collection before drawing
            if oldcol is not None:
                ax.collections.remove(oldcol)
    
            plt.pause(.01)
            
            
            
def write_record(opt, path, lab, input_text_file):
    
    textfile = '%s/%s.txt'%(path,lab)
    
    with open(textfile, "w") as f1:
        f1.write('===== Simulation values =====')
        f1.write('\n ')
        f1.write('\nPlanet:  %s'%(opt.planet.planet.name))
        f1.write('\nChannel:  %s'%(opt.channel.name))
        if opt.simulation.sim_mode.val == 3:
            f1.write('\n\nNoise option:  noise budget')          
        else:
            f1.write('\n\nNoise option:  %s'%(opt.noise_tag))     
        f1.write('\n ')

        f1.write('\nUse saturation time?:  %s'%(opt.observation.obs_use_sat.val))
        f1.write('\nSat time (to designated fraction of full well):  %s sec'%(opt.sat_time))
        f1.write('\nt_f:  %s sec'%(opt.t_f.value))
        f1.write('\nt_g:  %s sec'%(opt.t_g.value))
        f1.write('\nt_sim:  %s sec'%(opt.t_sim.value))
        
        f1.write('\nsaturation flag  :  %s'%(opt.sat_flag))
        f1.write('\nsaturation limit  :  %s electrons'%(opt.sat_limit.value))
        f1.write('\nnumber of saturated pixels per image  :  %s'%(opt.no_sat))
        f1.write('\nzero values applied to saturated pixels  :  %s'%(opt.pipeline.pipeline_bad_corr.val))
         
        f1.write('\nt_int:  %s sec'%(opt.t_int.value.item()))
        f1.write('\nt_cycle:  %s sec'%(opt.exposure_time.value.item()))  
        f1.write('\nprojected multiaccum (n groups):  %s'%(opt.projected_multiaccum))
        f1.write('\neffective multiaccum (n groups):  %s'%(opt.effective_multiaccum))
        f1.write('\nnumber of NDRs simulated:  %s'%(opt.n_ndr) )
        f1.write('\nnumber of integration cycles:  %s'%(opt.n_exp) )          
        f1.write('\n ')
        if opt.simulation.sim_output_type.val == 1: # excludes fits files
            f1.write('\nApFactor:  %s'%(opt.pipeline.pipeline_ap_factor.val) )
            f1.write('\nAperture shape:  %s'%(opt.pipeline.pipeline_ap_shape.val) )
            f1.write('\nSpectral binning:  %s '%(opt.pipeline.pipeline_binning.val) )
            if opt.pipeline.pipeline_binning.val == 'R-bin':
                        f1.write('\nBinned R power:  %s '%(opt.channel.pipeline_params.pipeline_R.val) )
            else:
                      f1.write('\nBin size (pixels):  %s '%(opt.channel.pipeline_params.pipeline_bin_size.val) )
        f1.write('\nWavelength: %s %s'%(opt.channel.pipeline_params.wavrange_lo.val, opt.channel.pipeline_params.wavrange_hi.val) )
    
    with open(textfile, "a") as f1:
        f1.write('\n ')
        f1.write('\n ')
        f1.write('===== Copy of input parameters file used =====')
        f1.write('\n ')
        f1.write('\n ')
      
    with open(input_text_file) as f:
        with open(textfile, "a") as f1:
            for line in f:       
                f1.write(line)

    f1.close()
  