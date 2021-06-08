"""
exosim_n
 
Backgrounds modules

"""


from exosim_n.lib.exosim_n_lib import exosim_n_msg, exosim_n_plot
from exosim_n.lib import exosim_n_lib, backgrounds_lib 
from exosim_n.classes.sed import Sed
import numpy  as np
from astropy import units as u
import scipy.interpolate
import scipy.signal
import copy

from exosim_n.classes.options import Options 
import matplotlib.pyplot as plt

def run(opt):
      DEBUG= Options.DEBUG
      if DEBUG:
        oldDiagnostic=opt.diagnostics
        opt.diagnostics=Options.showBackground
      exosim_n_msg('Running backgrounds module ...\n ', opt.diagnostics)
      
      opt.zodi, zodi_transmission  = backgrounds_lib.zodical_light(opt)
      opt.emission = backgrounds_lib.emission_from_file(opt)

      exosim_n_lib.sed_propagation(opt.star.sed, zodi_transmission)
 
      ch = opt.channel
      Omega_pix = 2.0*np.pi*(1.0-np.cos(np.arctan(0.5/ch.wfno())))*u.sr
      Apix = ((ch.detector_pixel.pixel_size.val).to(u.m))**2           
      opt.zodi.sed *= Apix*Omega_pix*opt.Re *u.electron/u.W/u.s 
      opt.emission.sed *= Apix*Omega_pix*opt.Re *u.electron/u.W/u.s 
      exosim_n_lib.sed_propagation(opt.zodi, opt.total_transmission)

      idx = np.argwhere(opt.x_wav_osr>0).T[0]
      wav = opt.x_wav_osr[idx]
      if np.argmin(wav) > np.argmax(wav):
          temp_zodi = Sed(opt.zodi.wl[::-1], opt.zodi.sed[::-1])
          temp_zodi.rebin(opt.x_wav_osr[::-1])
          opt.zodi.sed = temp_zodi.sed[::-1]
          opt.zodi.wl = opt.x_wav_osr
          temp_emission = Sed(opt.emission.wl[::-1], opt.emission.sed[::-1])
          temp_emission.rebin(opt.x_wav_osr[::-1])
          opt.emission.sed = temp_emission.sed[::-1]
          opt.emission.wl = opt.x_wav_osr
      else: 
          opt.zodi.rebin(opt.x_wav_osr)
          opt.emission.rebin(opt.x_wav_osr)  
      
      exosim_n_plot('zodi spectrum', opt.diagnostics, xdata = opt.zodi.wl, ydata=opt.zodi.sed)
      exosim_n_plot('emission spectrum', opt.diagnostics, xdata = opt.emission.wl, ydata=opt.emission.sed)   

      opt.zodi.sed     *= opt.d_x_wav_osr 
      opt.emission.sed *= opt.d_x_wav_osr
      zodi_photons_sed = copy.deepcopy(opt.zodi.sed)
      emission_photons_sed = copy.deepcopy(opt.emission.sed)

      # need to work on this to make more accurate : filter wavelength range might not correspond to wl sol on the detector, i.e. some wavelengths not being included that should be
      if opt.channel.is_spec.val == True:
          if ch.slit_width.val == 0:
              slit_size =  opt.fpn[1]*2  # to ensure all wavelengths convolved onto all pixels in slitless case
          else:
              slit_size =  ch.slit_width.val   
          opt.zodi.sed     = scipy.signal.convolve(opt.zodi.sed, 
    		      np.ones(np.int(slit_size*ch.osf())), 
    		      'same') * opt.zodi.sed.unit
          opt.emission.sed = scipy.signal.convolve(opt.emission.sed.value, 
    		      np.ones(np.int(slit_size*ch.osf())), 
    		      'same') * opt.emission.sed.unit
      
      elif opt.channel.is_spec.val == False:
          wav = np.repeat(opt.phot_wav,opt.fp.shape[1])
      
          opt.zodi.sed = np.repeat(opt.zodi.sed.sum(),
    					    wav.size)
          opt.zodi.wl = wav
          opt.emission.sed = np.repeat(opt.emission.sed.sum(),
    						wav.size)
          opt.emission.wl = wav
      
      opt.zodi_sed_original = copy.deepcopy(opt.zodi.sed)
      opt.emission_sed_original = copy.deepcopy(opt.emission.sed)    
      
      exosim_n_plot('zodi spectrum', opt.diagnostics, xdata = opt.zodi.wl, ydata=opt.zodi.sed)
      exosim_n_plot('emission spectrum', opt.diagnostics, xdata = opt.emission.wl, ydata=opt.emission.sed)   
      if DEBUG:
           if opt.diagnostics:
              plt.show()
           opt.diagnostics=oldDiagnostic
      return opt
      

    
  