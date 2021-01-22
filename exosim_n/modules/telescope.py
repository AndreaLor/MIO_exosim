"""
exosim_n

Telescope module

"""

from exosim_n.classes.sed import Sed
from exosim_n.lib import exosim_n_lib, instrument_lib
from exosim_n.lib.exosim_n_lib import exosim_n_msg
import numpy as np
from astropy import units as u
import copy 

def run(opt):
    
      opt.osf         = np.int(opt.channel.osf())
      opt.offs        = np.int(opt.channel.pix_offs())    
      fpn = opt.channel.array_geometry.val.split(',')
      opt.fpn = [int(fpn[0]), int(fpn[1])]  
      opt.fp  = np.zeros(( int(opt.fpn[0]*opt.channel.osf.val),
                            int(opt.fpn[1]*opt.channel.osf.val) ))    
      opt.fp_signal  = np.zeros(( int(opt.fpn[0]*opt.channel.osf.val),
                            int(opt.fpn[1]*opt.channel.osf.val) ))            
      opt.fp_delta = opt.channel.detector_pixel.pixel_size.val / opt.channel.osf.val                     
#      opt.x_wav_osr, opt.x_pix_osr, opt.y_pos_osr = instrument_lib.usePoly(opt)
      if opt.channel.is_spec.val == True:
          opt.x_wav_osr, opt.x_pix_osr, opt.y_pos_osr = instrument_lib.useInterp(opt)
      elif opt.channel.is_spec.val == False:
          opt.x_wav_osr, opt.x_pix_osr,opt.y_pos_osr, opt.phot_wav  =  instrument_lib.getPhotometerWavelengths(opt)
          
      exosim_n_msg('tel check 1: %s'%(opt.star.sed.sed.max()), opt.diagnostics)
      opt.star_sed = copy.deepcopy(opt.star.sed) #copy of star flux at telescope
      opt.star_sed2 = copy.deepcopy(opt.star.sed) #copy of star flux at telescope
      opt.Aeff = 0.25*np.pi*opt.common_optics.telescope_effective_diameter()**2 
      opt.star.sed.sed*= opt.Aeff
      exosim_n_msg('tel check 2: %s'%(opt.star.sed.sed.max()), opt.diagnostics)
      tr_ =np.array([1.]*len(opt.x_wav_osr))*u.dimensionless_unscaled
      opt.common_optics.optical_surface = opt.common_optics.optical_surface \
          if isinstance(opt.common_optics.optical_surface, list) \
          else [opt.common_optics.optical_surface]    
      for op in opt.common_optics.optical_surface:
          dtmp=np.loadtxt(op.transmission.replace('__path__', opt.__path__), delimiter=',')
          tr = Sed(dtmp[:,0]*u.um,dtmp[:,1]*u.dimensionless_unscaled)
          tr.rebin(opt.x_wav_osr)  
          tr_ *=tr.sed       
      opt.telescope_transmission = Sed(opt.x_wav_osr, tr_) 
                  
      opt.star.sed.rebin(opt.x_wav_osr)  
      exosim_n_msg('tel check 2a: %s'%(opt.star.sed.sed.max()), opt.diagnostics)
  
      opt.planet.sed.rebin(opt.x_wav_osr) 
      exosim_n_lib.sed_propagation(opt.star.sed, opt.telescope_transmission)
      exosim_n_msg('tel check 3: %s'%(opt.star.sed.sed.max()), opt.diagnostics)
      

      return opt
  
