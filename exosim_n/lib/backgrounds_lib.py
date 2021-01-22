"""
exosim_n

Backgrounds library

"""

import numpy as np
from exosim_n.lib import exosim_n_lib
from exosim_n.lib.exosim_n_lib import exosim_n_msg, planck
from astropy import units as u
from exosim_n.classes.sed import Sed


def zodical_light(opt):
   
    deg = opt.model_exosystem.ecliptic_lat.val.value
    wl = opt.x_wav_osr
    exosim_n_msg('Zodi model initiated...\n', opt.diagnostics)
    exosim_n_msg('Ecliptic latitude in deg: %s   \n'%(deg), opt.diagnostics)
    deg = abs(deg)
    if deg >=57.355:
        level = 1.0
    else:
        level = -0.22968868*(np.log10(deg+1))**7 + 1.12162927*(np.log10(deg+1))**6 - 1.72338015*(np.log10(deg+1))**5 + 1.13119022*(np.log10(deg+1))**4 - 0.95684987*(np.log10(deg+1))**3+ 0.2199208*(np.log10(deg+1))**2- 0.05989941*(np.log10(deg+1))  +2.57035947
    exosim_n_msg('Zodi model coefficient... %s\n'%(level), opt.diagnostics)
    spectrum = level*(3.5e-14*planck(wl, 5500*u.K) +
                      planck(wl, 270*u.K) * 3.58e-8)
    sed = Sed(wl, spectrum )
    transmission = Sed(wl, np.ones(wl.size)*u.dimensionless_unscaled)
    zodi = [sed, transmission]
    return zodi
    
def emission_from_file(opt):
    
  instrument_emission  = Sed(opt.x_wav_osr, 
                             np.zeros(len(opt.x_wav_osr))* 
                             u.W/u.m**2/u.um/u.sr)
  instrument_transmission = Sed(opt.x_wav_osr, np.ones(len(opt.x_wav_osr)))
  opt.common_optics.optical_surface = opt.common_optics.optical_surface \
          if isinstance(opt.common_optics.optical_surface, list) \
          else [opt.common_optics.optical_surface]    
  for op in opt.common_optics.optical_surface:
      dtmp=np.loadtxt(op.transmission.replace('__path__', opt.__path__), delimiter=',')
      tr = Sed(dtmp[:,0]*u.um,dtmp[:,1]*u.dimensionless_unscaled)     
      tr.rebin(opt.x_wav_osr) 
      TR =np.ones(len(opt.x_wav_osr))
      TR =TR*tr.sed
      if opt.diagnostics==1:
          import matplotlib.pyplot as plt
          plt.figure('Total transmission after each optical surface added, OTE + channel %s'%(opt.channel.name))
          plt.plot(opt.x_wav_osr,TR)   
          plt.ylim(0,1.2)
        
          plt.figure('Transmission of each optical element added in the OTE %s'%(opt.channel.name))
          plt.plot(tr.wl, tr.sed) 
          plt.ylim(0,1.2)     
    
      em = Sed(dtmp[:,0]*u.um,dtmp[:,2]*u.dimensionless_unscaled)
      em.rebin(opt.x_wav_osr)
          
      exosim_n_lib.sed_propagation(instrument_emission, tr, emissivity=em, temperature=op())     
      instrument_transmission.sed = instrument_transmission.sed*tr.sed

  for op in opt.channel.optical_surface:
      dtmp=np.loadtxt(op.transmission.replace(
          '__path__', opt.__path__), delimiter=',')
      tr = Sed(dtmp[:,0]*u.um, dtmp[:,1]*u.dimensionless_unscaled)
      tr.rebin(opt.x_wav_osr)
      TR *=tr.sed
      if opt.diagnostics==1:
          plt.figure('Total transmission after each optical surface added, OTE + channel %s'%(opt.channel.name))
          plt.plot(opt.x_wav_osr,TR)   
          plt.ylim(0,1.2)
          
          plt.figure('Transmission of each optical element added in the channel %s'%(opt.channel.name))
          plt.plot(tr.wl, tr.sed) 
          plt.ylim(0,1.2) 
           
      em = Sed(dtmp[:,0]*u.um, dtmp[:,2]*u.dimensionless_unscaled)
      em.rebin(opt.x_wav_osr)

      exosim_n_lib.sed_propagation(instrument_emission, tr, emissivity=em,temperature=op())
                        
            
  emission=  instrument_emission
      
  return emission
    
    
  
 