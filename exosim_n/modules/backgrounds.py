"""
JexoSim
2.0
Backgrounds modules
v1.0

"""

from exosim_n.classes.sed import Sed
from exosim_n.lib.exolib import exosim_msg, exosim_plot
from exosim_n.lib import exolib 
import numpy           as np
from astropy import units as u
import scipy.interpolate
import scipy.signal
import copy
import matplotlib.pyplot as plt

def zodical_light(opt):
   
    deg = opt.model_exosystem.ecliptic_lat.val.value
    wl = opt.x_wav_osr
    exosim_msg('Zodi model initiated...\n', opt.diagnostics)
    # exosim_msg('Ecliptic latitude in deg: %s   \n'%(deg), opt.diagnostics)
    # deg = abs(deg)
 
    # if deg >=57.355:
    #     level = 1.0
    # else:
    #     level = -0.22968868*(np.log10(deg+1))**7 + 1.12162927*(np.log10(deg+1))**6 - 1.72338015*(np.log10(deg+1))**5 + 1.13119022*(np.log10(deg+1))**4 - 0.95684987*(np.log10(deg+1))**3+ 0.2199208*(np.log10(deg+1))**2- 0.05989941*(np.log10(deg+1))  +2.57035947
    # exosim_msg('Zodi model coefficient... %s\n'%(level), opt.diagnostics)
    level = 1.0
    spectrum = level*(3.5e-14*exolib.planck(wl, 5500*u.K) +
                      exolib.planck(wl, 270*u.K) * 3.58e-8)
    sed = Sed(wl, spectrum )
    transmission = Sed(wl, np.ones(wl.size)*u.dimensionless_unscaled)
    
    
    zodi = [sed, transmission]

    return zodi
    
# def emission_estimation(N, temp, emissivity, transmission):
 
#     t_i = transmission.sed**(1./N) # estimate transmission of each optical surface 
#     for i in range(0,N): #loop through each surface i=0 to i = N-1
#         if i == 0 :
#              emission_sed =  emissivity*exolib.planck(transmission.wl, temp)*t_i**(N-1-i)
#         else:
#             emission_sed = emission_sed + emissivity*exolib.planck(transmission.wl, temp)*t_i**(N-1-i)
#         idx = np.argwhere(transmission.sed==0) 
#         emission_sed[idx] = 0.0* emission_sed.unit            
            
#     emission=  Sed(transmission.wl, emission_sed) 
      
#     return emission


# def optical_emission(opt):
#     #telescope     
#       N = np.int(opt.common_optics.emissions.optical_surface.no_surfaces)
#       temp = opt.common_optics.emissions.optical_surface.val
#       emissivity = np.float(opt.common_optics.emissions.optical_surface.emissivity) 
      
#       opt.telescope_emission = emission_estimation(N, temp, emissivity, opt.telescope_transmission)
#       exolib.sed_propagation(opt.telescope_emission,opt.channel_transmission)
                   
#     #channel
#       N = np.int(opt.channel.emissions.optical_surface.no_surfaces)
#       temp = opt.channel.emissions.optical_surface.val
#       emissivity = np.float(opt.channel.emissions.optical_surface.emissivity)  
#       opt.channel_emission = emission_estimation(N, temp, emissivity, opt.channel_transmission)
      
#       emission = Sed(opt.channel_transmission.wl, opt.telescope_emission.sed+ opt.channel_emission.sed )

#       return emission  
  
    

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
          plt.figure('Total transmission after each optical surface added, OTE + channel %s'%(opt.channel.name))
          plt.plot(opt.x_wav_osr,TR)   
          plt.ylim(0,1.2)
        
          plt.figure('Transmission of each optical element added in the OTE %s'%(opt.channel.name))
          plt.plot(tr.wl, tr.sed) 
          plt.ylim(0,1.2)     
    
      em = Sed(dtmp[:,0]*u.um,dtmp[:,2]*u.dimensionless_unscaled)
      em.rebin(opt.x_wav_osr)

      exolib.sed_propagation(instrument_emission, tr, emissivity=em, temperature=op())
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

      exolib.sed_propagation(instrument_emission, tr, emissivity=em,temperature=op())
                        
            
  emission=  instrument_emission
      
  return emission
    
    

def run(opt):
    
      exosim_msg('Running backgrounds module ...\n ', opt.diagnostics)
      
      opt.zodi, zodi_transmission  = zodical_light(opt)
      

      # opt.emission = optical_emission(opt)
      opt.emission = emission_from_file(opt)
      

           
      exolib.sed_propagation(opt.star.sed, zodi_transmission)
      
    
 

      ch = opt.channel

      Omega_pix = 2.0*np.pi*(1.0-np.cos(np.arctan(0.5/ch.wfno())))*u.sr
      Apix = ((ch.detector_pixel.pixel_size.val).to(u.m))**2 

      opt.zodi.sed *= Apix*Omega_pix*opt.Re *u.electron/u.W/u.s 
      opt.emission.sed *= Apix*Omega_pix*opt.Re *u.electron/u.W/u.s 

      exolib.sed_propagation(opt.zodi, opt.total_transmission)
      
      opt.zodi.rebin(opt.x_wav_osr)
      opt.emission.rebin(opt.x_wav_osr)   

      opt.zodi.sed     *= opt.d_x_wav_osr 
      opt.emission.sed *= opt.d_x_wav_osr
      
      zodi_photons_sed = copy.deepcopy(opt.zodi.sed)
      emission_photons_sed = copy.deepcopy(opt.emission.sed)
      
      if ch.slit_width.val == 0:
          slit_size =  opt.fpn[1]*2  # to ensure all wavelengths convolved onto all pixels in slitless case
      else:
          slit_size =  ch.slit_width.val
      
      # need to work on this to make more accurate : filter wavelength range might not correspond to wl sol on the detector, i.e. some wavelengths not being included that should be
      opt.zodi.sed     = scipy.signal.convolve(opt.zodi.sed, 
		      np.ones(np.int(slit_size*ch.osf())), 
		      'same') * opt.zodi.sed.unit
      opt.emission.sed = scipy.signal.convolve(opt.emission.sed.value, 
		      np.ones(np.int(slit_size*ch.osf())), 
		      'same') * opt.emission.sed.unit

      opt.zodi_sed_original = copy.deepcopy(opt.zodi.sed)
      opt.emission_sed_original = copy.deepcopy(opt.emission.sed)
      
      # exosim_plot('zodi spectrum', opt.diagnostics, xdata = opt.zodi.wl, ydata=opt.zodi.sed)
      # exosim_plot('emission spectrum', opt.diagnostics, xdata = opt.emission.wl, ydata=opt.emission.sed)
                

    
      return opt
      

    
  