"""
exosim_n

Instrument library

"""
import numpy as np
from exosim_n.lib import exosim_n_lib
from exosim_n.lib.exosim_n_lib import exosim_n_msg, exosim_n_plot, planck
from astropy import units as u
from astropy import constants as const
from astropy.io import fits
from exosim_n.classes.sed import Sed
import scipy
import matplotlib.pyplot as plt
import sys, os
from scipy import interpolate
import copy
import pickle

def useInterp(opt):    
      offset = opt.fpn[1]*opt.channel.detector_pixel.pixel_size()/2.0
 
      dtmp=np.loadtxt(opt.channel.dispersion.path.replace(
	  '__path__', opt.__path__), delimiter=',')
      
      if opt.channel.is_spec.val == 1:
          try:
              opt.channel.disp_adjust.val
          except AttributeError:
              opt.channel.disp_adjust.val = 0*u.um

      disp_adjust = opt.channel.disp_adjust.val
      print (disp_adjust)
   
      ld = scipy.interpolate.interp1d(dtmp[...,2]*u.um + offset.to(u.um) +disp_adjust, dtmp[...,0], bounds_error=False,kind='slinear',
					fill_value=0.0)          
      x_pix_osr = np.arange(opt.fp.shape[1]) * opt.fp_delta +  opt.fp_delta/2.0
      x_wav_osr = ld(x_pix_osr.to(u.um))*u.um # walength on each x pixel 
 

      pos_y = dtmp[...,1]     
      x_edge = np.arange(opt.fpn[1]+1) * opt.fp_delta*3
      opt.x_wav_edges = ld(x_edge.to(u.um))*u.um
      if pos_y[0] !=0:
          y_pos_osr = scipy.interpolate.interp1d(dtmp[...,0],  dtmp[...,1], bounds_error=False, fill_value=0.0)(x_wav_osr)
          # there must be a fill value as each x_wav_osr must have a corresponding y_pos_osr
          y_pos_osr = (y_pos_osr/((opt.fp_delta).to(u.um))).value
      else:
          y_pos_osr = []
      if opt.diagnostics == 1:
          import matplotlib.pyplot as plt
          plt.figure('wav sol test 1')
          plt.plot(x_wav_osr, 'ro')
          plt.figure('wav sol test 2')
          plt.plot(x_wav_osr, np.gradient(x_wav_osr), 'ro')
      # smooth any irregularities in the wavelength solution
      box = 7
      bbox = np.ones(box)/box
      idx = np.argwhere(x_wav_osr.value>0).T[0]
      idx0 = np.argwhere(x_wav_osr.value==0).T[0] 
      aa = x_wav_osr[idx]
      bb0 = aa[0:box]
      bb1 = aa[-box:]
      aa = np.convolve(aa,bbox,'same')
      aa[0:box] = bb0; aa[-box:]= bb1
      x_wav_osr[idx] = aa

      if opt.diagnostics == 1:
          plt.figure('wav sol test 1')
          plt.plot(x_wav_osr, 'bx')
          plt.figure('wav sol test 2')
          plt.plot(x_wav_osr, np.gradient(x_wav_osr), 'bo')    
    
      return  x_wav_osr, x_pix_osr.to(u.um), y_pos_osr
      

def usePoly(opt): #this needs updating   
      offset = opt.fpn[1]*opt.channel.detector_pixel.pixel_size()/2.0    
  
      dtmp=np.loadtxt(opt.channel.dispersion.path.replace(
	  '__path__', opt.__path__), delimiter=',')      
      pos = (dtmp[...,2]*u.um  + offset.to(u.um)) + (opt.channel.dispersion() ).to(u.um)
      wav = dtmp[...,0]
      pos_y = dtmp[...,1] 
      r = 7 # degree 7 works well
      z = np.polyfit(pos, wav, r)
      if pos_y[0] !=0:
          r_y = 7 # degree 7 works well
          z_y = np.polyfit(wav, pos_y, r_y)
      x_pix_osr = (np.arange(opt.fp.shape[1])*opt.fp_delta + opt.fp_delta/2).to(u.um)
      x_wav_osr = 0
      for i in range (0,r+1):
          x_wav_osr = x_wav_osr + z[i]*x_pix_osr.value**(r-i) 
      x_wav_osr =x_wav_osr*u.um           
      idx = np.argwhere(x_wav_osr<wav.min())
      x_wav_osr[idx] = 0*u.um 
      idx = np.argwhere(x_wav_osr>wav.max())
      x_wav_osr[idx] = 0*u.um               
      idx = np.argwhere(x_pix_osr.value < pos.min())
      x_wav_osr[idx] = 0*u.um
      idx = np.argwhere(x_pix_osr.value > pos.max())
      x_wav_osr[idx] = 0*u.um                   
      if pos_y[0] !=0:
          y_pos_osr =0
          for i in range (0,r_y+1):
              y_pos_osr = y_pos_osr + z_y[i]*x_wav_osr.value**(r_y-i)        
          y_pos_osr = (y_pos_osr/((opt.fp_delta).to(u.um))).value
      else:
          y_pos_osr = []         
      return x_wav_osr, x_pix_osr, y_pos_osr
            

      
def sanity_check(opt):
    wl = opt.x_wav_osr[1::3][::-1]
    wl0= opt.x_wav_osr[1::3]
    del_wl = abs(np.gradient(wl))
#    del_wl = opt.d_x_wav_osr[1::3]*3
    star_spec = opt.star_sed
    star_spec.rebin(wl)
    T = opt.planet.planet.star.T
    trans_sed = opt.total_transmission.sed*u.dimensionless_unscaled
    trans = Sed(opt.total_transmission.wl[::-1],trans_sed[::-1]) 
    trans.rebin(wl)
    QE = opt.qe_spec
    QE.sed = QE.sed[::-1]
    QE.wl = QE.wl[::-1]
    QE.rebin(wl)
    quantum_yield = 1
    Rs = (opt.planet.planet.star.R).to(u.m)
    D = (opt.planet.planet.star.d).to(u.m)
    n= quantum_yield* trans.sed*del_wl*np.pi*planck(wl,T)*(Rs/D)**2*opt.Aeff*QE.sed/(const.h.value*const.c.value/(wl*1e-6))
      # opt.Re = qe.sed * (qe.wl).to(u.m)/(const.c.value * const.h.value * u.m)
    n2= quantum_yield* trans.sed*del_wl*star_spec.sed*opt.Aeff*QE.sed/(const.h*const.c/(wl*1e-6))
    
    sim_sig = quantum_yield*opt.fp_signal[1::3,1::3].sum(axis=0)
    R = opt.pipeline.pipeline_R.val
    del_wav = wl/R
    opt.exp_sig  = opt.t_int*del_wav[::-1]*sim_sig/del_wl[::-1]
    if opt.diagnostics ==1:
        plt.figure('sanity check 1 - check focal plane signal')
        plt.plot(wl,n, 'b^', label='BB check')
        plt.plot(wl,n2, 'r+', label='Phoenix check')  # not convolved with PSF unlike exosim_n, so peak may be higher
        plt.plot(wl0, sim_sig, 'gx', label='exosim_n')
        plt.ylabel('e/s/pixel col'); plt.xlabel('pixel col wavelength (microns)')
        plt.legend(loc='best')
     ################      
        plt.figure('sanity check 2 - expected final star signal in R bin of %s'%((R)))
        plt.plot(wl, opt.exp_sig)
        plt.ylabel('e/bin'); plt.xlabel('Wavelength (microns)')
      ################      
        plt.figure('sanity check 3 - expected photon noise (sd) in R bin of %s'%((R)))
        plt.plot(wl, opt.exp_sig**0.5)    
        plt.ylabel('e/bin'); plt.xlabel('Wavelength (microns)')     


def get_psf(opt):  # there is an issue with PSF interpolation of very abberent PSFs due to aliasing - have to look at further
      if opt.simulation.sim_use_wfe.val == 0: 
            psf = exosim_n_lib.Psf(opt.x_wav_osr.value, opt.channel.wfno.val, opt.channel.wfno_y.val, opt.fp_delta.value, shape='airy')  
            psf[np.isnan(psf)] =0
            opt.psf_type = 'airy'
      elif opt.simulation.sim_use_wfe.val  == 1:
            # zfile = '%s/systems/%s/PSF/PSF_stack_%s.pickle'%(opt.__path__,opt.system,opt.channel.name)
            zfile = '%s/systems/%s/PSF/%s_psf_stack.fits'%(opt.__path__,opt.system,opt.channel.name)

            if opt.channel.is_spec.val == False:
                psf = exosim_n_lib.Psf_photometer(zfile, opt.fp_delta, opt.channel.osf(), opt.fpn, opt.x_wav_osr)   
            else:
                # psf =  Psf_spectrometer_file(zfile, opt)
                psf =  exosim_n_lib.Psf_spectrometer(zfile, opt.fp_delta, opt.channel.osf(), opt.fpn, opt.x_wav_osr)   
          
      opt.psf = psf
      exosim_n_msg("PSF shape %s, %s"%(opt.psf.shape[0],opt.psf.shape[1] ), opt.diagnostics)     
      exosim_n_plot('psf check', opt.diagnostics, image=True, image_data = psf[..., int(psf.shape[2]/2)])
    
      sum1=[]
      for i in range(psf.shape[2]):
          # print (psf[...,i].sum(), opt.x_wav_osr[i])
          sum1.append(psf[...,i].sum())
          if psf[...,i].sum()!=0:
               # psf[...,i] =psf[...,i]/psf[...,i].sum() 
                 if np.round(psf[...,i].sum(),3) !=1.0:    
                     exosim_n_msg('error... check PSF normalisation %s %s'%(psf[...,i].sum(), opt.x_wav_osr[i]), 1)
              
      exosim_n_plot('test7 - psf sum vs subpixel position (should be 1)', opt.diagnostics, 
                   xdata=opt.x_pix_osr, ydata = sum1, marker='bo')
    
      return psf      
  
def get_focal_plane(opt):
      if opt.channel.is_spec.val == True:
          # Populate focal plane with monochromatic PSFs
            buffer = int(opt.psf.shape[1]/2)   
            fp = np.zeros((opt.fpn[0]*3 , opt.fpn[1]*3 + buffer*2)    )       
            i0 = np.array([fp.shape[0]/2 - opt.psf.shape[0]/2]* fp.shape[1]).astype(np.int)                 
            i1 = i0 + opt.psf.shape[0]  
            for k in range (0,  fp.shape[1]-buffer*2,1):
                fp[i0[k]: i1[k] ,k:k+opt.psf.shape[1] ] += opt.psf[...,k]*opt.star.sed.sed[k].value   
            opt.fp_signal = fp[:,buffer:-buffer]  
            if opt.background.EnableSource.val == 1 or opt.background.EnableAll.val == 1:    
                 opt.fp = opt.fp_signal*1
            # exosim_n_plot('test3a', opt.diagnostics, xdata=opt.x_wav_osr, ydata = opt.fp.sum(axis=0), marker='bo')

      else: # for photometer
          j0 = np.repeat(int(opt.fp.shape[1]/2 - opt.psf.shape[1]/2), opt.x_wav_osr.size)
          j1 = j0 + opt.psf.shape[1]
          idx = np.where((j0>=0) & (j1 < opt.fp.shape[1]))[0]
          i0 = int(opt.fp.shape[0]/2 - opt.psf.shape[0]/2 + opt.offs)
          i0+=1         
          i1 = i0 + opt.psf.shape[0]  
          
          if opt.background.EnableSource.val == 1 or opt.background.EnableAll.val == 1:     
              for k in idx: # actual signal 
                  opt.fp[i0:i1, j0[k]:j1[k]] += opt.psf[...,k] * opt.star.sed.sed[k].value              
          for k in idx: # used for getting sat time, and sizing
              opt.fp_signal[i0:i1, j0[k]:j1[k]] += opt.psf[...,k] * opt.star.sed.sed[k].value  
      import matplotlib.pyplot as plt

      return opt.fp, opt.fp_signal
  



      
def get_planet_spectrum(opt):
    if opt.channel.is_spec.val == True:
        exosim_n_plot('planet sed 1', opt.diagnostics, 
                     xdata=opt.planet.sed.wl, ydata = opt.planet.sed.sed, 
                       ylim=[0,1])
        opt.planet.sed = opt.planet.sed
    else:
        pce = opt.PCE
        cr = opt.planet.sed.sed
        idx = np.argwhere(cr!=0).T[0]
        cr = cr[idx]
        pce= pce[idx]
        planet_cr  = np.sum((pce*cr)) / pce.sum()  # weighted average of the cr
        planet_cr = np.repeat(planet_cr, len(opt.x_wav_osr))
        planet_wl = np.repeat(opt.phot_wav, len(opt.x_wav_osr))
        opt.planet.sed = Sed(planet_wl, planet_cr*u.dimensionless_unscaled) 
         
 
    return opt.planet.sed
# alternate code that does not rely on multiplying by the star.sed first- not sure if this works correctly yet.    
#          planet_fp = copy.deepcopy(opt.fp)*0
#          i0 = np.array([opt.fp.shape[0]/2 - psf.shape[0]/2 + opt.offs]*len(j1))
#          i1 = i0 + psf.shape[0]
#          
#          for k in idx: 
#              planet_fp[i0[k]:i1[k], j0[k]:j1[k]] += psf[...,k] * opt.planet_sed_original[k] 
#
#  
#          planet_response = planet_fp.sum(axis=0)
#          plt.figure('planet sed 111')
#          plt.plot(opt.planet.sed.wl, planet_response)
#          
#          plt.figure('test')
#          plt.imshow(planet_fp)

def convolve_prf(opt):

      kernel, kernel_delta = exosim_n_lib.PixelResponseFunction(opt, 
        opt.psf.shape[0:2],
        7*opt.channel.osf(),   
        opt.channel.detector_pixel.pixel_size(),
        lx = opt.channel.detector_pixel.pixel_diffusion_length())
      exosim_n_msg ("kernel sum %s"%(kernel.sum()), opt.diagnostics) 
      exosim_n_msg ("check 3.9 - unconvolved FP max %s"%(opt.fp.max()) , opt.diagnostics)      
      exosim_n_plot('test3', opt.diagnostics, xdata=opt.x_wav_osr, ydata = opt.fp.sum(axis=0), marker='bo')
      exosim_n_plot('test4', opt.diagnostics, xdata=opt.x_pix_osr, ydata = opt.x_wav_osr, marker='bo')                  
      exosim_n_plot('test5', opt.diagnostics, xdata=opt.x_pix_osr, ydata = opt.fp.sum(axis=0), marker='bo')
      exosim_n_plot('test6', opt.diagnostics, xdata=opt.x_wav_osr, ydata = opt.star.sed.sed, marker='bo') 
      opt.fp = exosim_n_lib.fast_convolution(opt.fp, opt.fp_delta, kernel, kernel_delta)
      opt.fp_signal = exosim_n_lib.fast_convolution(opt.fp_signal, opt.fp_delta, kernel, kernel_delta)
      exosim_n_msg ("check 4 - convolved FP max %s"%(opt.fp.max()), opt.diagnostics)    
      opt.kernel = kernel
      opt.kernel_delta = kernel_delta     
      # Fix units
      opt.fp = opt.fp*opt.star.sed.sed.unit 
      opt.fp_signal = opt.fp_signal*opt.star.sed.sed.unit  
      return opt.fp, opt.fp_signal
  
def exposure_timing(opt):    
    
    ## Find count rate with diffuse radiation
      fp_copy = copy.deepcopy(opt.fp_signal)
      fp_count_no_bkg = 1*fp_copy[1::3,1::3] 
      fp_count = opt.zodi.sed[1::3] + opt.emission.sed[1::3]+ fp_copy[1::3,1::3] 
   
      fp_count += opt.channel.detector_pixel.Idc.val 
      
      exosim_n_msg ("check 4a - %s"%(fp_count.max()), opt.diagnostics)

      fp_count = fp_count.value
      fp_count_no_bkg = fp_count_no_bkg.value 

      exosim_n_msg ("check 5 - %s"%(fp_count.max()), opt.diagnostics)
      FW = opt.channel.detector_pixel.full_well.val       
      A,B = np.unravel_index(fp_count.argmax(), fp_count.shape)
      exosim_n_msg ("maximum index and count with all backgrounds %s %s %s"%(A,B, fp_count.max()), opt.diagnostics)      
      A,B = np.unravel_index(fp_count_no_bkg.argmax(), fp_count_no_bkg.shape)
      exosim_n_msg ("maximum index and count with no backgrounds %s %s %s"%(A,B, fp_count_no_bkg.max()), opt.diagnostics)   
      exosim_n_msg ("full well in electron %s"%(FW), opt.diagnostics)   
      exosim_n_msg ("saturation time assuming 100 percent full well with all backgrounds %s"%((FW / fp_count.max())), opt.diagnostics)     
      exosim_n_msg ("full well percentage chosen for saturation limit: %s"%(opt.observation.obs_fw_percent.val), opt.diagnostics)
      opt.sat_time = ((FW / fp_count.max()) *opt.observation.obs_fw_percent.val/100.0 )*u.s     
      opt.sat_time_no_bkg = ((FW / fp_count_no_bkg.max()) *opt.observation.obs_fw_percent.val/100.0 )*u.s
      opt.sat_limit = u.electron*FW*opt.observation.obs_fw_percent.val/100.0  
      opt.sat_limit_fw = u.electron*FW 
      exosim_n_msg ('saturation limit %s'%(opt.sat_limit), opt.diagnostics)   
      exosim_n_msg ("saturation time with all backgrounds %s"%(opt.sat_time), opt.diagnostics)     
      exosim_n_msg ("saturation time with no backgrounds or dc %s"%(opt.sat_time_no_bkg), opt.diagnostics)
      
      opt.t_g = opt.t_f = opt.channel.array_read_time.val
      opt.dead_time = (opt.timeline.nGND.val+ opt.timeline.nRST.val)* opt.t_g
      opt.zero_time = opt.timeline.nNDR0.val* opt.t_g 
      
      # n_full_ndrs are the number of ndrs not including the zeroth ndr
      if opt.observation.obs_use_sat.val == 1: 
          exosim_n_msg('Using saturation time to set n_groups', opt.diagnostics)
          n_full_ndrs = int((opt.sat_time-opt.zero_time)/opt.t_g) # does not include reset group (assume this is after final read so saturation in this period does not affect read counts)
          if n_full_ndrs <1:
              n_full_ndrs=1
      else:
          exosim_n_msg('Using user-defined n_groups', opt.diagnostics)
          n_full_ndrs = opt.observation.obs_n_groups.val -1        
           
      exosim_n_msg('t_f %s'%(opt.t_f), opt.diagnostics)
      exosim_n_msg('t_g %s'%(opt.t_g), opt.diagnostics)
      exosim_n_msg('dead time %s'%(opt.dead_time), opt.diagnostics)
      exosim_n_msg('zero time %s'%(opt.zero_time), opt.diagnostics)
      exosim_n_msg('n_groups %s'%(n_full_ndrs + 1), opt.diagnostics)              
      # t_sim is not currently used, and by default is set to the same value as t_f
      opt.t_sim = opt.t_f             
      opt.t_int =  (n_full_ndrs)*opt.t_g # this is the 'cds time': i.e. misses the zeroth ndr time
      opt.t_cycle = opt.t_int + opt.dead_time + opt.zero_time  
          
      if opt.t_cycle > opt.sat_time:
          exosim_n_msg ("\n****** Warning!!!!!!!!!  : some pixels will exceed saturation limit ******\n", opt.diagnostics  )
          opt.sat_flag = 1
          # sys.exit()      
      else:
          exosim_n_msg ("\n****** OK!!!!!!  Cycle time within saturation time ******\n", opt.diagnostics  )
          opt.sat_flag = 0
    ######  Set effective multiaccum
      if opt.simulation.sim_full_ramps.val == 0:
          exosim_n_msg ("Approximating ramps with corrected CDS method, so only 2 NDRs simulated", 1)
          opt.effective_multiaccum = 2 # effective multiaccum is what is implemented in sim
          opt.projected_multiaccum = int(n_full_ndrs+1)
      else:
          opt.effective_multiaccum = int(n_full_ndrs+1)
          opt.projected_multiaccum = int(n_full_ndrs+1)       
      exosim_n_msg ("projected multiaccum: %s"%(opt.projected_multiaccum), opt.diagnostics)
      exosim_n_msg ("effective multiaccum: %s"%(opt.effective_multiaccum), opt.diagnostics)                                
      opt.exposure_time = (opt.t_int + opt.dead_time + opt.zero_time) #same as t_cycle
      return opt
    
    
def getPhotometerWavelengths(opt):
    
      wav_lo = opt.channel.pipeline_params.wavrange_lo.val
      wav_hi = opt.channel.pipeline_params.wavrange_hi.val
      idx = np.argwhere((opt.star.sed.wl.value>=wav_lo)&(opt.star.sed.wl.value<=wav_hi)).T[0]
      wav = opt.star.sed.wl[idx]
      tr_ =np.array([1.]*len(wav))*u.dimensionless_unscaled
      opt.common_optics.optical_surface = opt.common_optics.optical_surface \
          if isinstance(opt.common_optics.optical_surface, list) \
          else [opt.common_optics.optical_surface]    
      for op in opt.common_optics.optical_surface:
          dtmp=np.loadtxt(op.transmission.replace('__path__', opt.__path__), delimiter=',')
          tr = Sed(dtmp[:,0]*u.um,dtmp[:,1]*u.dimensionless_unscaled)
          tr.rebin(wav)  
          tr_ *=tr.sed
      opt.channel.optical_surface = opt.channel.optical_surface if isinstance(opt.channel.optical_surface, list) else \
                                                             [opt.channel.optical_surface]                                                              
      for op in opt.channel.optical_surface:
          dtmp=np.loadtxt(op.transmission.replace('__path__', opt.__path__), delimiter=',')
          tr = Sed(dtmp[:,0]*u.um, dtmp[:,1]*u.dimensionless_unscaled)
          tr.rebin(wav)
          tr_ *=tr.sed
      plt.figure('photometer wavelength range')
      plt.plot(wav, tr_)
      dtmp=np.loadtxt(opt.channel.qe().replace('__path__', opt.__path__), delimiter=',')
      qe = Sed(dtmp[:,0]*u.um, dtmp[:,1]*u.dimensionless_unscaled)
      qe.rebin(wav)
      tr_ *=qe.sed      
      phot_wav  = (wav*tr_).sum() / tr_.sum()  # weighted average of the wavelength
      print ('photometer central wavelength:', phot_wav)
 
      plt.figure('photometer wavelength range')
      plt.plot(wav, tr_)

      x_wav_osr = np.linspace(wav[0].value, wav[-1].value, opt.fp.shape[1])*u.um
      # x_wav_osr = wav
      
      x_pix_osr = np.arange(len(x_wav_osr)) # irrelevant here     
      y_pos_osr = []
      return x_wav_osr, x_pix_osr,y_pos_osr, phot_wav 
       
