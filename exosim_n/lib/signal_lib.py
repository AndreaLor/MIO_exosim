"""
exosim_n

Signal library

"""
import numpy as np
from exosim_n.lib.exosim_n_lib import exosim_n_msg, exosim_n_plot
from exosim_n.lib import exosim_n_lib, noise_lib 
from astropy import units as u

def fast_method(opt):
       
  idx0, idx1, idxA, idxB = noise_lib.crop_array(opt) 
  # apply x crop to 1D arrays    
  opt.x_wav_osr0 = opt.x_wav_osr*1
  opt.x_wav_osr = opt.x_wav_osr[idxA*3:idxB*3]
  opt.x_pix_osr = opt.x_pix_osr[idxA*3:idxB*3]
  opt.cr_wl = opt.cr_wl[idxA:idxB]
  opt.cr = opt.cr[idxA:idxB]
  if opt.timeline.apply_lc.val ==1:
       opt.ldc = opt.ldc[:,idxA:idxB]
       opt.lc  = opt.lc[idxA:idxB]   
  opt.zodi.sed = opt.zodi.sed[idxA*3:idxB*3]
  opt.emission.sed = opt.emission.sed[idxA*3:idxB*3]
  opt.quantum_yield.sed = opt.quantum_yield.sed[idxA*3:idxB*3]
  opt.qy_zodi  = opt.qy_zodi[idxA*3:idxB*3]
  opt.qy_emission = opt.qy_emission[idxA*3:idxB*3]
  opt.exp_sig = opt.exp_sig[idxA:idxB] #expected final star signal in R bin from detector module sanity check
  # apply x and y crop to 2D arrays     
  opt.qe = opt.qe[idx0:idx1][:, idxA:idxB]
  opt.qe_uncert = opt.qe_uncert[idx0:idx1][:, idxA:idxB]
  opt.fp = opt.fp[idx0*3:idx1*3][:, idxA*3:idxB*3]
  opt.fp_signal = opt.fp_signal[idx0*3:idx1*3][:, idxA*3:idxB*3]
  

  # obtain proxy bkg signal for pipeline background subtraction step   
  opt.bkg_signal = noise_lib.proxy_bkg(opt)
  exosim_n_msg ("fast method used, focal place image reduced from  %s x %s to %s x %s"%(opt.fp_original.shape[0]/3, opt.fp_original.shape[1]/3 , opt.fp.shape[0]/3, opt.fp.shape[1]/3), opt.diagnostics)

  return opt

def initiate_signal(opt):
  opt.signal_only = noise_lib.obtain_signal_only(opt) 
 
  return opt

def initiate_noise(opt):
    opt.combined_noise  = np.zeros_like(opt.signal)
    opt.star_signal  = opt.signal*1 
    return opt

def apply_lc(opt):
    if opt.timeline.apply_lc.val ==0:
      exosim_n_msg ("OMITTING LIGHT CURVE...", opt.diagnostics)
      return opt
    else:
      exosim_n_msg ("APPLYING LIGHT CURVE...", opt.diagnostics)
      print (opt.lc.max(), opt.lc.min())
      print (opt.signal.max(), opt.signal.min())
      signal = opt.signal*1
      lc = opt.lc*1
      signal = signal*lc
      opt.signal = signal 
      return opt
    
def apply_systematic(opt):
    if opt.simulation.sim_use_systematic_model.val == 0: 
      return opt
    else:
      exosim_n_msg ("Applying systematic model...", opt.diagnostics)
      signal = opt.signal*1
      syst = opt.syst_grid*1
      signal = signal*syst
      opt.signal = signal 
      return opt  
   
def apply_non_stellar_photons(opt):
  blank_fp_shape = np.ones((int(opt.fp.shape[0]/3), int(opt.fp.shape[1]/3), len(opt.frames_per_ndr)))
  opt.zodi_signal = np.zeros_like(blank_fp_shape)*u.electron
  opt.emission_signal = np.zeros_like(blank_fp_shape)*u.electron
  if opt.background.EnableZodi.val == 1:
      exosim_n_msg ("APPLYING ZODI..." , opt.diagnostics)     
      zodi_signal = blank_fp_shape* opt.zodi.sed[opt.offs::opt.osf][np.newaxis, :, np.newaxis]      
      opt.zodi_signal = zodi_signal * opt.frame_time * opt.frames_per_ndr
      opt.signal = opt.signal + opt.zodi_signal
  else:
      exosim_n_msg ("NOT APPLYING ZODI..." , opt.diagnostics)  
  if opt.background.EnableEmission.val == 1:
      exosim_n_msg ( "APPLYING EMISSION..."  ,  opt.diagnostics  )         
      emission_signal = blank_fp_shape* opt.emission.sed[opt.offs::opt.osf][np.newaxis, :, np.newaxis]       
      opt.emission_signal = emission_signal * opt.frame_time * opt.frames_per_ndr
      opt.signal = opt.signal + opt.emission_signal     
  else:
      exosim_n_msg ("NOT APPLYING EMISSION..." , opt.diagnostics) 
  return opt

def apply_prnu(opt):
    opt.qe_grid = opt.qe # used for flat field in pipeline 
    exosim_n_msg ("mean and standard deviation of flat field: should match applied flat in data reduction %s %s"%(opt.qe.mean(), opt.qe.std()),  opt.diagnostics  )
    exosim_n_msg ("standard deviation of flat field uncertainty %s"%(opt.qe_uncert.std()), opt.diagnostics  )     
    applied_qe= opt.qe*opt.qe_uncert 
    if opt.noise.ApplyPRNU.val == 1: 
        exosim_n_msg ("PRNU GRID BEING APPLIED", opt.diagnostics) 
        opt.signal= (opt.signal.transpose() * applied_qe.transpose() ).transpose()  
    else:
        exosim_n_msg ("PRNU GRID NOT APPLIED...", opt.diagnostics) 
    return opt

def apply_dc(opt):
  blank_fp_shape = np.ones((int(opt.fp.shape[0]/3), int(opt.fp.shape[1]/3), len(opt.frames_per_ndr)))
  opt.dc_signal = np.zeros_like(blank_fp_shape)*u.electron
  if opt.background.EnableDC.val ==1:
      exosim_n_msg ("DARK CURRENT being added...%s"%(opt.channel.detector_pixel.Idc.val ), opt.diagnostics  )
      opt.dc_signal = blank_fp_shape*  opt.channel.detector_pixel.Idc() * opt.frame_time * opt.frames_per_ndr                 
      opt.signal = opt.signal + opt.dc_signal     
  else:
      exosim_n_msg ("DARK CURRENT.. not... being added...", opt.diagnostics)
  return opt
     
# def apply_quantum_yield(opt):     
#    # Apply weighted quantum yield to signal
#    exosim_n_msg ("Applying quantum yield", opt.diagnostics)
#    qy_emission = opt.qy_emission[1::3][np.newaxis, :, np.newaxis]
#    qy_zodi = opt.qy_zodi[1::3][np.newaxis, :, np.newaxis]    
#    qy_star = opt.quantum_yield.sed[1::3][np.newaxis, :, np.newaxis]
#    signal_total = opt.star_signal + opt.zodi_signal + opt.emission_signal + opt.dc_signal
#    signal_total = np.where(signal_total<=0, 1e-10*u.electron, signal_total)
#    opt.quantum_yield_total = (qy_star*opt.star_signal + qy_zodi*opt.zodi_signal + qy_emission*opt.emission_signal + 1*opt.dc_signal) / signal_total
#    opt.quantum_yield_total = np.where(opt.quantum_yield_total<1, 1, opt.quantum_yield_total)
#    exosim_n_msg(f'max, min quantum yield applied: {opt.quantum_yield_total.max()}, {opt.quantum_yield_total.min()}', opt.diagnostics)
#    opt.signal = opt.signal*opt.quantum_yield_total
#    opt.combined_noise = opt.combined_noise*opt.quantum_yield_total

#    return opt

def apply_combined_noise(opt):
    opt.signal = opt.signal + opt.combined_noise
    return opt

def apply_utr_correction(opt):
  if opt.noise.EnableShotNoise.val == 1:
       if opt.simulation.sim_full_ramps.val == 0 and opt.simulation.sim_use_UTR_noise_correction.val == 1:
             exosim_n_msg ("applying correction to photon noise for UTR read", opt.diagnostics  )
             n = opt.projected_multiaccum
             opt.combined_noise, scale = noise_lib.poission_noise_UTR_correction(n, opt.combined_noise.value)              
             opt.combined_noise = opt.combined_noise * u.electron
  return opt

def make_ramps(opt):
  exosim_n_msg ("check point 6.1 %s"%(opt.signal.max()) , opt.diagnostics  )
  # Make ramps 
  for i in range(0, opt.n_ndr, opt.effective_multiaccum):
      opt.signal[...,i:i+opt.effective_multiaccum] = np.cumsum(opt.signal[...,i:i+opt.effective_multiaccum], axis=2) 
  return opt

# def apply_ipc(opt):
#   # Apply IPC: should be before read noise
#   if opt.simulation.sim_use_ipc.val == 1 and opt.channel.instrument.val !='MIRI': 
#        exosim_n_msg('Applying IPC...', opt.diagnostics)    
#        ipc_kernel = np.load(opt.channel.detector_array.ipc_kernel.val.replace('__path__', '%s/%s'%(opt.exosim_n_path, 'exosim_n')))
#        opt.signal =  noise_lib.apply_ipc(opt.signal, ipc_kernel)
#   return opt

def apply_read_noise(opt): 
  if opt.noise.EnableReadoutNoise.val == 1:
      exosim_n_msg ("READ NOISE... being added...", opt.diagnostics)
      opt.signal = noise_lib.apply_read_noise(opt.signal,opt)         
  else:
      exosim_n_msg ("READ NOISE...not... being added..." , opt.diagnostics  )                
  exosim_n_msg ("check point 7 %s %s"%(opt.signal.max(), opt.signal.min()) , opt.diagnostics  )
  return opt
  
def apply_jitter(opt):   
   if opt.noise.EnableSpatialJitter.val  ==1 or opt.noise.EnableSpectralJitter.val==1:
       opt = noise_lib.simulate_jitter(opt)
   else:
       opt = noise_lib.jitterless(opt)
   if opt.noise.EnableSpatialJitter.val== 0:
        opt.pointing_timeline[:,2] =  opt.pointing_timeline[:,2]*0
   if opt.noise.EnableSpectralJitter.val== 0:
        opt.pointing_timeline[:,1] =  opt.pointing_timeline[:,1]*0   
   opt.signal[np.isnan(opt.signal)] = 0    
   return opt 
       
      

def apply_poisson_noise(opt):    
    if opt.noise.EnableShotNoise.val == 1:    
       exosim_n_msg ("Applying Poisson noise", opt.diagnostics)
       opt.signal = np.where(opt.signal >= 0.0, opt.signal, 0)
       opt.combined_noise = noise_lib.calc_poission_noise(opt.signal.value) *u.electron
       
    else:
       exosim_n_msg ("Poisson noise not being applied...", opt.diagnostics)
    return opt

 