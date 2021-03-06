"""
exosim_n 
 
Timeline module
 
"""
import numpy as np
from astropy import units as u
from exosim_n.lib.exosim_n_lib import exosim_n_msg

def run(opt):

  planet = opt.planet
  opt.T14 =   planet.T14
  
  opt.time_at_transit      = opt.T14*(0.5+  opt.observation.obs_frac_t14_pre_transit.val )
   
  opt.frame_time = opt.t_f  
  opt.T14 = (opt.T14).to(u.s)

  opt.multiaccum     = opt.effective_multiaccum    # Number of NDRs per exposure
  opt.allocated_time = (opt.timeline.nGND()+
		      opt.timeline.nNDR0()+
		      opt.timeline.nRST()) * opt.frame_time
    
  exosim_n_msg ('exposure_time %s'%(opt.exposure_time), opt.diagnostics)
  exosim_n_msg ('allocated_time %s'%(opt.allocated_time) ,  opt.diagnostics)
  exosim_n_msg ('effective multiaccum %s'%(opt.effective_multiaccum),  opt.diagnostics)
  NDR_time       = (opt.exposure_time-opt.allocated_time)/(opt.multiaccum-1)
        
     
  exosim_n_msg ("initial NDR time %s"%(NDR_time ),  opt.diagnostics)
   
  exosim_n_msg ("overheads %s %s %s"%(opt.timeline.nGND()*opt.frame_time, opt.timeline.nNDR0()*opt.frame_time, opt.timeline.nRST()*opt.frame_time),  opt.diagnostics)
    
  # find number of frames in each non-zeroth NDR  
  nNDR = np.round(NDR_time/opt.frame_time).astype(np.int).take(0)
    
  base = [opt.timeline.nGND.val, opt.timeline.nNDR0.val]
  for x in range(opt.multiaccum-1): base.append(nNDR)
  base.append(opt.timeline.nRST.val)  
    
 # Calcaulate exposure time (=total cycle time) and estimates how many exposures are needed
  opt.exposure_time = sum(base)*opt.frame_time
  opt.frames_per_exposure = sum(base)
  
  if opt.timeline.apply_lc.val ==1:
      exosim_n_msg ("Since light curve is implemented, observing time is set to 2x T14",  opt.diagnostics)
      opt.timeline.use_T14.val = 1
  
  if opt.timeline.use_T14.val ==1:
      total_observing_time = opt.T14*(1.0+opt.observation.obs_frac_t14_pre_transit.val+opt.observation.obs_frac_t14_post_transit.val)
      number_of_exposures = np.ceil((total_observing_time.to(u.s)/opt.exposure_time.to(u.s))).astype(np.int)
      number_of_exposures = number_of_exposures.value
  elif opt.timeline.use_T14.val ==0:   
      number_of_exposures = int(opt.timeline.n_exp.val)    
      if opt.timeline.obs_time.val >0:
            number_of_exposures = int(opt.timeline.obs_time.val.to(u.s) / opt.exposure_time) 
  opt.n_exp = number_of_exposures
 
  opt.total_observing_time = opt.exposure_time*opt.n_exp
  
  exosim_n_msg ("number of integrations %s"%(number_of_exposures),  opt.diagnostics)
  exosim_n_msg ("number of NDRs %s"%(number_of_exposures*opt.multiaccum),  opt.diagnostics)
  exosim_n_msg ("total observing time (hrs) %s"%((number_of_exposures*opt.exposure_time/3600).value),  opt.diagnostics)
  exosim_n_msg ("T14 %s"%(opt.T14),  opt.diagnostics)
  
  opt.frame_sequence=np.tile(base, number_of_exposures)  
  opt.time_sequence = opt.frame_time * opt.frame_sequence.cumsum()  
  
  # End time of each NDR
  opt.ndr_end_time = np.dstack([opt.time_sequence[1+i::len(base)] \
      for i in range(opt.multiaccum)]).flatten()
                  
  exosim_n_msg ("actual NDR time %s"%(opt.ndr_end_time[1]-opt.ndr_end_time[0]  ),  opt.diagnostics) 
  exosim_n_msg ("final exposure time %s"%(opt.exposure_time ) , opt.diagnostics) 
      
    # Number of frames contributing to each NDR
  opt.frames_per_ndr = np.dstack([opt.frame_sequence[1+i::len(base)] \
      for i in range(opt.multiaccum)]).flatten()
           
  opt.ndr_end_frame_number = np.round(opt.ndr_end_time/opt.frame_time).astype(int).value
  opt.duration_per_ndr  = opt.frames_per_ndr*opt.frame_time 
  opt.n_ndr = number_of_exposures*opt.multiaccum
  opt.ndr_list = np.arange(0,opt.n_ndr,1)
           
  return opt

