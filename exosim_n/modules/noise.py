"""
ExoSim_N
Noise module

"""

import numpy as np
from exosim_n.lib.exolib import exosim_msg, exosim_plot
from exosim_n.lib import exolib 
from astropy import units as u
import copy
import time
import sys
import scipy
from numba import jit, prange
     

def poission_noise (signal):    
    noise = np.random.poisson(signal.value) - signal.value
    return noise

def poission_noise_UTR_correction (n, noise):    
    alpha = 6.0*(n**2+1)/(5.0*n*(n+1))
    noise = noise*np.sqrt(alpha) #scale the noise
    return noise, alpha   

def read_noise (signal, sigma_ro):   
 
    signal = signal + np.random.normal(scale = sigma_ro.value, size = signal.shape)  
    return signal

def read_noise_UTR_correction (n, noise, sigma_ro):   
    # correction for utr noise: achieves same noise ar utr but using cds      
    sigma_ro_prime = np.sqrt(12.0*(n-1)/(n*(n+1)))*sigma_ro /np.sqrt(2)  
    noise = noise + np.random.normal(scale = sigma_ro_prime, size = noise.shape)
    return noise 

def combine_noise(noise1, noise2):
    sign = np.sign(noise1 + noise2)
    combined_noise = sign*(abs(noise1**2 + noise2**2))**0.5
    return combined_noise
    

def apply_ipc (data, ipc_kernel):
    data = data.value
    for i in range(data.shape[2]):
        data[...,i] = scipy.signal.fftconvolve(data[...,i], ipc_kernel, 'same')
    data = data*u.electron
    return data  
 
@jit(nopython=True) 
def create_jitter_noise(fp, osf, frames_per_ndr, frame_osf, jitter_x, jitter_y):		
    n_ndr = len(frames_per_ndr)
    fp_whole = fp[int(osf/2.)::osf, int(osf/2.)::osf]
    ndr_count = np.zeros((fp_whole.shape[0], fp_whole.shape[1], n_ndr)) 
    ct = 0
    prog0 = 0        
    for j in range(n_ndr):
        prog =  float(j)/float(n_ndr) 
        a =int(prog*10)  
        if a>prog0:
            print (int(prog*100), '% complete')
        prog0=int(prog*10)

        count = np.zeros_like(fp[int(osf/2.)::osf, int(osf/2.)::osf])   
        
        print ('ccc', count.size)
        for i in range(frames_per_ndr[j]): # frames per ndr
            for k in range(frame_osf): # jitter subframes per frame
      
                off_x = jitter_x[ct]
                off_y = jitter_y[ct]
                count0 = fp[int(osf/2.)+off_y::osf, int(osf/2.)+off_x::osf]
                print ('ccc0', count0.size)
                count += count0      
                ct=ct+1                
        ndr_count[...,j] = count /  (frames_per_ndr[j]*frame_osf)         
    return ndr_count
  
@jit(nopython=True, parallel=True)
def jitter(fp0,  osf, frames_per_ndr, frame_osf, jitter_x , jitter_y, start_list, buffer, fp):
    
    buffer_offset = buffer*osf

    n_ndr = len(frames_per_ndr)
    fp_whole = fp[int(osf/2.)::osf, int(osf/2.)::osf]
    ndr_count = np.zeros((fp_whole.shape[0], fp_whole.shape[1], n_ndr)) 
   
    for j in prange(n_ndr):

        count = np.zeros_like(fp[int(osf/2.)::osf, int(osf/2.)::osf])   
        for i in prange(frames_per_ndr[j]): # frames per ndr
            for k in prange(frame_osf): # jitter subframes per frame                                              
                start = start_list[j] 
                ct = start + i*frame_osf + k # for accurate parallelisation ct must be defined in terms of j, k and i
                off_x = jitter_x[ct] 
                off_y = jitter_y[ct] 
                start_x = buffer_offset + int(osf/2.)+ off_x 
                start_y = buffer_offset + int(osf/2.)+ off_y
                end_x = start_x + fp_whole.shape[1]*osf 
                end_y = start_y + fp_whole.shape[0]*osf 
                
                count0 = fp0[start_y: end_y: osf, start_x: end_x : osf]
                count += count0      
                ct=ct+1                
        ndr_count[...,j] = count /  (frames_per_ndr[j]*frame_osf)                

    return ndr_count

 
def jitter2(fp0,  osf, frames_per_ndr, frame_osf, jitter_x , jitter_y, start_list, buffer, fp):
      
    buffer_offset = buffer*osf

    n_ndr = len(frames_per_ndr)
    fp_whole = fp[int(osf/2.)::osf, int(osf/2.)::osf]
    ndr_count = np.zeros((fp_whole.shape[0], fp_whole.shape[1], n_ndr)) 
   
    for j in range(n_ndr):

        count = np.zeros_like(fp[int(osf/2.)::osf, int(osf/2.)::osf])   
        for i in range(frames_per_ndr[j]): # frames per ndr
            for k in range(frame_osf): # jitter subframes per frame                                              
                start = start_list[j] 
                ct = start + i*frame_osf + k # for accurate parallelisation ct must be defined in terms of j, k and i
                off_x = jitter_x[ct] 
                off_y = jitter_y[ct] 
                start_x = buffer_offset + int(osf/2.)+ off_x 
                start_y = buffer_offset + int(osf/2.)+ off_y
                
                end_x = start_x + fp_whole.shape[1]*osf 
                end_y = start_y + fp_whole.shape[0]*osf 
                
                count0 = fp0[start_y: end_y: osf, start_x: end_x : osf]
                # print (off_x, off_y, count0.shape)
                count += count0      
                ct=ct+1                
        ndr_count[...,j] = count /  (frames_per_ndr[j]*frame_osf)                

    return ndr_count


def simulate_jitter(opt):
    
  jitter_x = opt.osf*(opt.yaw_jitter_effective/opt.channel.plate_scale_x())
  jitter_y = opt.osf*(opt.pitch_jitter_effective/opt.channel.plate_scale_y())
   
  exosim_msg ("RMS jitter in pixel units x axis %s"%(np.std(jitter_x)/opt.osf), opt.diagnostics)
  exosim_msg ("RMS jitter in pixel units y axis %s"%(np.std(jitter_y)/opt.osf), opt.diagnostics)

  fp_units = opt.fp.unit
  fp  	   = opt.fp.value
  osf      = np.int32(opt.osf)
  offs     = np.int32(opt.offs)
    
  magnification_factor = np.ceil( max(3.0/jitter_x.std(), 3.0/jitter_y.std()) )
           
  if (magnification_factor > 1):
        try:
          mag = np.int(magnification_factor.item()) | 1
        except:
          mag = np.int(magnification_factor) | 1        
        
  exosim_msg ("mag %s"%(mag) , opt.diagnostics) 
  fp = exolib.oversample(fp, mag)
  osf *= mag
  offs = mag*offs + mag//2
  jitter_x *= mag
  jitter_y *= mag
  
  exosim_msg ("spatial osf %s"%(osf) , opt.diagnostics) 
    
  if opt.noise.EnableSpatialJitter.val ==0: 
          jitter_y *= 0.0
  if opt.noise.EnableSpectralJitter.val ==0:
          jitter_x *= 0.0
    
  jitter_x = np.round(jitter_x)
  jitter_y = np.round(jitter_y) 
    
  #important for speed of jitter module that these not have units
  jitter_x = jitter_x.value
  jitter_y = jitter_y.value
    
  exosim_msg('Applying jitter', 1)
     
  n_ndr = len(opt.frames_per_ndr)
  start_list =[]
  for j in range(n_ndr):
        if j ==0:
            start_list.append(0)
        else:            
            start_list.append(np.cumsum(opt.frames_per_ndr[0:j])[-1]*opt.frame_osf)
  start_list = np.array(start_list)

  aa = time.time()
  
  # cond =0
  # if int(fp.shape[0]/osf)*int(fp.shape[1]/osf) > 30000:
  #     print (int(fp.shape[0]/osf)*int(fp.shape[1]/osf))
  #     cond = 1 # splits array for jitter if above a certain size
 
  cond =1 # default to splitting as improved speed even for smaller arrays
  
  # edge buffer in whole pixels: needed to protect from errors arising when where jitter offset > half a pixel causing sampling to fall outside image area
  if jitter_y.max() > jitter_x.max():        
      buffer = int(np.ceil(jitter_y.max()/osf))    
  else:
      buffer = int(np.ceil(jitter_x.max()/osf)) 
  
  if cond == 0: # no spitting                   
      fp0 = np.zeros((fp.shape[0] + 2* buffer*osf, fp.shape[1] + 2* buffer*osf))
      fp0[buffer*osf: buffer*osf + fp.shape[0], buffer*osf: buffer*osf + fp.shape[1] ] = fp
                            
      noise =   jitter(fp0.astype(np.float32), 
        		 		osf.astype(np.int32),
        		 		opt.frames_per_ndr.astype(np.int32), 
        		 		opt.frame_osf.astype(np.int32),
        		 		jitter_x.astype(np.int32), 
        		 		jitter_y.astype(np.int32),
                        start_list.astype(np.int32),
                        buffer, 
                        fp.astype(np.float32)
                        ).astype(np.float32)  
 
       # noise =   create_jitter_noise(fp.astype(np.float32), 
       #     				osf.astype(np.int32),
       #     				opt.frames_per_ndr.astype(np.int32), 
       #     				opt.frame_osf.astype(np.int32),
       #     				jitter_x.astype(np.int32), 
       #     				jitter_y.astype(np.int32)).astype(np.float32)   
      
  if cond == 1: 
      exosim_msg('Splitting array for jitter code', 1)

      overlap = 5*osf
      x_size= int(5000/ (fp.shape[0]/osf) ) *osf
  
      n_div = int(fp.shape[1]/x_size)
      exosim_msg(f'Number of divisions {n_div}', opt.diagnostics)
      ct = 0
      
      # exosim_msg(f'x_size {x_size} osf {osf} fp shape {fp.shape}', opt.diagnostics)
      
      for i in range(n_div): 
          bb = time.time()
          
          ct = i*x_size
          
          if i == 0 and i == n_div-1:
             x0 = 0
             x1 = None
          elif i == 0:        
             x0 = ct
             x1 = ct + x_size + overlap
          elif i == n_div-1:
             x0 = ct - overlap
             x1 = None     
          else:
             x0 = ct - overlap
             x1 = ct + x_size + overlap
             
          # exosim_msg(f'x0 {x0} x1 {x1}', opt.diagnostics)

          div_fp= fp[:, x0:x1]
          
          # exosim_msg(f'div_fp shape {div_fp.shape}', opt.diagnostics)
          
          div_fp0 = np.zeros((div_fp.shape[0] + 2* buffer*osf, div_fp.shape[1] + 2* buffer*osf))
          div_fp0[buffer*osf: buffer*osf + div_fp.shape[0], buffer*osf: buffer*osf + div_fp.shape[1] ] = div_fp

          div_noise =   jitter(div_fp0.astype(np.float32), 
        		 		osf.astype(np.int32),
        		 		opt.frames_per_ndr.astype(np.int32), 
        		 		opt.frame_osf.astype(np.int32),
        		 		jitter_x.astype(np.int32), 
        		 		jitter_y.astype(np.int32),
                        start_list.astype(np.int32),
                        buffer, 
                        div_fp.astype(np.float32)
                        ).astype(np.float32)  
 
          if i == 0 and i == n_div-1:
              div_noise_stack = div_noise  
          elif i==0:
              div_noise = div_noise[:, 0:-int(overlap/osf), : ]
              div_noise_stack = div_noise
          elif i == n_div-1:
              div_noise = div_noise[:, int(overlap/osf): , : ]
              div_noise_stack = np.hstack((div_noise_stack,div_noise)) 
          else:
              div_noise = div_noise[:, int(overlap/osf):-int(overlap/osf), : ]
              div_noise_stack = np.hstack((div_noise_stack,div_noise))  
         
          # exosim_msg(f'div_noise_stack shape {div_noise_stack.shape}', opt.diagnostics)
          
          cc = time.time() -bb
          print(f'time remaining to complete jitter code {np.round((cc*(n_div-i)),2)} seconds', flush=True)
 
      exosim_msg(f'noise stack shape after split jitter code {div_noise_stack.shape}', opt.diagnostics)  
      noise  = div_noise_stack
                  
  qq = opt.frames_per_ndr* fp_units*opt.frame_time
  noise = noise*qq 

  exosim_msg('Time to run jitter code %s'%(time.time() - aa), opt.diagnostics) 
  
  return  noise 

#==============================================================================
#  Generates pointing timeline
#==============================================================================

def create_pointing_timeline(opt):
    # frames_per_ndr gives the number of frames in each ndr in sequence
    # ndr_end_frame_number gives the number of frames that have passed by the end of each NDR: take away the frames per ndr from ndr sequence to give a start index for the ndr in terms of frames
    # maccum = int(opt.effective_multiaccum)
 
    total_frames = int(np.sum(opt.frames_per_ndr))
    n_ndr = len(opt.frames_per_ndr)
    pointingArray = np.zeros((total_frames*opt.frame_osf, 3))    
    ct =0
    # for i in range (0, maccum*opt.n_exp, 1 ):
    for i in range (0, n_ndr, 1 ):
     
        idx = ct
        idx2 = int(idx + opt.frames_per_ndr[i]*opt.frame_osf)
        ct = idx2
    
        pointingArray[:,0][idx:idx2] = i        
        #This bit is very important for accuracy of jitter timeline > skips frames that fall in reset time
        start = int((opt.ndr_end_frame_number[i]- opt.frames_per_ndr[i])*opt.frame_osf)
        end =  int(start + opt.frames_per_ndr[i]*opt.frame_osf ) 
 
        pointingArray[:,1][idx:idx2] = opt.yaw_jitter[start:end]
        pointingArray[:,2][idx:idx2] = opt.pitch_jitter[start:end]
    

    return  pointingArray

 
 

def noise_simulator(opt): 

#==============================================================================
# Generates noisy image stack and returns pointing timeline   
#============================================================================== 

    #==============================================================================
    # set up jitterless or jittered pathways
    #==============================================================================

  fp = opt.fp[1::3,1::3] 
      
  if opt.noise.EnableSpatialJitter.val  ==1 or opt.noise.EnableSpectralJitter.val  ==1:
          exosim_msg ("using jittered array" , opt.diagnostics)
          opt.jitter_psd_file = opt.simulation.sim_pointing_model.val.replace('__path__', opt.__path__)
 
          exosim_msg ("psd file %s"%(opt.jitter_psd_file), opt.diagnostics)
          exosim_msg ("running jitter code", opt.diagnostics)
          try:
              opt.use_external_jitter
          except AttributeError:
             opt.use_external_jitter=0
          if opt.use_external_jitter==0:
              exosim_msg ("generating new jitter timeline...",  opt.diagnostics)
              opt.yaw_jitter, opt.pitch_jitter, opt.frame_osf = exolib.pointing_jitter(opt)
             
          elif opt.use_external_jitter==1:
              exosim_msg ("using external jitter timeline...", opt.diagnostics)
              opt.yaw_jitter, opt.pitch_jitter, opt.frame_osf = opt.input_yaw_jitter, opt.input_pitch_jitter, opt._input_frame_osf
       
          exosim_msg ("RMS jitter %s %s"%(np.std(opt.yaw_jitter), np.std(opt.pitch_jitter)  ) , opt.diagnostics)     
      
          pointing_timeline = create_pointing_timeline(opt) # this is also important to skip sections of jitter timeline that fall outside NDRs
              # the following takes into account skipped sections of jitter timeline due to reset groups
          opt.yaw_jitter_effective = pointing_timeline[:,1]*u.deg
          opt.pitch_jitter_effective = pointing_timeline[:,2]*u.deg              
          
          signal = simulate_jitter(opt)
            
  else:
     exosim_msg ("using jitterless array", opt.diagnostics)
     jitterless = np.ones((fp.shape[0], fp.shape[1], len(opt.frames_per_ndr)))
     jitterless =  np.rollaxis(jitterless,2,0)
     jitterless = jitterless*fp
     jitterless =  np.rollaxis(jitterless,0,3)  
     jitterless = jitterless*opt.frames_per_ndr*opt.frame_time      
     signal = jitterless
     pointing_timeline= np.zeros((opt.ndr_end_frame_number[-1], 3))

  if opt.timeline.apply_lc.val ==0:
     exosim_msg ("OMITTING LIGHT CURVE...", opt.diagnostics)
  else:
     exosim_msg ("APPLYING LIGHT CURVE...", opt.diagnostics)
     signal *= opt.lc 
     
     
  if opt.simulation.sim_use_systematic_model.val ==1:
 
     exosim_msg ("APPLYING SYSTEMATIC GRID...", opt.diagnostics)
     signal *= opt.syst_grid      


    #==============================================================================
    # add PRNU  APPPLY BEFOre PHOTON NOISE ADDED!!
    #==============================================================================
  if opt.noise.sim_prnu_rms.val ==0:
      opt.noise.ApplyPRNU.val =0

  opt.qe_grid = opt.qe # used for flat field     
  exosim_msg ("mean and standard deviation of flat field: should match applied flat in data reduction %s %s"%(opt.qe.mean(), opt.qe.std()),  opt.diagnostics  )
  exosim_msg ("standard deviation of flat field uncertainty %s"%(opt.qe_uncert.std()), opt.diagnostics  )     
  applied_qe= opt.qe*opt.qe_uncert 
  
  if  opt.noise.ApplyPRNU.val == 1:        
      exosim_msg ("PRNU GRID BEING APPLIED", opt.diagnostics  ) 
  else:
      exosim_msg ("PRNU GRID NOT APPLIED...", opt.diagnostics  )
      
    #  qe_expand  = qe.repeat(3,axis=0)
    #  qe_expand  = qe_expand.repeat(3,axis=1)
    #  qe_expand /=9.0
    #  c_qe_expand = exolib.fast_convolution(qe_expand, 6e-06*pq.m, opt.kernel, opt.kernel_delta) 
    #  qe_cross_talk = c_qe_expand[1::3,1::3]
      
       # exosim_msg ("Applied QE (PRNU) grid std (includes uncertainty) %s"%(applied_qe.std()), opt.diagnostics)
       # exosim_msg ("APPLYING PRNU GRID...", opt.diagnostics  )
       # signal =  np.rollaxis(signal,2,0)
       # signal = signal*applied_qe
       # signal =  np.rollaxis(signal,0,3)

    #==============================================================================
    # add backgrounds
    #==============================================================================
  blank_fp_shape = np.ones((fp.shape[0], fp.shape[1], len(opt.frames_per_ndr)))
  zodi_signal  = blank_fp_shape*0*u.electron
  emission_signal = blank_fp_shape*0*u.electron
  zodi_combined_noise = blank_fp_shape*0
  emission_combined_noise = blank_fp_shape*0
  combined_noise = blank_fp_shape*0
  
  if opt.background.EnableZodi.val == 1:
        exosim_msg ("APPLYING ZODI..." , opt.diagnostics)
        
        zodi_signal = opt.zodi.sed[opt.offs::opt.osf]
        zodi_signal = blank_fp_shape* zodi_signal[np.newaxis, :, np.newaxis]      
        zodi_signal = zodi_signal * opt.frame_time * opt.frames_per_ndr
          
        if  opt.noise.ApplyPRNU.val == 1:
            zodi_signal= ( zodi_signal.transpose() * applied_qe.transpose() ).transpose()
        
        zodi_signal = np.where(zodi_signal >= 0.0, zodi_signal, 0)
        
  else:
      exosim_msg ("NOT APPLYING ZODI..." , opt.diagnostics)

    
  if opt.background.EnableEmission.val == 1:
        exosim_msg ( "APPLYING EMISSION..."  ,  opt.diagnostics  )         
        emission_signal = opt.emission.sed[opt.offs::opt.osf]
        emission_signal = blank_fp_shape* emission_signal[np.newaxis, :, np.newaxis]       
        emission_signal = emission_signal * opt.frame_time * opt.frames_per_ndr
           
        if  opt.noise.ApplyPRNU.val == 1: 
            emission_signal= ( emission_signal.transpose() * applied_qe.transpose() ).transpose()
        
        emission_signal = np.where(emission_signal >= 0.0, emission_signal, 0)

  else:
      exosim_msg ("NOT APPLYING EMISSION..." , opt.diagnostics)
 
           
#  plt.figure('zodi_test')
#  plt.plot(opt.zodi.sed[opt.offs::opt.osf], 'bx')
#
#  plt.figure('emission_test')
#  plt.plot(opt.emission.sed[opt.offs::opt.osf], 'bx', markersize = 10)

    
# ==============================================================================
# Now do star signal
#============================================================================== 
  # I'm not sure if this is the right point to apply PRNU or whether it should be after all the poisson noise has been added.
  if  opt.noise.ApplyPRNU.val == 1:
      signal= (signal.transpose() * applied_qe.transpose() ).transpose()
  # signal is the star signal at this point witj jitter and lc, not emission or zodi photons
  signal = np.where(signal >= 0.0, signal, 0)
  
#==============================================================================
#   combine signal and noise upto this point
#==============================================================================

  signal = signal + zodi_signal + emission_signal

#==============================================================================
# add dark current and obtain dark current noise 
#==============================================================================
                 
  if opt.background.EnableDC.val ==1:
      exosim_msg ("DARK CURRENT being added...%s"%(opt.channel.detector_pixel.Idc.val ), opt.diagnostics  )
 
      dc_signal = blank_fp_shape*  opt.channel.detector_pixel.Idc() * opt.frame_time * opt.frames_per_ndr         
        
      signal = signal + dc_signal
      
  else:
      exosim_msg ("DARK CURRENT.. not... being added...", opt.diagnostics  )
  
#==============================================================================
# add combined shot noise to signal and account for ramps
#============================================================================== 
 
  if opt.noise.EnableShotNoise.val == 1:
      
       combined_noise = poission_noise(signal) # noise not signal+noise

       if opt.simulation.sim_full_ramps.val == 0 and opt.simulation.sim_use_UTR_noise_correction.val == 1:
             exosim_msg ("applying correction to photon noise for UTR read", opt.diagnostics  )
             n = opt.projected_multiaccum
             combined_noise, scale = poission_noise_UTR_correction(n, combined_noise)
         
  signal = signal.value + combined_noise 
  signal *=u.electron   
    
  signal = np.where(signal >= 0.0, signal, 0)
  
  exosim_msg ("check point 6.1 %s"%(signal.max()) , opt.diagnostics  )
  
    #==============================================================================
    # make ramps
    #==============================================================================

  for i in range(0, opt.n_ndr, opt.effective_multiaccum):
      signal[...,i:i+opt.effective_multiaccum] = np.cumsum(signal[...,i:i+opt.effective_multiaccum], axis=2)
 
  
    #==============================================================================
    # add read noise
    #==============================================================================

  if opt.noise.EnableReadoutNoise.val == 1:  
        exosim_msg ("READ NOISE... being added...", opt.diagnostics)
        if opt.simulation.sim_full_ramps.val == 0 and opt.simulation.sim_use_UTR_noise_correction.val == 1:             
            n = opt.projected_multiaccum
            exosim_msg ("applying correction to read noise for UTR read", opt.diagnostics  )
            signal = read_noise_UTR_correction(n, signal.value , opt.channel.detector_pixel.sigma_ro.val)
            signal = signal*u.electron    
        else:    
            signal = read_noise(signal.value, opt.channel.detector_pixel.sigma_ro.val)
            signal = signal*u.electron    

  else:
       exosim_msg ("READ NOISE...not... being added..." , opt.diagnostics  )
                
  exosim_msg ("check point 7 %s %s"%(signal.max(), signal.min()) , opt.diagnostics  )
  

  return signal, pointing_timeline

#==============================================================================
# generates a signal only stack (containing noiseless star signal)
#==============================================================================
def signal_simulator(opt):
    
  fp_signal = opt.fp_signal[1::3,1::3]
  
  signal = np.ones((fp_signal.shape[0], fp_signal.shape[1], len(opt.frames_per_ndr)))
  signal =  np.rollaxis(signal,2,0)
  signal = signal*fp_signal
  signal =  np.rollaxis(signal,0,3)

  signal = signal*opt.frames_per_ndr*opt.frame_time
       
  if opt.timeline.apply_lc.val ==1:
      signal *= opt.lc 
  
  signal = np.where(signal >= 0.0, signal, 0)   
   
  exosim_msg ("generating seperate noiseless signal array",  opt.diagnostics  )
   
    #==============================================================================
    # make ramps
    #==============================================================================

  for i in range(0, opt.n_ndr, opt.effective_multiaccum):
      signal[...,i:i+opt.effective_multiaccum] = np.cumsum(signal[...,i:i+opt.effective_multiaccum], axis=2)
 
          
  signal = signal*u.electron 
 
  return signal

#==============================================================================
# run code
#==============================================================================

def run(opt):
 
  
  opt.fp =   copy.deepcopy(opt.fp_original) # needed to reuse if cropped fp is used and monte carlo mode used
  opt.fp_signal =  copy.deepcopy(opt.fp_signal_original) # "
  
  opt.zodi.sed =  copy.deepcopy(opt.zodi_sed_original) # "
  opt.emission.sed =  copy.deepcopy(opt.emission_sed_original) # "

  opt.lc = copy.deepcopy(opt.lc_original)
  opt.syst_grid = copy.deepcopy(opt.syst_grid_original)
    
  # import matplotlib.pyplot as plt
  # plt.figure('syst_lc_ex2')
  # plt.plot(opt.lc[int(len(opt.x_wav_osr[1::3])/2)], 'b-')
  
  # xxxx
  
  opt.ldc = copy.deepcopy(opt.ldc_original)
  opt.cr_wl = copy.deepcopy(opt.cr_wl_original) 
  opt.cr = copy.deepcopy(opt.cr_original)  
  opt.x_wav_osr = copy.deepcopy(opt.x_wav_osr_original)
  opt.x_pix_osr = copy.deepcopy(opt.x_pix_osr_original)
  opt.qe = copy.deepcopy(opt.qe_original)
  opt.qe_uncert = copy.deepcopy(opt.qe_uncert_original)
  

  opt.data , opt.pointing_timeline = noise_simulator(opt)
  if opt.pipeline.useSignal.val == 1:
      opt.data_signal_only = signal_simulator(opt) 
   
 
 
#==============================================================================
# this is very important! 
  if opt.noise.EnableSpatialJitter.val== 0:
        opt.pointing_timeline[:,2] =  opt.pointing_timeline[:,2]*0
  if opt.noise.EnableSpectralJitter.val== 0:
        opt.pointing_timeline[:,1] =  opt.pointing_timeline[:,1]*0   
 
        
#==============================================================================    

  exosim_plot('focal plane check1', opt.diagnostics, image=True, 
                  image_data=opt.fp_signal[1::3,1::3], aspect='auto', interpolation = None,
                  xlabel = 'x \'spectral\' pixel', ylabel = 'y \'spatial\' pixel')
  
  exosim_plot('test - check NDR0', opt.diagnostics,
               image=True,  image_data = opt.data[...,0])
  exosim_plot('test - check NDR1', opt.diagnostics,
               image=True,  image_data = opt.data[...,1])
 

  return opt
   
  
