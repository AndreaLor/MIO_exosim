"""
exosim_n
 
Noise module

"""

import numpy as np
from exosim_n.lib.exosim_n_lib import exosim_n_msg, exosim_n_plot
from exosim_n.lib import exosim_n_lib, noise_lib, signal_lib 
from astropy import units as u
from exosim_n.classes.options import Options
import copy

def run(opt):
  DEBUG= Options.DEBUG
  if DEBUG:
        oldDiagnostic=opt.diagnostics
        opt.diagnostics=Options.showNoise

  opt.fp =   copy.deepcopy(opt.fp_original) # needed to reuse if cropped fp is used and monte carlo mode used
  opt.fp_signal =  copy.deepcopy(opt.fp_signal_original) # " 
  opt.zodi.sed =  copy.deepcopy(opt.zodi_sed_original) # "
  opt.emission.sed =  copy.deepcopy(opt.emission_sed_original) # "
  opt.lc = copy.deepcopy(opt.lc_original)    
  opt.ldc = copy.deepcopy(opt.ldc_original)
  opt.cr_wl = copy.deepcopy(opt.cr_wl_original) 
  opt.cr = copy.deepcopy(opt.cr_original)  
  opt.x_wav_osr = copy.deepcopy(opt.x_wav_osr_original)
  opt.x_pix_osr = copy.deepcopy(opt.x_pix_osr_original)
  opt.qe = copy.deepcopy(opt.qe_original)
  opt.qe_uncert = copy.deepcopy(opt.qe_uncert_original)

  opt.syst_grid = copy.deepcopy(opt.syst_grid_original)
   
  opt = signal_lib.initiate_signal(opt) 
 
  opt = signal_lib.apply_jitter(opt)
  
  print (opt.signal.max())
  import matplotlib.pyplot as plt
  plt.figure(333)
  plt.imshow(opt.signal[...,1].value)
  opt = signal_lib.apply_lc(opt)
  opt = signal_lib.apply_systematic(opt)
  
  opt = signal_lib.initiate_noise(opt)
  opt = signal_lib.apply_non_stellar_photons(opt)
  opt = signal_lib.apply_prnu(opt)
  opt = signal_lib.apply_dc(opt)
  print ('www', opt.signal.max())
  
  plt.figure('signal pixel level')
  plt.plot(opt.x_wav_osr[1::3], opt.signal[...,1].sum(axis=0), 'b-')
  plt.figure('signal ')
  plt.imshow(opt.signal[...,1].value)
  opt = signal_lib.apply_poisson_noise(opt)
  
  print ('non-jitter noise',opt.combined_noise.max()) 
  print ('After poisson', opt.signal.max())
  import matplotlib.pyplot as plt
  n= opt.combined_noise.sum(axis=0)
  n = n.std(axis=1)
  plt.figure('noise pixel level')
  plt.plot(opt.x_wav_osr[1::3], n, 'b-')
  plt.legend()
  print ('non-jitter noise',opt.combined_noise.max())
  n= opt.combined_noise.sum(axis=0)
  n = n.std(axis=1)
  plt.figure('noise pixel level')
  plt.plot(opt.x_wav_osr[1::3], n, 'r-')
  
  n= opt.combined_noise.sum(axis=0)
  n = n.std(axis=1)
  plt.figure('noise pixel level')
  plt.plot(opt.x_wav_osr[1::3], n, 'g-')

  print ('non-jitter noise',opt.combined_noise.max())

  print ('signal max', opt.signal.max())
  opt = signal_lib.apply_utr_correction(opt)
  
  plt.figure('noise pixel level')
  n= opt.combined_noise.sum(axis=0)
  n = n.std(axis=1)
  plt.plot(opt.x_wav_osr[1::3], n, 'r-')

  print ('signal max', opt.signal.max())
  print ('adding in non-jitter noise')
  opt = signal_lib.apply_combined_noise(opt)
  print ('signal max', opt.signal.max())
  
  opt = signal_lib.make_ramps(opt)

  opt = signal_lib.apply_read_noise(opt)

  opt.data = opt.signal
  opt.data_signal_only = opt.signal_only

  exosim_n_plot('focal plane check1', opt.diagnostics, image=True, 
                  image_data=opt.fp_signal[1::3,1::3], aspect='auto', interpolation = None,
                  xlabel = 'x \'spectral\' pixel', ylabel = 'y \'spatial\' pixel')
  
  exosim_n_plot('test - check NDR0', opt.diagnostics,
               image=True,  image_data = opt.data[...,0])
  exosim_n_plot('test - check NDR1', opt.diagnostics,
               image=True,  image_data = opt.data[...,1])
  if DEBUG:
      if opt.diagnostics:
        plt.show()
      else:
        plt.close(fig='all')
      opt.diagnostics=oldDiagnostic
  
  return opt
   
  
