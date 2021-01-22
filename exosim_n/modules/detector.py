"""
exosim_n
 
Detector module

"""
from exosim_n.lib import instrument_lib
from exosim_n.lib.exosim_n_lib import exosim_n_msg, exosim_n_plot
import numpy           as np
import matplotlib.pyplot as plt
import copy
 
def run(opt):
    
      opt.observation_feasibility = 1  # this variable is currently not changed to zero under any circumstances
      opt.psf = instrument_lib.get_psf(opt)
      opt.fp, opt.fp_signal = instrument_lib.get_focal_plane(opt)

      opt.planet.sed = instrument_lib.get_planet_spectrum(opt)
      opt.fp, opt.fp_signal = instrument_lib.convolve_prf(opt) #this step can gen some small negative values on fp
      opt.fp = np.where(opt.fp<0,0, opt.fp)
      opt.fp_signal = np.where(opt.fp_signal<0,0, opt.fp_signal)

      opt = instrument_lib.exposure_timing(opt)

      exosim_n_msg ("Integration time - zeroth read %s"%(opt.t_int), opt.diagnostics)  
      exosim_n_msg ("Estimated integration time incl. zeroth read %s"%(opt.t_int  + opt.zero_time), opt.diagnostics)
      exosim_n_msg ("Estimated TOTAL CYCLE TIME %s"%(opt.exposure_time), opt.diagnostics)         
      exosim_n_msg ("CDS time estimate %s"%(opt.t_int), opt.diagnostics)         
      exosim_n_msg ("FP max %s"%(opt.fp[1::3,1::3].max()), opt.diagnostics)
      exosim_n_msg ("DC SWITCH.......%s"%(opt.background.EnableDC.val), opt.diagnostics)
      exosim_n_msg ("DISABLE ALL SWITCH.......%s"%(opt.background.DisableAll.val), opt.diagnostics)   
      exosim_n_msg ("DARK CURRENT %s"%(opt.channel.detector_pixel.Idc.val) , opt.diagnostics)
      exosim_n_plot('focal plane check', opt.diagnostics, image=True, 
                  image_data=opt.fp[1::3,1::3], aspect='auto', interpolation = None,
                  xlabel = 'x \'spectral\' pixel', ylabel = 'y \'spatial\' pixel')
      if opt.diagnostics ==1:
          plt.figure('focal plane check')
          cbar = plt.colorbar()
          cbar.set_label(('Count (e$^-$/s)'), rotation=270, size=15,labelpad=20)
          cbar.ax.tick_params(labelsize=15) 
          ax = plt.gca()
          for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
              item.set_fontsize(15)

      opt.fp_original = copy.deepcopy(opt.fp)
      opt.fp_signal_original = copy.deepcopy(opt.fp_signal)  
      opt.x_wav_osr_original = copy.deepcopy(opt.x_wav_osr)
      opt.x_pix_osr_original = copy.deepcopy(opt.x_pix_osr) 

      opt.zodi_sed_original = copy.deepcopy(opt.zodi.sed) # needed here due to possible cropping above for subarrays
      opt.emission_sed_original = copy.deepcopy(opt.emission.sed) 

      exosim_n_plot('final wl solution on subarray', opt.diagnostics,
                   ydata=opt.x_wav_osr[1::3],
                   xlabel = 'x \'spectral\' pixel', ylabel = 'y \'spatial\' pixel',
                   grid=True)  
      if opt.diagnostics ==1:         
          wl = opt.x_wav_osr  
          wav = opt.x_wav_osr[1::3]         
          
          fp_1d = opt.fp_signal[1::3,1::3].sum(axis=0)         
          plt.figure('Focal plane pixel count rate')
          plt.plot(wav, fp_1d, 'r-')
          plt.grid()
          plt.figure('Final PCE and transmission on subarray')          
          plt.plot(wl, opt.PCE, '-') 
          plt.plot(wl, opt.TotTrans, '--')        
 
 
      instrument_lib.sanity_check(opt)
 
      return opt
