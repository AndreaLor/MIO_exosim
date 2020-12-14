"""
ExoSim-N
Detector module

"""

from exosim_n.classes.sed import Sed
from exosim_n.lib import exolib 
from exosim_n.lib.exolib import exosim_msg, exosim_plot, planck
import numpy           as np
from astropy import units as u
import matplotlib.pyplot as plt
import copy, os
from scipy import interpolate
import sys
 
 
def run(opt):
    

#==============================================================================
#     get PSFs  
#==============================================================================
      if opt.simulation.sim_use_wfe.val == 0: 
            psf = exolib.Psf(opt.x_wav_osr.value, opt.channel.camera.wfno.val, opt.channel.wfno_y.val, opt.fp_delta.value, shape='airy')  
            psf[np.isnan(psf)] =0
            opt.psf_type = 'airy'
      elif opt.simulation.sim_use_wfe.val  == 1: 
            if opt.system == 'Ariel':
                zfile = '%s/../archive/PSF/%s/%s_psf_stack.fits'%(opt.__path__,opt.system,opt.channel.name)
                if opt.channel.is_spec.val == False:
                    psf = exolib.Psf_photometer(zfile, opt.fp_delta, opt.channel.osf(), opt.fpn, opt.x_wav_osr)   
                else:
                    psf = exolib.Psf_spectrometer(zfile, opt.fp_delta, opt.channel.osf(), opt.fpn, opt.x_wav_osr)  
            else:
                if os.path.exists('%s/../archive/PSF/%s_psf_stack.npy'%(opt.__path__,opt.channel.name)):
                  psf_stack = np.load('%s/../archive/PSF/%s_psf_stack.npy'%(opt.__path__,opt.channel.instrument.val))
                  psf_stack_wl =  1e6*np.load('%s/../archive/PSF/%s_psf_stack_wl.npy'%(opt.__path__,opt.channel.instrument.val))
                  psf = interpolate.interp1d(psf_stack_wl, psf_stack, axis=2,bounds_error=False, fill_value=0.0, kind='linear')(opt.x_wav_osr.value)      
                  psf = np.rot90(psf) 
                  opt.psf_type = 'wfe'
 

      opt.psf = psf 
      opt.channel.name


      exosim_msg("PSF shape %s, %s"%(opt.psf.shape[0],opt.psf.shape[1] ), opt.diagnostics)     
      exosim_plot('psf check', opt.diagnostics, image=True, image_data = psf[..., int(psf.shape[2]/2)])
    
      sum1=[]
      for i in range(psf.shape[2]):
          # print (psf[...,i].sum(), opt.x_wav_osr[i])
          sum1.append(psf[...,i].sum())
          if psf[...,i].sum()!=0:
               # psf[...,i] =psf[...,i]/psf[...,i].sum() 
                 if np.round(psf[...,i].sum(),3) !=1.0:    
                     exosim_msg('error... check PSF normalisation %s %s'%(psf[...,i].sum(), opt.x_wav_osr[i]), 1)
              
      exosim_plot('test7 - psf sum vs subpixel position (should be 1)', opt.diagnostics, 
                   xdata=opt.x_pix_osr, ydata = sum1, marker='bo')
    
#==============================================================================    
      # Populate focal plane with monochromatic PSFs
#============================================================================== 

      j0 = np.arange(opt.fp.shape[1]) - int(opt.psf.shape[1]/2)
      j1 = j0 + opt.psf.shape[1]
      idx = np.where((j0>=0) & (j1 < opt.fp.shape[1]))[0]
      i0 = np.array([opt.fp.shape[0]/2 - psf.shape[0]/2 + opt.offs]*len(j1)).astype(np.int)     
      i0+=1     
      
      # variable y position (applies only to NIRISS)
      if opt.y_pos_osr != []:
          opt.y_pos_osr = np.where(opt.y_pos_osr<0,0, opt.y_pos_osr).astype(np.int) +opt.fp.shape[0]/2 -opt.psf.shape[0]/2 + opt.offs 
          i0 = (opt.y_pos_osr).astype(np.int)
          
      i1 = i0 + psf.shape[0]
      

      if opt.background.EnableSource.val == 1 or opt.background.EnableAll.val == 1:  
          
          for k in idx: # actual signal 
              opt.fp[i0[k]:i1[k], j0[k]:j1[k]] += psf[...,k] * opt.star.sed.sed[k].value 
              
      for k in idx: # used for getting sat time, and sizing
          opt.fp_signal[i0[k]:i1[k], j0[k]:j1[k]] += psf[...,k] * opt.star.sed.sed[k].value 


#==============================================================================
#     Now deal with the planet
#==============================================================================
      a=0
      if a==1:
          opt.planet.sed =  Sed(opt.x_wav_osr, opt.planet_sed_original)
      else:              
          i0p = np.unravel_index(np.argmax(opt.psf.sum(axis=2)), opt.psf[...,0].shape)[0]
          planet_response = np.zeros((opt.fp.shape[1]))         
          for k in idx: 
              planet_response[j0[k]:j1[k]] += psf[i0p,:,k] * opt.planet.sed.sed[k].value          
          pl_sed = np.zeros((opt.fp.shape[1]))
          for i in range (len(planet_response)): 
               pl_sed[i] = planet_response[i]/(1e-30+ opt.fp[:,i][(i0[i]+i1[i])//2])                   
          opt.planet.sed =  Sed(opt.x_wav_osr, pl_sed*u.dimensionless_unscaled)
                
      exosim_plot('planet sed 1', opt.diagnostics, 
                   xdata=opt.planet.sed.wl, ydata = opt.planet.sed.sed, 
                   ylim=[0,1])

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
          
#==============================================================================
#     Allocate pixel response function and convolve with focal plane
#==============================================================================
 
      kernel, kernel_delta = exolib.PixelResponseFunction(opt, 
        opt.psf.shape[0:2],
        7*opt.channel.osf(),   
        opt.channel.detector_pixel.pixel_size(),
        lx = opt.channel.detector_pixel.pixel_diffusion_length.val)

      exosim_msg ("kernel sum %s"%(kernel.sum()), opt.diagnostics) 
      exosim_msg ("check 3.9 - unconvolved FP max %s"%(opt.fp.max()) , opt.diagnostics)      
      exosim_plot('test3', opt.diagnostics, 
                   xdata=opt.x_wav_osr, ydata = opt.fp.sum(axis=0), marker='bo')
      exosim_plot('test4', opt.diagnostics, 
                   xdata=opt.x_pix_osr, ydata = opt.x_wav_osr, marker='bo')                  
      exosim_plot('test5', opt.diagnostics, 
                   xdata=opt.x_pix_osr, ydata = opt.fp.sum(axis=0), marker='bo')
      exosim_plot('test6', opt.diagnostics, 
                   xdata=opt.x_wav_osr, ydata = opt.star.sed.sed, marker='bo') 
      
      opt.fp = exolib.fast_convolution(opt.fp, opt.fp_delta, kernel, kernel_delta)
      opt.fp_signal = exolib.fast_convolution(opt.fp_signal, opt.fp_delta, kernel, kernel_delta)

      exosim_msg ("check 4 - convolved FP max %s"%(opt.fp.max()), opt.diagnostics)    #FPCOPY0 = exolib.fast_convolution(FPCOPY[1::3,1::3], 18e-6*pq.m, kernel, kernel_delta)    

      # exosim_msg ("check 4 - convolved FP max %s %s"%(opt.fp.max(), FPCOPY[1::3,1::3].max()) , opt.diagnostics)    #FPCOPY0 = exolib.fast_convolution(FPCOPY[1::3,1::3], 18e-6*pq.m, kernel, kernel_delta)    
      opt.kernel = kernel
      opt.kernel_delta = kernel_delta
      
      # Fix units
      opt.fp = opt.fp*opt.star.sed.sed.unit 
      opt.fp_signal = opt.fp_signal*opt.star.sed.sed.unit  


#==============================================================================
#     Find saturation time
#==============================================================================
    ## Find count rate with diffuse radiation
      FPCOPY = copy.deepcopy(opt.fp_signal)
      FPCOUNT_no_bkg = FPCOPY[1::3,1::3] 
      FPCOPY += opt.zodi.sed  + opt.emission.sed 
      FPCOUNT = FPCOPY[1::3,1::3] 

      FPCOUNT += opt.channel.detector_pixel.Idc.val   
      FPCOUNT = FPCOUNT.value
      FPCOUNT_no_bkg = FPCOUNT_no_bkg.value
      
            
      exosim_msg ("check 5 - %s"%(FPCOUNT.max()), opt.diagnostics)

      FW = opt.channel.detector_pixel.full_well.val 
      
      A,B = np.unravel_index(FPCOUNT.argmax(), FPCOUNT.shape)
      exosim_msg ("maximum index and count with all backgrounds %s %s %s"%(A,B, FPCOUNT.max()), opt.diagnostics)
      
      A,B = np.unravel_index(FPCOUNT_no_bkg.argmax(), FPCOUNT_no_bkg.shape)
      exosim_msg ("maximum index and count with no backgrounds %s %s %s"%(A,B, FPCOUNT_no_bkg.max()), opt.diagnostics)
      
      exosim_msg ("full well in electron %s"%(FW), opt.diagnostics)   
      exosim_msg ("saturation time assuming 100 percent full well with all backgrounds %s"%((FW / FPCOUNT.max())), opt.diagnostics)     
      exosim_msg ("full well percentage chosen for saturation limit: %s"%(opt.observation.obs_fw_percent.val), opt.diagnostics)

      opt.sat_time = ((FW / FPCOUNT.max()) *opt.observation.obs_fw_percent.val/100.0 )*u.s
     
      opt.sat_time_no_bkg = ((FW / FPCOUNT_no_bkg.max()) *opt.observation.obs_fw_percent.val/100.0 )*u.s

      opt.sat_limit = u.electron*FW*opt.observation.obs_fw_percent.val/100.0  
      opt.sat_limit_fw = u.electron*FW
      
      exosim_msg ('saturation limit %s'%(opt.sat_limit), opt.diagnostics)
        
      exosim_msg ("saturation time with all backgrounds %s"%(opt.sat_time), opt.diagnostics)     
      exosim_msg ("saturation time with no backgrounds or dc %s"%(opt.sat_time_no_bkg), opt.diagnostics)
      
      opt.t_g = opt.t_f = opt.channel.array_read_time.val
      opt.dead_time = (opt.timeline.nGND.val+ opt.timeline.nRST.val)* opt.t_g
      opt.zero_time = opt.timeline.nNDR0.val* opt.t_g 
    
    
      if opt.observation.obs_use_sat.val == 1: 
          exosim_msg('Using saturation time to set n_groups', opt.diagnostics)
          n_groups = int(opt.sat_time/opt.t_g) # does not include reset group (assume this is after final read so saturation in this period does not affect read counts)
          if n_groups <2:
              n_groups=2
      else:
          exosim_msg('Using user-defined n_groups', opt.diagnostics)
          n_groups = opt.observation.obs_n_groups.val
          
      exosim_msg('t_f %s'%(opt.t_f), opt.diagnostics)
      exosim_msg('dead time %s'%(opt.dead_time), opt.diagnostics)
      exosim_msg('zero time %s'%(opt.zero_time), opt.diagnostics)
      exosim_msg('maccum %s'%(n_groups), opt.diagnostics)        
      
      # t_sim is not currently used, and by default is set to the same value as t_f
      opt.t_sim = opt.t_f             
      opt.t_int =  (n_groups-1)*opt.t_g   
      opt.t_cycle = n_groups*opt.t_g+ opt.dead_time
        
      if n_groups*opt.t_g > opt.sat_time:
          exosim_msg ("\n****** Warning!!!!!!!!!  : some pixels will exceed saturation limit ******\n", opt.diagnostics  )
          opt.sat_flag = 1
          # sys.exit()
           
      else:
          exosim_msg ("\n****** OK!!!!!!  Cycle time within saturation time ******\n", opt.diagnostics  )
          opt.sat_flag = 0
  
#==============================================================================
# 10.  Set effective multiaccum
#==============================================================================
      if opt.simulation.sim_full_ramps.val == 0:
          exosim_msg ("Approximating ramps with corrected CDS method, so only 2 NDRs simulated", 1)
          opt.effective_multiaccum = 2 # effective multiaccum is what is implemented in sim
          opt.projected_multiaccum = int(n_groups)
      else:
          opt.effective_multiaccum = int(n_groups)
          opt.projected_multiaccum = int(n_groups)
          
      exosim_msg ("projected multiaccum: %s"%(opt.projected_multiaccum), opt.diagnostics)
      exosim_msg ("effective multiaccum: %s"%(opt.effective_multiaccum), opt.diagnostics)
                                  
      opt.exposure_time = (opt.t_int + opt.dead_time + opt.zero_time) #same as t_cycle
      
      exosim_msg ("Integration time - zeroth read %s"%(opt.t_int), opt.diagnostics)  
      exosim_msg ("Estimated integration time incl. zeroth read %s"%(opt.t_int  + opt.zero_time), opt.diagnostics)
      exosim_msg ("Estimated TOTAL CYCLE TIME %s"%(opt.exposure_time), opt.diagnostics)
           
      exosim_msg ("CDS time estimate %s"%(opt.t_int), opt.diagnostics)  
        
      exosim_msg ("FP max %s"%(opt.fp[1::3,1::3].max()), opt.diagnostics)
      exosim_msg ("DC SWITCH.......%s"%(opt.background.EnableDC.val), opt.diagnostics)
      exosim_msg ("DISABLE ALL SWITCH.......%s"%(opt.background.DisableAll.val), opt.diagnostics)   
      exosim_msg ("DARK CURRENT %s"%(opt.channel.detector_pixel.Idc.val) , opt.diagnostics)
  
      exosim_plot('focal plane check', opt.diagnostics, image=True, 
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
  
      exosim_plot('final wl solution on subarray', opt.diagnostics,
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


      sanity_check(opt)
      
   
      return opt
# =============================================================================
#       
# =============================================================================
      
def sanity_check(opt):
    import scipy.constants as spc
    
    wl = opt.x_wav_osr[1::3]
    del_wl = abs(np.gradient(wl))
#    del_wl = opt.d_x_wav_osr[1::3]*3
    star_spec = opt.star_sed
    star_spec.rebin(wl)
    T = opt.planet.planet.star.T
    trans_sed = opt.total_transmission.sed*u.dimensionless_unscaled
    trans = Sed(opt.total_transmission.wl,trans_sed)
    trans.rebin(wl)
    QE = opt.qe_spec
    QE.rebin(wl)
    quantum_yield = 1
    Rs = (opt.planet.planet.star.R).to(u.m)
    D = (opt.planet.planet.star.d).to(u.m)
    n= quantum_yield* trans.sed*del_wl*np.pi*planck(wl,T)*(Rs/D)**2*opt.Aeff*QE.sed/(spc.h*spc.c/(wl*1e-6))

    n2= quantum_yield* trans.sed*del_wl*star_spec.sed*opt.Aeff*QE.sed/(spc.h*spc.c/(wl*1e-6))
    
    jex_sig = quantum_yield*opt.fp_signal[1::3,1::3].sum(axis=0)
    R = opt.channel.pipeline_params.pipeline_R.val
    del_wav = wl/R
    opt.exp_sig  = opt.t_int*del_wav*jex_sig/del_wl
    
    if opt.diagnostics ==1:
        plt.figure('sanity check 1 - check focal plane signal')
        plt.plot(wl,n, 'b^', label='BB check')
        plt.plot(wl,n2, 'r+', label='Phoenix check')  # not convolved with PSF unlike JexoSim, so peak may be higher
        plt.plot(wl, jex_sig, 'gx', label='JexoSim')
        plt.ylabel('e/s/pixel col'); plt.xlabel('pixel col wavelength (microns)')
        plt.legend(loc='best')
        
        plt.figure('sanity check 2 - expected final star signal in R bin of %s'%((R)))
        plt.plot(wl, opt.exp_sig)
        plt.ylabel('e/bin'); plt.xlabel('Wavelength (microns)')
            
        plt.figure('sanity check 3 - expected photon noise (sd) in R bin of %s'%((R)))
        plt.plot(wl, opt.exp_sig**0.5)    
        plt.ylabel('e/bin'); plt.xlabel('Wavelength (microns)')  
 
