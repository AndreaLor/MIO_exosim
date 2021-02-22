import numpy as np
from scipy import signal
from scipy import interpolate
from scipy.integrate import cumtrapz
import scipy.special
from astropy import units as u
from astropy import constants as const
from astropy.io import fits
import sys, os
import matplotlib.pyplot as plt
import pandas as pd
import pickle

 
def animate(Data):
    
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
 
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    wframe = None
    
    for j in range(0,Data.shape[2]):
    
            oldcol = wframe

            X = np.arange(0, Data.shape[1])
            Y = np.arange(0, Data.shape[0])
            X, Y = np.meshgrid(X, Y)
            
           
        
            Z = Data[...,j]
#            print Z.sum()
            
            wframe = ax.plot_wireframe(X, Y, Z, rstride=2, cstride=2)
            
#            ax.set_zlim(0,20000)
#            ax.set_title(j)
#            
     
        

    
        # Remove old line collection before drawing
            if oldcol is not None:
                ax.collections.remove(oldcol)
    
            plt.pause(.01)
            
            
def write_record(opt, path, lab, input_text_file):
    
    textfile = '%s/%s.txt'%(path,lab)
    
    with open(textfile, "w") as f1:
        f1.write('===== Simulation values =====')
        f1.write('\n ')
        f1.write('\nPlanet:  %s'%(opt.planet.planet.name))
        f1.write('\nChannel:  %s'%(opt.channel.name))
        if opt.channel.is_spec.val == True:    
            f1.write('\nType:  Spectrometer')
        else:
            f1.write('\nType:  Photometer')
        if opt.simulation.sim_mode.val == 3:
            f1.write('\n\nNoise option:  noise budget')          
        else:
            f1.write('\n\nNoise option:  %s'%(opt.noise_tag))     
        f1.write('\n ')

        f1.write('\nUse saturation time?:  %s'%(opt.observation.obs_use_sat.val))
        f1.write('\nSat time (to designated fraction of full well):  %s sec'%(opt.sat_time))
        f1.write('\nt_f:  %s sec'%(opt.t_f.value))
        f1.write('\nt_g:  %s sec'%(opt.t_g.value))
        f1.write('\nt_sim:  %s sec'%(opt.t_sim.value))
        
        f1.write('\nsaturation flag  :  %s'%(opt.sat_flag))
        f1.write('\nsaturation limit  :  %s electrons'%(opt.sat_limit.value))
        f1.write('\nnumber of saturated pixels per image  :  %s'%(opt.no_sat))
        f1.write('\nzero values applied to saturated pixels  :  %s'%(opt.pipeline.pipeline_bad_corr.val))
         
        f1.write('\nt_int:  %s sec'%(opt.t_int.value.item()))
        f1.write('\nt_cycle:  %s sec'%(opt.exposure_time.value.item()))  
        f1.write('\nprojected multiaccum (n groups):  %s'%(opt.projected_multiaccum))
        f1.write('\neffective multiaccum (n groups):  %s'%(opt.effective_multiaccum))
        f1.write('\nnumber of NDRs simulated:  %s'%(opt.n_ndr) )
        f1.write('\nnumber of integration cycles:  %s'%(opt.n_exp) )          
        f1.write('\n ')
        if opt.simulation.sim_output_type.val == 1: # excludes fits files
            f1.write('\nApFactor:  %s'%(opt.pipeline.pipeline_ap_factor.val) )
            f1.write('\nAperture shape:  %s'%(opt.pipeline.pipeline_ap_shape.val) )
            f1.write('\nSpectral binning:  %s '%(opt.pipeline.pipeline_binning.val) )
            if opt.pipeline.pipeline_binning.val == 'R-bin':
                        f1.write('\nBinned R power:  %s '%(opt.pipeline.pipeline_R.val) )
            else:
                      f1.write('\nBin size (pixels):  %s '%(opt.pipeline.pipeline_bin_size.val) )
        f1.write('\nWavelength: %s %s'%(opt.channel.pipeline_params.wavrange_lo.val, opt.channel.pipeline_params.wavrange_hi.val) )
    
    with open(textfile, "a") as f1:
        f1.write('\n ')
        f1.write('\n ')
        f1.write('===== Copy of input parameters file used =====')
        f1.write('\n ')
        f1.write('\n ')
      
    with open(input_text_file) as f:
        with open(textfile, "a") as f1:
            for line in f:       
                f1.write(line)

    f1.close()
    
    
def write_record_no_pipeline(opt, path, lab, input_text_file):
    
    textfile = '%s.txt'%(lab)
    
    with open(textfile, "w") as f1:
        f1.write('===== Simulation values =====')
        f1.write('\n ')
        f1.write('\nPlanet:  %s'%(opt.planet.planet.name))
        f1.write('\nChannel:  %s'%(opt.channel.name))
        if opt.channel.is_spec.val == True:    
            f1.write('\nType:  Spectrometer')
        else:
            f1.write('\nType:  Photometer')
        if opt.simulation.sim_mode.val == 3:
            f1.write('\n\nNoise option:  noise budget')          
        else:
            f1.write('\n\nNoise option:  %s'%(opt.noise_tag))     
        f1.write('\n ')

        f1.write('\nUse saturation time?:  %s'%(opt.observation.obs_use_sat.val))
        f1.write('\nSat time (to designated fraction of full well):  %s sec'%(opt.sat_time))
        f1.write('\nt_f:  %s sec'%(opt.t_f.value))
        f1.write('\nt_g:  %s sec'%(opt.t_g.value))
        f1.write('\nt_sim:  %s sec'%(opt.t_sim.value))
     
        f1.write('\nsaturation limit  :  %s electrons'%(opt.sat_limit.value))
        
        f1.write('\nt_int:  %s sec'%(opt.t_int.value.item()))
        f1.write('\nt_cycle:  %s sec'%(opt.exposure_time.value.item()))  
        f1.write('\nprojected multiaccum (n groups):  %s'%(opt.projected_multiaccum))
        f1.write('\neffective multiaccum (n groups):  %s'%(opt.effective_multiaccum))
        f1.write('\nnumber of NDRs simulated:  %s'%(opt.n_ndr) )
        f1.write('\nnumber of integration cycles:  %s'%(opt.n_exp) )          
        f1.write('\n ')
 
    with open(textfile, "a") as f1:
        f1.write('\n ')
        f1.write('\n ')
        f1.write('===== Copy of input parameters file used =====')
        f1.write('\n ')
        f1.write('\n ')
      
    with open(input_text_file) as f:
        with open(textfile, "a") as f1:
            for line in f:       
                f1.write(line)

    f1.close()
    
