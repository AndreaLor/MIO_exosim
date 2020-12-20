"""
ExoSim_N
Systematics module

"""

import numpy as np
from exosim_n.lib import exolib
from exosim_n.lib.exolib import exosim_msg
import copy
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
import os, sys



def add_systematic_grid(opt):
    
    wl = opt.x_wav_osr[1::3]
    
    time = opt.ndr_end_time
    
    syst_grid_name = opt.simulation.sim_systematic_model_name.val
    syst_grid_folder = '%s/exosim_n/data/systematic_models/%s'%(opt.exosim_path, syst_grid_name)
   
    # print ('%s/syst.txt'%(syst_grid_folder))
    
    # xxxxx
    if os.path.exists('%s/syst.txt'%(syst_grid_folder)):        
        syst_grid0 = np.loadtxt('%s/syst.txt'%(syst_grid_folder))     
        syst_wl = np.loadtxt('%s/wl.txt'%(syst_grid_folder))
        syst_time = np.loadtxt('%s/time.txt'%(syst_grid_folder))
    elif os.path.exists('%s/syst.npy'%(syst_grid_folder)):      
        syst_grid0 = np.load('%s/syst.npy'%(syst_grid_folder))  +1.0  
        syst_wl = np.load('%s/wl.npy'%(syst_grid_folder)) 
        syst_time = np.load('%s/time.npy'%(syst_grid_folder))
    else:
        exosim_msg('ERROR Systematics module 1: No systematic files found')
        sys.exit()
   
    f = interp2d(syst_time, syst_wl[::-1], syst_grid0, kind='linear')\
        
    # time = syst_time
    syst_grid = f(time, wl)
    
    plt.figure('syst_grid0')
    for i in range (len(wl)):
        if i == 0:
            plt.plot(time, syst_grid[i], 'r-', alpha = 0.5, markersize = 1,label = 'interpolated to ExoSim grid')
        else:
            plt.plot(time, syst_grid[i], 'r-', markersize = 1, alpha = 0.5)
    plt.figure('syst_grid0')
    for i in range (len(syst_wl)):
        if i == 0:
            plt.plot(syst_time, syst_grid0[i], 'b-', markersize = 1, alpha = 0.5, label = 'input file')
        else:
            plt.plot(syst_time, syst_grid0[i], 'b-', markersize = 1, alpha = 0.5)
        plt.xlabel('time (s)')
        plt.ylabel('fractional change')
    plt.legend()
    
    
    # plt.figure('syst_grid1')
    # for i in range (len(wl)):
    #     if i == 0:
    #         plt.plot(wl[i], syst_grid[i].std(), 'rx',  markersize = 1,label = 'interpolated to ExoSim grid')
    #     else:
    #         plt.plot(wl[i], syst_grid[i].std(), 'rx',  alpha = 0.5)
    # plt.figure('syst_grid2')
    std = syst_grid0.max(axis = 1)
    plt.figure('syst_grid1')
    plt.plot(syst_wl, std, 'bx',  alpha = 0.5)
    std2 = syst_grid.max(axis = 1)
    plt.figure('syst_grid1')
    plt.plot(wl, std2, 'rx',  alpha = 0.5)

    
  
    opt.syst_grid = syst_grid
    
    opt.syst_grid_original = copy.deepcopy(opt.syst_grid)
 
  

def run(opt):
    
    if opt.noise.ApplyRandomPRNU.val == 1:
        opt.qe = np.random.normal(1, 0.01*opt.noise.sim_prnu_rms.val, opt.fp_original[1::3,1::3].shape) # for random uncertainty
        opt.qe_uncert = np.random.normal(1, 0.01*opt.noise.sim_flat_field_uncert.val, opt.fp_original[1::3,1::3].shape) # for random uncertainty
        exosim_msg ("RANDOM PRNU GRID SELECTED...",  opt.diagnostics)
    else:
        opt.qe = np.load('%s/data/JWST/PRNU/qe_rms.npy'%(opt.__path__))[0:opt.fp_original [1::3,1::3].shape[0],0:opt.fp_original[1::3,1::3].shape[1]]
        opt.qe_uncert = np.load('%s/data/JWST/PRNU/qe_uncert.npy'%(opt.__path__))[0:opt.fp_original[1::3,1::3].shape[0],0:opt.fp_original[1::3,1::3].shape[1]]
        exosim_msg ("PRNU GRID SELECTED FROM FILE...", opt.diagnostics)

    opt.qe_original = copy.deepcopy(opt.qe)
    opt.qe_uncert_original = copy.deepcopy(opt.qe_uncert)

    if opt.simulation.sim_use_systematic_model.val == 1:
        add_systematic_grid(opt)
        
    return opt

      
    
