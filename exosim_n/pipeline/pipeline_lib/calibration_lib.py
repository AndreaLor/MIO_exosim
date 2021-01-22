# -*- coding: utf-8 -*-
"""
Created on Tue May 15 10:46:14 2018

@author: c1341133

Calibration

"""
import numpy as np
from exosim_n.lib.exosim_n_lib import exosim_n_msg
from astropy import units as u
import matplotlib.pyplot as plt
import pandas as pd

def fluxconvert(data, metadata):  
    cds_time = metadata['CDS_time']  
    file_path =  metadata['opt'].flux_conversion_file.replace('__path__', metadata['opt'].__path__)
    df = np.array(pd.read_csv(file_path,header=None))
    wl = df[:,0]
    X = df[:,1]
    
    data = data/cds_time
    
    for i in range(data.shape[2]):
        data[...,i] = data[...,i]*X
           
    return data

    # data now converted from e to W/m2/um     


def lincorr(data, metadata):
    
    print ( "correcting non-linearity...")
 
    aa = np.array(pd.read_csv('/Users/user1/Downloads/nonlin.csv', header=None))
    x = aa[:,0]*1e4
    y = aa[:,1]*1e4
  
    r = 3 # degree 7 works well        
    z1 = np.polyfit(y, x, r)
 
    x0 =0
    for i in range (0,r+1):
        x0 = x0 + z1[i]*data**(r-i) 
       
    aa = np.array(pd.read_csv('/Users/user1/Downloads/lin.csv', header=None))
    x = aa[:,0]*1e4
    y = aa[:,1]*1e4
    r = 1        
    z = np.polyfit(y, x, r)
     
    y0 =0
    for i in range (0,r+1):
        y0 = y0 + z[i]*x0**(r-i) 
        
    data = y0
    
    return data 


def calcAp(wav0, opt):
    #calculate aperture radius at wavlengths wav assuming 91% encircled energy - use only for Nirspec
    
    if opt.ch== 'NIR Spec':
#==============================================================================
#  91 % EE       
#==============================================================================
#        ApFactor = [4.41, 3.70, 3.36 , 3.04 , 2.86]
#==============================================================================
#  86% EE  
#==============================================================================      
        ApFactor = [5.45, 4.18, 3.32, 3.04, 2.58, 2.31 , 2.06 , 1.88]
        wl = [0.55, 0.70, 0.90, 1.00, 1.24 ,1.48 ,1.71 ,1.95]

#        r = 4 #  
#        z = np.polyfit(wl, ApFactor, r)
#        wavrange = np.linspace(0.9,2.0,100)
#        ApFactor_fit = 0
#        for i in range (0,r+1):
#            ApFactor_fit = ApFactor_fit + z[i]*wavrange**(r-i)  
#        
#        ApFactor0 = 0
#        for i in range (0,r+1):
#            ApFactor0 = ApFactor0 + z[i]*wav0**(r-i)      
#        wav = wav0
        
        log_wl = np.log10(wl)
        log_ApFactor = np.log10(ApFactor)
        r = 4 #  
        z = np.polyfit(log_wl, log_ApFactor, r)
        log_wavrange = np.linspace(np.log10(0.5),np.log10(2.0),100)
        log_ApFactor_fit = 0
        for i in range (0,r+1):
            log_ApFactor_fit = log_ApFactor_fit + z[i]*log_wavrange**(r-i) 
        wavrange = 10**log_wavrange
        ApFactor_fit = 10**log_ApFactor_fit

        log_wav = np.log10(wav0)        
        log_ApFactor0 = 0
        for i in range (0,r+1):
            log_ApFactor0 = log_ApFactor0 + z[i]*log_wav**(r-i) 
        ApFactor0 = 10**log_ApFactor0
        wav = 10**log_wav            
        if opt.diagnostics==1:             
            plt.figure('Ap factor calculation')
            plt.plot(wav, ApFactor0, 'gx' )
            plt.plot(wl, ApFactor,'ro' )
            plt.plot(wavrange, ApFactor_fit, 'b-' )

    else:
        
#==============================================================================
# if 91% EE        
#==============================================================================
        
#        wl = np.array([1.95, 3.0, 3.9, 5.9, 7.8])
#        ApFactor = np.array([2.99, 2.66, 2.57,2.48, 2.44])
#        
#        log_wl = np.log10(wl)
#        log_ApFactor = np.log10(ApFactor)
#        r = 4 #  
#        z = np.polyfit(log_wl, log_ApFactor, r)
#        log_wavrange = np.linspace(np.log10(1.9),np.log10(8.0),100)
#        log_ApFactor_fit = 0
#        for i in range (0,r+1):
#            log_ApFactor_fit = log_ApFactor_fit + z[i]*log_wavrange**(r-i) 
#        wavrange = 10**log_wavrange
#        ApFactor_fit = 10**log_ApFactor_fit
#                
#        log_wav = np.log10(wav0)        
#        log_ApFactor0 = 0
#        for i in range (0,r+1):
#            log_ApFactor0 = log_ApFactor0 + z[i]*log_wav**(r-i) 
#        ApFactor0 = 10**log_ApFactor0
#        wav = 10**log_wav
        
#==============================================================================
#   if 86 % EE
#==============================================================================

        wl = np.array([1.95, 3.0, 3.9, 5.9, 7.8])
        ApFactor = np.array([2.03, 1.66, 1.57, 1.48, 1.44])

        log_wl = np.log10(wl)
        log_ApFactor = np.log10(ApFactor)
        r = 4 #  
        z = np.polyfit(log_wl, log_ApFactor, r)
        log_wavrange = np.linspace(np.log10(1.9),np.log10(8.0),100)
        log_ApFactor_fit = 0
        for i in range (0,r+1):
            log_ApFactor_fit = log_ApFactor_fit + z[i]*log_wavrange**(r-i) 
        wavrange = 10**log_wavrange
        ApFactor_fit = 10**log_ApFactor_fit
        
        log_wav = np.log10(wav0)        
        log_ApFactor0 = 0
        for i in range (0,r+1):
            log_ApFactor0 = log_ApFactor0 + z[i]*log_wav**(r-i) 
        ApFactor0 = 10**log_ApFactor0
        wav = 10**log_wav             
        if opt.plots==1:
            plt.figure('Ap factor calculation')
            plt.plot(wav, ApFactor0, 'gx' )
            plt.plot(wl, ApFactor,'ro' )
            plt.plot(wavrange, ApFactor_fit, 'b-' )
     
    return ApFactor0
    
 

def loadData(opt):
    return opt    
    
def dqInit(data, opt):  
    dq_array = np.zeros(data.shape)
    for ii in range(data.shape[2]):
        dq_array[...,ii] = dq_array[...,ii] + opt.init_pix
    # flag 1 = initially known 'inoperable pixels'
#    plt.figure('dq array')
#    plt.imshow(dq_array[...,0])        
    opt.dq_array = dq_array
    return opt.dq_array     
    

def satFlag(data, opt): 
    
    margin = 0.1
    exosim_n_msg ("flagging pixels more than %s percent above saturation limit..."%(margin*100), opt.diagnostics)
    # margin accounts for fact that noise may increase the counts to above the designated sat limit

    sat_limit = opt.sat_limit.value + margin* opt.sat_limit.value

    idx = np.argwhere(data.value > sat_limit)

    dq_array = opt.dq_array
    for i in range (len(idx)):
        dq_array[idx[i][0]][idx[i][1]][idx[i][2]] = 2      
    # flag 2 = 'saturated pixels' (overmaps flag 1)

    exosim_n_msg ("number of saturated pixels over all NDRs %s"%(len(idx)),  opt.diagnostics)
    opt.dq_array = dq_array      
     
    return opt.dq_array     
 

def badCorr(data, opt): 
    # exosim_n_msg ("correcting bad pixels...", opt.diagnostics)
    exosim_n_msg ("applying zero value to saturated pixel timelines...", opt.diagnostics)
    opt.data_pre_flag = data*1
    dq_array = opt.dq_array
    bad_map = dq_array.sum(axis=2)
    idx = np.argwhere(bad_map > 0)
  
    exosim_n_msg('number of saturated pixels per image %s'%(len(idx )), opt.diagnostics)

    bad_map = np.where(bad_map == 0, 0, 1)  
    if opt.pipeline.pipeline_bad_corr.val == 1:
        for i in range (len(idx)):
            data[idx[i][0]][idx[i][1]] = 0  # make the entire time series of a saturated pixel = 0   
    
    opt.exp_image = data[...,0]
    opt.bad_map = bad_map
    opt.no_sat = len(idx)
    opt.data = data
    
    return opt      
    
      
def subZero(opt):
    exosim_n_msg ("subtracting zeroth read...", opt.diagnostics)
    multiaccum = opt.effective_multiaccum
    
    new_data = np.zeros((opt.data.shape[0],opt.data.shape[1],opt.data.shape[2]))
      
    for i in range (0, opt.data.shape[2], multiaccum):
        for j in range (0,multiaccum):
            new_data[...,i+j] = opt.data[...,i+j] - opt.data[..., i]
    opt.data = new_data        
    return opt   
        
def subDark(opt):    
    exosim_n_msg ("subtracting dark signal...", opt.diagnostics)
    
    dc_time =  opt.channel.detector_pixel.Idc.val* opt.duration_per_ndr
    new_data = opt.data - dc_time   
    opt.data = new_data
    return opt
        
    
def flatField(opt):    
    exosim_n_msg ("applying flat field...", opt.diagnostics)
    
    QE_grid = opt.qe_grid   
    opt.data =  np.rollaxis(opt.data,2,0)
    opt.data = opt.data/QE_grid
    opt.data =  np.rollaxis(opt.data,0,3) 
    
    exosim_n_msg("std of flat field...%s"%(QE_grid.std()), opt.diagnostics)
    
    return opt
    
 

def subBackground(opt) : 
    exosim_n_msg ("subtracting background...", opt.diagnostics)
   
    border_pix = 5
    if opt.data.shape[0]<20:
           border_pix = 1
    background_1 = opt.data[0:border_pix]
    background_2 = opt.data[-border_pix:]   
    background = np.vstack ( (background_1, background_2) )     
    aa = background.sum(axis=0)/background.shape[0]

    opt.data = opt.data- aa

    return opt


    
def doUTR(data, opt):
       
    multiaccum = int(opt.effective_multiaccum)
    n_exp = opt.n_exp
    time =  opt.ndr_end_time.value
    
    if multiaccum==2 or opt.simulation.sim_full_ramps.val == 0:
        exosim_n_msg ("doing CDS...", opt.diagnostics)

        new_data = np.zeros((data.shape[0],data.shape[1], n_exp))
        ct = 0        
        for i in range (multiaccum-1, data.shape[2], multiaccum):
            new_data[...,ct] = data[...,i]-data[...,i-multiaccum+1]
            ct = ct+1
        return new_data
    else:
        exosim_n_msg ("fitting ramps...", opt.diagnostics)
        utr_data = np.zeros(( data.shape[0], data.shape[1] ,n_exp))
#        utr_error = np.zeros(( data.shape[0], data.shape[1] ,n_exp))
   
        x = data.shape[1]/2
        y = data.shape[0]/2
        
        t_ndr = time[0:multiaccum]
    
        ct=0
        for i in range (0, n_exp*multiaccum, multiaccum):
              
            a = data[...,i:i+multiaccum]
            
#            Mx = t_ndr.mean()
#            My = a.mean(axis=2)
            Sx = t_ndr.std()
            Sy = a.std(axis=2)
            R = (len(t_ndr)*(t_ndr*a).sum(axis=2) - t_ndr.sum()*a.sum(axis=2))  / \
            (np.sqrt( (len(t_ndr)*np.sum(t_ndr**2) - (np.sum(t_ndr))**2) * ( len(t_ndr)*(a**2).sum(axis=2) - (a.sum(axis=2))**2) ))
             
            m = R*Sy/Sx
#            c = My - m*Mx
#            m_ = np.ones((m.shape[0], m.shape[1], len(t_ndr)))
#            c_ = np.ones((c.shape[0], c.shape[1], len(t_ndr))) 
#            
#            for j in range(len(t_ndr)):
#                m_[...,j] = m
#                c_[...,j] = c
#                
#            y_model = m_*t_ndr+c_
#            
#            residuals = ((y_model-a)**2).sum(axis=2)
#            
#            n = len(t_ndr)
#            D = sum(t_ndr**2) - 1./n * sum(t_ndr)**2
#            x_bar = np.mean(t_ndr)
#            dm_squared = 1./(n-2)*residuals/D
#            dc_squared = 1./(n-2)*(D/n + x_bar**2)*residuals/D
#            dm = np.sqrt(dm_squared)
#            dc = np.sqrt(dc_squared)
         
#             plot central pixel ramp and fit
#            print m[y][x], c[y][x], dm[y][x], dc[y][x],  R[y][x]
            
#            plt.figure('slope')
#            plt.plot(t_ndr, a[y][x], 'ro')
#            plt.plot(t_ndr, y_model[y][x], 'k--')
#            plt.grid(True)
#            plt.xlabel('time (s)')
#            plt.ylabel('electrons')
            
            utr_data[...,ct] = m
#            utr_error[...,ct] = dm
            ct +=1
        # if zero background used, then nan created with slopes - messes up later steps.
        utr_data[np.isnan(utr_data)] = 0
        utr_data[np.isinf(utr_data)] = 0  

#        utr_error[np.isnan(utr_error)] = 0 
        
        utr_data *= opt.t_int.value

#        return utr_data, utr_error   
          
        return utr_data    
  