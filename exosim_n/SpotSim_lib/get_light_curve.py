# -*- coding: utf-8 -*-
"""
Created on Wed May  2 07:20:47 2018

@author: subi


Generate array of light curves for a channel from SpotSim

"""
import numpy as np
import spotsim2
import matplotlib.pyplot as plt
from scipy import interpolate
 

def run(opt, planet, channel, key):
    
    RstarPx = 200.5  # radius of the star in pixels
    ImHalfSizePx = RstarPx+100 # = halfframe -0.5
    
   #------- set up z array for spotsim     
    
    print "impact parameter", opt.b               
    xArr= np.linspace(0, ImHalfSizePx*2-1,ImHalfSizePx*2)
    yArr = opt.b*RstarPx
    y4z =  yArr
    x4z = xArr-(ImHalfSizePx-0.5)
    zArr_s = ((x4z**2+y4z**2)**0.5)/(RstarPx)
    
   #------- calculate T_spot and T_fac     
    
    Tstar = planet.planet.star.T.magnitude 
    CONTspot = np.round(-(1.343E-04*Tstar**2 - 6.849E-01*Tstar + 1.180E+03))
    CONTfac = np.round(-(86.8/663) *CONTspot)
    # from Muenier - 87.4 is the average over all angles of the fac temp and 663 is the spot contrast in their model.      
    Tspot = Tstar + CONTspot; Tfac = Tstar + CONTfac
    print "Star temp", Tstar
    print "Spot and faculae contrast", CONTspot, CONTfac

   #------- calculate Q factor and set Spot-umbra ratio     

    SUR = 6.5 # ratio of total spot size to umbral size 
       
#        Y = 8.297E+01*ppm_ff**(7.720E-01) # our power law fit to Chapman - fails at high filling factors
#        Q = np.round((Y/ppm_ff),1)
#
#        Q = 10**-0.13087492*(filling/100)**-0.43475708 # our power law fit to Chapman + assumes Q of 1 at 50%)
#        Q = np.round(Q,1)
#        
    Q = 10**-0.36368632*(opt.SpotSim_params['filling']/100.)**-0.51594621
    Q = 0.4710*(opt.SpotSim_params['filling']/100.)**(-0.4915) #updated based on new plot in our paper (28/9/19)
    Q = np.round(Q,1)    
    print "Q factor selected", Q 
 
   #------- instantiate spotsim abd create spotted disk   
 
    Spotsim = spotsim2.SpotSim(Tstar, Tspot, Tfac, Q, ImHalfSizePx, RstarPx, SUR, opt)
      
    for i in range(10):     
#        Spotsim.size_dist(opt.SpotSim_params['filling'], SUR, opt.SpotSim_params['sizeDist'])
        Spotsim.size_dist_new(opt.SpotSim_params['filling'], SUR, opt.SpotSim_params['sizeDist'])
        
        
        
        Spotsim.spatial_dist(RstarPx, opt.SpotSim_params['spatDist']) 
        FF = Spotsim.create_spotted_disc()
        if FF < 1.0*opt.SpotSim_params['filling']/100 and FF > 0.70*opt.SpotSim_params['filling']/100:
                break
        elif opt.SpotSim_params['useFileSpots'] == 1:
            break
        elif i == 9:
            print "sim failure >>>>>>>>>>>>>>>"
            xxxxxx
        else:
            print "Repeating spot simulation >>>>>>>>>>>>>>>"


   #------- simulate transit and obtain light curve   

    cr    =  channel[key].planet.sed[channel[key].offs::channel[key].osf]
    cr_wl =  (channel[key].planet.wl[channel[key].offs::channel[key].osf]).magnitude
    limb_coeff = channel[key].limb_coeff
    zArr  = channel[key].z2
    high_timegrid= channel[key].high_timegrid 
    lc_exosim  = channel[key].lc2
      
    wl0, lc, lc_i, sim, simI, FF = Spotsim.calcLC(cr_wl, cr.max(), limb_coeff[1], limb_coeff[2], xArr)
    wl_s = cr_wl

        
    ##------ interpolate to exosim z and wl
        
    x_zArr = (zArr**2- zArr.min()**2) **0.5
    x_zArr_s = (zArr_s**2 - zArr_s.min()**2 )**0.5    
    
    idx = np.argmin(x_zArr)  
    if x_zArr[idx+1]- x_zArr[idx] > x_zArr[idx-1]- x_zArr[idx]:
        x_zArr[:idx] = -x_zArr[:idx]
    else:
        x_zArr[:idx+1] = -x_zArr[:idx+1]        
    ImHalfSizePx = np.int(ImHalfSizePx)                      
    x_zArr_s[:ImHalfSizePx] = -x_zArr_s[:ImHalfSizePx]
    f = interpolate.interp1d(x_zArr, high_timegrid)
    time_z = f(x_zArr_s)        
    
    if channel.keys()[0]=='FGS Prime' or channel.keys()[0]=='FGS Red' or channel.keys()[0]=='NIR Phot':
        f1 = interpolate.interp1d(x_zArr_s, lc[0], kind = 'linear',bounds_error=False, fill_value=0) 
        lc_spots = f1(x_zArr)                          
        f1 = interpolate.interp1d(x_zArr_s,lc_i[0], kind = 'linear',bounds_error=False, fill_value=0) 
        lc_spot_free = f1(x_zArr)    
        
        # fix to ensure outer points have a value
        for j in range(len(lc_spots)):
            if lc_spots[j] > 0:
                idx = j
                break                    
        idx0 = np.argwhere(lc_spots==0)
        idx1 = np.argwhere(lc_spot_free==0)                
        lc_spots[idx0] = lc_spots[idx]
        lc_spot_free[idx1] = lc_spot_free[idx]        
        
    else:

        f1 = interpolate.interp1d(x_zArr_s, lc, axis = 1, kind = 'linear', bounds_error=False, fill_value=0) 
        lc_spots = f1(x_zArr)
        f1 = interpolate.interp1d(x_zArr_s, lc_i, axis = 1, kind = 'linear', bounds_error=False, fill_value=0) 
        lc_spot_free = f1(x_zArr)
        for i in range(lc_spots.shape[0]):  # no correction where no cr
            if lc_exosim[i].min() == 1:     
                lc_spots[i] = lc_spot_free[i] = lc_spots[i]*0 + 1

        # fix to ensure outer points have a value               
        for i in range(lc_spots.shape[0]):
                for j in range(len(lc_spots[i])):
                    if lc_spots[i][j] > 0:
                        idx = j
                        break                    
                idx0 = np.argwhere(lc_spots[i]==0)
                idx1 = np.argwhere(lc_spot_free[i]==0)                
                lc_spots[i][idx0] = lc_spots[i][idx]
                lc_spot_free[i][idx1] = lc_spot_free[i][idx]
                
    lc_corr = (lc_spots/lc_spot_free) # must not normalise te curves since will lost OOT change in siganl due to background spots
        
    lc_final = lc_exosim*lc_corr

#    plt.figure ('check')  
#    plt.plot(lc_spots[10])
#    plt.plot(lc_spot_free[10])    
    return lc_final
#    
