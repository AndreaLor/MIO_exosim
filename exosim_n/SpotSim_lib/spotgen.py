# -*- coding: utf-8 -*-
"""
Created on Wed May  2 06:36:18 2018

@author: subi


Generate spots and faculae pre-sims
and store in array
"""
import time
import numpy as np
import quantities as pq
import exosim_a as exosim
import spotsim, spotsim2
import matplotlib.pyplot as plt
from scipy import interpolate
import os
import gc


def run(opt, planet, var):

    RstarPx = 200.5  # radius of the star in pixels
    ImHalfSizePx = RstarPx+100 # = halfframe -0.5    
    
    Tstar = planet.planet.star.T.magnitude 
    CONTspot = np.round(-(1.343E-04*Tstar**2 - 6.849E-01*Tstar + 1.180E+03))
    CONTfac = np.round(-(87.4/663) *CONTspot)
    # from Muenier - 87.4 is the average over all angles of the fac temp and 663 is the spot contrast in their model.      
    Tspot = Tstar + CONTspot; Tfac = Tstar + CONTfac
    print "Star temp", Tstar
    print "Spot and faculae contrast", CONTspot, CONTfac    
    
    SUR = 6.5  
    opt.b = 0
    print "star temp", Tstar
    print "Spot and faculae contrast", CONTspot, CONTfac               
    Filling0 = opt.SpotSim_params['filling']
    

    variation = (var/100.)*Filling0
    print "Variation to filling factor std of", var, "%"
    opt.SpotSim_params['filling'] = np.random.normal(Filling0, variation)
    print "Input FF", Filling0, "Variation FF", opt.SpotSim_params['filling']
    Q = 10**-0.36368632*(opt.SpotSim_params['filling']/100.)**-0.51594621
    Q = np.round(Q,1)
    print "Q factor selected", Q  


    print "initiating spotsim"  
    Spotsim = spotsim2.SpotSim(Tstar, Tspot, Tfac, Q, ImHalfSizePx, RstarPx, SUR, opt)
    for i in range(10):
        Spotsim.size_dist(opt.SpotSim_params['filling'], SUR, opt.SpotSim_params['sizeDist'])
        Spotsim.spatial_dist(RstarPx, opt.SpotSim_params['spatDist']) 
        FF = Spotsim.create_spotted_disc()
        print "INITIAL FF",  FF
        if FF < 1.0*opt.SpotSim_params['filling']/100 and FF > 0.70*opt.SpotSim_params['filling']/100:
                break
        elif opt.SpotSim_params['useFileSpots'] == 1:
            break
        elif i == 9:
            print "sim failure >>>>>>>>>>>>>>>"
            xxxxxx
        else:
            print "Repeating spot simulation >>>>>>>>>>>>>>>"
            if FF <= 0.70*opt.SpotSim_params['filling']/100:
                print "FF too small"
            else:
                print "FF too big"
 
         
    spots, facs = Spotsim.calcStellarModel() 
#    plt.figure(4444)
#    plt.imshow(spots+facs*2)
#    xx
    
    return spots, facs


      
