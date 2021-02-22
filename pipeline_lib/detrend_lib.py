# -*- coding: utf-8 -*-
"""
Created on Tue May 15 16:16:45 2018

@author: c1341133

Detrending

"""

import numpy as np
from scipy import optimize
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy.optimize as opt
from photutils import CircularAperture as CircAp
from photutils import aperture_photometry as ApPhot
import sys
import uncertainties as uc
import copy

from scipy.optimize import minimize
from scipy.optimize import leastsq


import numpy as np
from scipy import interpolate, optimize
import matplotlib.pyplot as plt
from scipy import interpolate
import sys
import pytransit
# import pylightcurve
from scipy.optimize import minimize
from astropy import units as u
from pytransit import QuadraticModel, SwiftModel
from exosim_n.lib.exosim_n_lib import exosim_n_msg, exosim_n_plot


#==============================================================================
# Processing of transit light curves: fitting and extraction of transit depths
#==============================================================================


class Detrend():
    def __init__(self, binnedLC, binnedWav, opt):  
        self.opt = opt
        self.binnedLC = binnedLC
        self.binnedWav = binnedWav
        self.time = opt.exp_end_time_grid.value
        detrended_binnedLC = self.detrend_poly()
        
    def detrend_poly(self):
        for i in range(len(self.binnedWav)):
            plt.plot(self.time, self.binnedLC[:,i])
            r = 4 
                 
            z = np.polyfit(self.time,  self.binnedLC[:,i], r)
            p= np.poly1d(z)
            # yhat = p(wav)
            # ybar = sum(p_std)/len(p_std)
            # SST = sum((p_std - ybar)**2)
            # SSreg = sum((yhat - ybar)**2)
            # R2 = SSreg/SST  
            y =0
            for i in range (0,r+1):
                    y = y + z[i]*self.time**(r-i) 
                
            plt.plot(self.time, y, '--', color='r', linewidth=2) 
            plt.grid(True)
            
        xxx
        
       