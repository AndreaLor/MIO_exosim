# -*- coding: utf-8 -*-
"""
Created on Tue May 15 16:16:45 2018

@author: c1341133

Masking and binning

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


class processLC():
    def __init__(self, LC, binnedWav, binnedGamma, opt):
        self.opt = opt
        
        self.timegrid = opt.z_params[0]
        self.t0 = opt.z_params[1]
        self.per = opt.z_params[2]
        self.ars = opt.z_params[3]
        self.inc = opt.z_params[4]
        self.ecc = opt.z_params[5]
        self.omega = opt.z_params[6] 
 
#        self.opt.pipeline.useReduced.val =1
     
        try: 
            self.opt.LDClaw
        except AttributeError:
            self.opt.LDClaw ='quad'
            
        if self.opt.LDClaw !='claret4':       
           self.tm = QuadraticModel(interpolate=False)
        else:
           self.tm = SwiftModel(interpolate=False)          
        self.tm.set_data(self.timegrid)
                     

        self.LC = LC
        self.binnedWav = binnedWav
        self.binnedGamma = binnedGamma
        self.transitDepths = []
        self.transitDepths_err =[]
        self.model_gamma1 = []
        self.model_gamma2 = []
        self.model_gamma3 = []
        self.model_gamma4 = []    
        self.model_f = []
        self.finalWav = []
        self.fitModelLC()
        
        
        # if self.opt.pipeline.useReduced.val == 1:
        #     self.binnedWav = self.finalWav
               

    def fitModelLC_phot(self):

            self.dataLC = self.LC.reshape(1, self.LC.shape[0])[0]
            
            self.final_fit = self.light_curve_fit(0)
            p = self.final_fit[0]
            gamma1 = self.final_fit[2]
            gamma2 = self.final_fit[3]
            f = self.final_fit[1]
            self.transitDepths.append(p)
            self.model_gamma1.append(gamma1)
            self.model_gamma2.append(gamma2)
            self.model_f.append(f)
            
            if self.metadata['opt'].LDClaw == 'claret4':
                            gamma3 = self.final_fit[4]
                            gamma4 = self.final_fit[5]
                            self.model_gamma3.append(gamma3)
                            self.model_gamma4.append(gamma4)            
            
            print (p)


 
    def fitModelLC(self):
        
        nWLs = self.LC.shape[1]  # how many steps in the loop
    
        # Progress Bar setup:
        ProgMax = 100    # number of dots in progress bar
        if nWLs<ProgMax:   ProgMax = nWLs   # if less than 20 points in scan, shorten bar
        print ("|" +    ProgMax*"-"    + "|     Fitting model light curves")
        sys.stdout.write('|'); sys.stdout.flush();  # exosim_n_msg start of progress bar
        nProg = 0   # fraction of progress   

        step = 1
        # if self.opt.pipeline.useReduced.val == 1:
        #     exosim_n_msg ("using reduced number of R-bins for speed..." , self.opt.diagnostics)
        #     step = 5        
        
             
        for i in range(0, self.LC.shape[1], step):
            self.finalWav.append(self.binnedWav[i])
            
            if ( i >= nProg*nWLs/ProgMax ):
                    sys.stdout.write('*'); sys.stdout.flush();  
                    nProg = nProg+1
            if ( i >= nWLs-1 ):
                    sys.stdout.write('|     done  \n'); sys.stdout.flush(); 

            self.dataLC = self.LC[:,i]
            
            #fit light curve
            self.final_fit = self.light_curve_fit(i)
            
            p = self.final_fit[0]**2 #to give (Rp/Rs)^2
            f = self.final_fit[1]
            # p_err = self.final_fit[-1]
            self.transitDepths.append(p)  
            # self.transitDepths_err.append(p_err)
            self.model_f.append(f)   
            
            gamma1 = self.final_fit[2]
            gamma2 = self.final_fit[3]
            self.model_gamma1.append(gamma1)
            self.model_gamma2.append(gamma2)
            
            if self.opt.LDClaw == 'claret4':
                            gamma3 = self.final_fit[4]
                            gamma4 = self.final_fit[5]
                            self.model_gamma3.append(gamma3)
                            self.model_gamma4.append(gamma4)
                
            exosim_n_msg ('%s %s %s %s'%(i, len(self.transitDepths), p, self.binnedWav[i]),  self.opt.diagnostics)      

    
        
    def getModelLC_integrated(self, k, gamma):    
    
        modelLC_ = self.tm.evaluate(k=k, ldc=gamma, t0=self.t0, p=self.per, a=self.ars, i=self.inc, e=self.ecc, w=self.omega)   

        # now 'integrate' to match simulation
        lc = np.hstack((modelLC_, [0]))
        idx = np.hstack(([0],np.cumsum(self.opt.frames_per_ndr).astype(int)))   
        lc0 = np.add.reduceat(lc,idx)[:-1] /self.opt.frames_per_ndr
        lc0[lc0==np.inf] = 0
        # now 'do CDS' to match simulation
        lc_ndr_zero = lc0[0: len(lc0) : self.opt.effective_multiaccum]
        lc_ndr_final = lc0[self.opt.effective_multiaccum-1: len(lc0) : self.opt.effective_multiaccum]
        modelLC = lc_ndr_final - lc_ndr_zero
        
        return modelLC
        
    def getModelLC_instant(self, k, gamma):  

        
        modelLC_ = self.tm.evaluate(k=k, ldc=gamma, t0=self.t0, p=self.per, a=self.ars, i=self.inc, e=self.ecc, w=self.omega)   

        # now 'do CDS' to match simulation

        if self.opt.effective_multiaccum == 2:       
            lc0 = modelLC_
            lc0[lc0==np.inf] = 0
            modelLC =  modelLC_[1::2]
            
        else: 

            lc0 = modelLC_*1
            lc0[lc0==np.inf] = 0            
            
            n_exp = len(lc0)/self.opt.effective_multiaccum
            a = np.cumsum(lc0.reshape(n_exp, self.opt.effective_multiaccum), axis=1)
            t_ndr = self.opt.ndr_end_time[0:self.opt.effective_multiaccum]
            Sx = t_ndr.std()
            Sy = a.std(axis=1)
            R = (len(t_ndr)*(t_ndr*a).sum(axis=1) - t_ndr.sum()*a.sum(axis=1))  / \
            (np.sqrt( (len(t_ndr)*np.sum(t_ndr**2) - (np.sum(t_ndr))**2) * ( len(t_ndr)*(a**2).sum(axis=1) - (a.sum(axis=1))**2) ))  
            m = R*Sy/Sx
            modelLC   = m/m[0]          
                     
        return modelLC    
        
    def chi_sq (self, X):
        p = X[0]
        F = X[1]
        g1 = X[2]
        g2 = X[3]
      
        
        if self.opt.timeline.apply_lo_dens_LC.val == 1:
            model = self.getModelLC_instant(p, [g1,g2])* F 
        else:
            model = self.getModelLC_integrated(p, [g1,g2])* F  

        return np.sum(((model-self.dataLC))**2)
    
    def chi_sq_no_gamma (self, X, g1, g2):
        p = X[0]
        F = X[1]      
        
        if self.opt.timeline.apply_lo_dens_LC.val == 1:
            model = self.getModelLC_instant(p, [g1,g2])* F
        else:
            model = self.getModelLC_integrated(p, [g1,g2])* F  
        
        return np.sum(((model-self.dataLC))**2)
                
    def chi_sq_4 (self,X):
        p = X[0]
        F = X[1]
        g1 = X[2]
        g2 = X[3]      
        g3 = X[4]
        g4 = X[5]   
        
        if self.opt.timeline.apply_lo_dens_LC.val == 1:
            model = self.getModelLC_instant(p, [g1,g2, g3, g4])* F 
        else:
            model = self.getModelLC_integrated(p, [g1,g2, g3, g4])* F  

        return np.sum(((model-self.dataLC))**2)        
        
    def chi_sq_no_gamma_4 (self, X, g1, g2, g3, g4):
        p = X[0]
        F = X[1]    
        
        if self.opt.timeline.apply_lo_dens_LC.val == 1:
            model = self.getModelLC_instant(p, [g1,g2, g3, g4])* F 
        else:
            model = self.getModelLC_integrated(p, [g1,g2, g3, g4])* F  

        return np.sum(((model-self.dataLC))**2)        
        
        
    def light_curve_fit(self, i): 
        
        ex = self.dataLC
        
        err = np.std(np.hstack((ex[0:np.int(0.2*len(ex))], ex[np.int(0.8*len(ex)):])))
        oot_est = np.mean(np.hstack((ex[0:np.int(0.2*len(ex))], ex[np.int(0.8*len(ex)):])))
        it_est = np.mean(ex  [  int(len(ex)/2 -len(ex)/8)   :  int(len(ex)/2 +len(ex)/8)   ]     )
        cr_est = (oot_est - it_est ) / oot_est
        p_est = cr_est**0.5
        if err == 0:
            err = 1e-8      
         
        if self.opt.pipeline.fit_gamma.val == 1: 
            if self.opt.LDClaw=='claret4':
                fit_init = [p_est, oot_est, 0.0, 0.0, 0.0, 0.0]
                fit  = minimize(self.chi_sq_4, fit_init, args=(), method='Nelder-Mead', jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=None, callback=None, options=None) 
                final_fit = [fit['x'][0], fit['x'][1], fit['x'][2], fit['x'][3],  fit['x'][4], fit['x'][5]]  
            else:      
                fit_init = [p_est, oot_est, 0.0, 0.0]
                fit  = minimize(self.chi_sq, fit_init, args=(), method='Nelder-Mead', jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=None, callback=None, options=None) 
                final_fit = [fit['x'][0], fit['x'][1], fit['x'][2], fit['x'][3]]  
                 
        else:
            if self.opt.LDClaw=='claret4':
                fit_init = [p_est, oot_est]
                fit  = minimize(self.chi_sq_no_gamma_4, fit_init, args=(self.binnedGamma[0][i], self.binnedGamma[1][i], self.binnedGamma[2][i], self.binnedGamma[3][i]), method='Nelder-Mead', jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=None, callback=None, options=None) 
                final_fit = [fit['x'][0], fit['x'][1], self.binnedGamma[0][i], self.binnedGamma[1][i],self.binnedGamma[2][i], self.binnedGamma[3][i]]      
            
            else:                       
                fit_init = [p_est, oot_est]
                fit  = minimize(self.chi_sq_no_gamma, fit_init, args=(self.binnedGamma[0][i], self.binnedGamma[1][i]), method='Nelder-Mead', jac=None, hess=None, hessp=None, bounds=None, constraints=(), tol=None, callback=None, options=None) 
                final_fit = [fit['x'][0], fit['x'][1], self.binnedGamma[0][i], self.binnedGamma[1][i]]      
                            
  
        p_err = np.sqrt(2)* (err/oot_est)/(np.sqrt(self.opt.n_exp/2.)) #optional err estimate on p
        final_fit.append(p_err)

        if self.opt.diagnostics == 1:
            if i == int(self.LC.shape[1]/2):
        
                if self.opt.timeline.apply_lo_dens_LC.val == 1:
                    if self.opt.LDClaw=='claret4':
                        model_curvefit   =   self.getModelLC_instant(final_fit[0], [final_fit[2], final_fit[3], final_fit[4], final_fit[5]])*final_fit[1]      
                    else: 
                        model_curvefit   =   self.getModelLC_instant(final_fit[0], [final_fit[2], final_fit[3]])*final_fit[1]          
                else:
                    if self.opt.LDClaw=='claret4':
                        model_curvefit   =   self.getModelLC_integrated(final_fit[0], [final_fit[2], final_fit[3], final_fit[4], final_fit[5]])*final_fit[1]  
                    else:
                        model_curvefit   =   self.getModelLC_integrated(final_fit[0], [final_fit[2], final_fit[3]])*final_fit[1]
                   
                
                time = self.opt.ndr_end_time[self.opt.effective_multiaccum-1::self.opt.effective_multiaccum]
    
                plt.figure('LC fit light curve %s microns'%(np.round(self.binnedWav[i],2)))    
                plt.ylabel('count in spectral bin (electrons)')
                plt.xlabel('time (sec)')
                plt.plot(time, ex, 'r+-')
                plt.plot(time, model_curvefit, 'bx-')
                plt.figure('LC fit normalised light curve %s microns'%(np.round(self.binnedWav[i],2)))
                plt.ylabel('count in spectral bin (electrons)')
                plt.xlabel('time (sec)')
                plt.ylim(1-0.015,1.001)
                plt.plot(time, ex/oot_est, 'r+-')
                plt.plot(time, model_curvefit/final_fit[1], 'bx-')
     
        return final_fit
        


#==============================================================================
# Extraction of spectrum and binning
#==============================================================================

class extractSpec():
    
    def __init__(self, data, opt, diff, ApFactor, final_ap):
        self.data = data
        self.opt = opt
        self.diff = diff # note this must not be substituted by opt.diff, since will vary for signal only, noise etc
        self.ApFactor = ApFactor
        self.final_ap = final_ap
             
    #==============================================================================
    # Apply a mask based on wavelength and F and extract 1-D spectrum
    #==============================================================================
    def applyMask_extract_1D(self):
        
        wl = self.opt.cr_wl.value
        F = self.opt.channel.wfno.val
        pixSize = (self.opt.channel.detector_pixel.pixel_size.val).to(u.um).value
        ApFactor = self.ApFactor
        ApShape = self.opt.pipeline.pipeline_ap_shape.val
        wl_max = self.opt.channel.pipeline_params.wavrange_hi.val
        
        exosim_n_msg ("ap factor %s"%(ApFactor) , self.opt.diagnostics)
        exosim_n_msg ("ap shape %s"%(ApShape ),      self.opt.diagnostics)
        exosim_n_msg ("max wl %s"%(wl_max), self.opt.diagnostics)
     
#====================plot mask centre and edges========================================================== 

        w_list =[]
        for  i in range(len(wl)):
            if ApShape =='wav':                  
                    w = ApFactor*F*wl[i] 
            elif ApShape =='rect':  
                    w = ApFactor*F*wl_max  #defaults to rectangular
            w_list.append(w) # in distance units

        if self.diff !=1 :           
    #        self.opt.aa = self.data_signal_only.sum(axis =2)
            self.opt.aa = self.data.sum(axis =2)
             # 1) find max position of image
            indices = np.where(self.opt.aa == self.opt.aa.max())  
            y_max = indices[0].max() # max excludes rare situation of having 2 or more indices
            x_max = indices[1].max()
            ydata = self.opt.aa[:,x_max]
            xdata = np.arange(0,self.opt.aa.shape[0])            
            
            fitfunc  = lambda p, x: p[0]*np.exp(-0.5*((x-p[1])/p[2])**2)
            errfunc  = lambda p, x, y: (y - fitfunc(p, x))
            init  = [self.opt.aa[y_max,x_max], y_max , 2.0]
            out   = optimize.leastsq(errfunc, init, args=(xdata, ydata))[0]
    #                ymodel = out[0]*np.exp(-0.5*((xdata-out[1])/out[2])**2) 
            Cen0 = out[1]   # in pixel units     pixel unit = (distance unit/pixsize)  -0.5      
        else:
            Cen0 = (self.data.shape[0]/2.)-0.5
        
        exosim_n_msg ("Cen0 %s"%(Cen0),  self.opt.diagnostics)
                           
        X1 = Cen0 - np.array(w_list)/pixSize
        X2 = Cen0 +  np.array(w_list)/pixSize
        
        X0 = [Cen0]* len(w_list)
        
        if self.final_ap == 1:
            self.opt.cen = Cen0
        

        if self.final_ap ==2 or self.final_ap ==1:   #excludes signal only and n_pix runs 
            exosim_n_plot('y position of mask centre (pixels) vs x pixel', self.opt.diagnostics,
                         image = True,
                         image_data = self.data.sum(axis=2), aspect='auto', interpolation = None)    
        if self.final_ap ==2:   #2 = test ap chosen in noisy data 
            exosim_n_plot('y position of mask centre (pixels) vs x pixel', self.opt.diagnostics,
                         ydata = X1, marker = 'b--')
            exosim_n_plot('y position of mask centre (pixels) vs x pixel', self.opt.diagnostics,
                         ydata = X2, marker = 'b--')   
        if self.final_ap ==1:   #1 = final ap chosen in noisy data   
            exosim_n_plot('y position of mask centre (pixels) vs x pixel', self.opt.diagnostics,
                         ydata = X0, marker = 'c--',  linewidth =3)  
            exosim_n_plot('y position of mask centre (pixels) vs x pixel', self.opt.diagnostics,
                         ydata = X1, marker = 'w--',  linewidth =3)              
            exosim_n_plot('y position of mask centre (pixels) vs x pixel', self.opt.diagnostics,
                         ydata = X2, marker = 'w--',  linewidth =3)              


#==============================================================================

        nWLs = self.data.shape[2]  # how many steps in the loop
        # Progress Bar setup:
        ProgMax = 100    # number of dots in progress bar
        if nWLs<ProgMax:   ProgMax = nWLs   # if less than 20 points in scan, shorten bar
        
        if self.final_ap == 1:
            print ("|" +    ProgMax*"-"    + "|     Applying variable position mask: progress")
            sys.stdout.write('|'); sys.stdout.flush();  # exosim_n_msg start of progress bar
        nProg = 0   # fraction of progress   
               
    
        Cen_list=[]
        w_list =[]
        
        for i in range(self.data.shape[2]):
            if self.final_ap == 1:
                if ( i >= nProg*nWLs/ProgMax ):
     
                        sys.stdout.write('*'); sys.stdout.flush();  
                        nProg = nProg+1
                if ( i >= nWLs-1 ):
     
                        sys.stdout.write('|     done  \n'); sys.stdout.flush();                  
            inImage = self.data[...,i]
             # 1) find max position of image
            indices = np.where(inImage == inImage.max())  
            y_max = indices[0].max() # max excludes rare situation of having 2 or more indices
            x_max = indices[1].max()
            ydata = inImage[:,x_max]
            xdata = np.arange(0,inImage.shape[0])

#            if self.diff ==0:
#                fitfunc  = lambda p, x: p[0]*np.exp(-0.5*((x-p[1])/p[2])**2)
#                errfunc  = lambda p, x, y: (y - fitfunc(p, x))
#                init  = [inImage[y_max,x_max], y_max , 2.0]
#                out   = optimize.leastsq(errfunc, init, args=(xdata, ydata))[0]
##                ymodel = out[0]*np.exp(-0.5*((xdata-out[1])/out[2])**2) 
#                Cen = out[1]   #in pixel coords              
#            elif self.diff ==1:
#                Cen=Cen0*1.0
##                Cen = self.opt.aa.shape[0] /2. -0.5  
            
            Cen = Cen0*1.0
                            
            Cen_list.append(Cen)
            
            signal_list = [] 
            for j in range(inImage.shape[1]):                
                sig = inImage[: ,j]
                if ApShape =='wav':                  
                    w = ApFactor*F*wl[j] 
                elif ApShape =='rect':     
                    w = ApFactor*F*wl_max  #defaults to rectangular
                
                if i ==0:
                    w_list.append(2*w)
   
                X1 = Cen - w/pixSize
                X2 = Cen +  w/pixSize
                
                pixA = int(X1+0.5)
                pixB = int(X2+0.5)
                
                
                if pixA!= pixB:
                    wholepix = np.sum(sig[pixA+1:pixB])
                    RS = sig[pixB]*(X2+0.5-pixB)
                    LS = sig[pixA]*(pixA+0.5-X1)
                    in_mask_signal = wholepix +RS+LS
                    signal_list.append(in_mask_signal)
                                                   

                else:
                    in_mask_signal= sig[pixA]*(X2-X1)
                    signal_list.append(in_mask_signal)
                    
                
            if i == 0:
                signal_stack = signal_list
            else:
                signal_stack = np.vstack((signal_stack, signal_list))
 


        exosim_n_plot('Centre of mask (pixel units) vs exposure', self.opt.diagnostics, 
             ydata=Cen_list)    
        
        exosim_n_plot('width of mask (microns) vs pixel column', self.opt.diagnostics, 
             ydata=w_list, marker='bo')  
        
        exosim_n_plot('width of mask (microns) vs wavelength of pixel column', self.opt.diagnostics, 
             xdata=wl, ydata = w_list, marker = 'bo' )     
        
  
 
        self.spectra=signal_stack


        
    #==============================================================================
    # Extract 1-D spectrum
    #==============================================================================
        
    def extract1DSpectra(self):   
        
        spectra = np.zeros((self.data.shape[2],self.data.shape[1]))  
        # extract spectrum
        for i in range(self.data.shape[2]):
            spectra[i] = self.data[...,i].sum(axis = 0)
        
        self.spectra=spectra


    #==============================================================================
    # Binning of 1-D spectra into R- or fixed bins  
    #==============================================================================
        
    def binSpectra(self):    
        
        wl = self.opt.cr_wl.value
        R = self.opt.pipeline.pipeline_R.val
        wavrange_lo = self.opt.channel.pipeline_params.wavrange_lo.val
        wavrange_hi = self.opt.channel.pipeline_params.wavrange_hi.val
        wavrange = [wavrange_lo, wavrange_hi]
        x_wav_osr = self.opt.x_wav_osr
        x_pix_osr = self.opt.x_pix_osr
        pixSize = (self.opt.channel.detector_pixel.pixel_size.val).to(u.um).value
        bin_size = self.opt.pipeline.pipeline_bin_size.val

        cond=1 #directs to new method
#       cond=0 # previous method


#        1) find the bin sizes in wavelength space
        if self.opt.pipeline.pipeline_binning.val == 'R-bin' or  self.opt.pipeline.pipeline_binning.val == 'R' or  self.opt.pipeline.pipeline_binning.val == 'R_bin':
            exosim_n_msg ('binning spectra into R-bins...',  self.opt.diagnostics)
            # a) REMOVE ZEROS FROM EACH END OF WL SOLUTION
            for i in range (len(wl)):
                if wl[i]>0:
                    idx0 = i
                    break
            for i in range (len(wl)-1,0,-1):
                if wl[i]>0:
                    idx1 = i
                    break
            wl0 = wl[idx0:idx1]
       
            # xxxx
            # b) find w0, the starting wavelength  
            if wl0[-1] < wl0[0]:
                w0 = wl0[-1] 
            else:
                w0 = wl0[0]
                
            # c) calculate the sizes of each bin in microns of wavelength    
            dw = w0/(R-0.5)       
            bin_sizes=[dw]        
            for i in range(1000):
                dw2 = (1+1/(R-0.5))*dw
                bin_sizes.append(dw2)
                dw = dw2
                if np.sum(bin_sizes) > wavrange[1]-w0:
                    break
            bin_sizes = np.array(bin_sizes)   
                     
            
#         2) find the edges of each bin in wavelength units
            wavcen = w0+np.cumsum(bin_sizes)  # the central wavelength of each bin
      
            wavedge1 = wavcen-bin_sizes/2.   #  front edges
            wavedge2 =  wavcen+bin_sizes/2.  #  back edges
            wavedge = np.hstack((wavedge1[0],((wavedge1[1:]+wavedge2[:-1])/2.), wavedge1[-1])) # obtain an average value for each edge 
            # length of wavedge = length of wavcen +1 
            # print (len(wavedge), len(wavcen))
            # print (wavedge)
         
#         3)  find the bin edges in spatial units, in microns where 0 is at the left edge and centre of pixel is pixsize/2
#           # a) translate wl to x (microns)            
            wl_osr = x_wav_osr
            x_osr =  np.arange(pixSize/3./2., (pixSize/3.)*(len(x_wav_osr)), pixSize/3.)
            # this is the same as x_pix_osr but need the above if using cropped fp   
            
            # b) convert wavedge to xedge
            xedge = interpolate.interp1d(wl_osr,x_osr, kind ='linear', bounds_error = False)(wavedge)            
       
            
#           c) lay out pixel edge positions 
            x = np.arange(pixSize, pixSize*(len(wl)+1), pixSize)  # position of EDGES of pixels in microns starting at 0 on left side
                   
#           d) invert depending on wavelength solution 
            if wl0[-1]<wl0[0]:          
                xedge = xedge[::-1]  
                wavcen = wavcen[::-1]
                
            # print (len(wavcen), len(xedge))
         
            xedge = xedge[1:] #make len(xedge) = len(wavcen)
            # print (len(wavcen), len(xedge))     
 
         # d) remove nans
         
          
            idx  =  np.argwhere(np.isnan(xedge))
            xedge0 = np.delete(xedge, idx)           
            wavcen0 = np.delete(wavcen, idx) 
            
            # print (len(wavcen0), len(xedge0))  

#            So now we have A) edges of bins in wavelength units, B) edges of bins distance units : xedge0
#            C) edges of pixels in distance units : x 
            
            
 

#==============================================================================
# old code
#==============================================================================
            if cond == 0: 
                b_pix =[]
                wavcen_list=[]
                ct = 0
                

                for i in range(len(x)):
                     cond2 =0
                    
 
                     if ct == len(xedge0)-1:
                         break
                     if xedge0[ct+1]- xedge0[ct] < pixSize:
                         exosim_n_msg ("breaking since.. bin size < 1 pixel", self.opt.diagnostics)              
                         cond2 = 1               
                         break
                     
                     dist = x[i]
  
                     if dist > xedge0[ct]:
                         b_pix.append(i)
                         wavcen_list.append(wavcen0[ct])
                         ct = ct+1
                if cond2 ==1:        
                    exosim_n_msg ("min wav before subpixel binning starts: %s"%(np.min(wavcen_list)),   self.opt.diagnostics)            
                    
                for exp in range(self.spectra.shape[0]): 
                                              
                    spec  = self.spectra[exp]
                                    
                    count =[] 
                    count0 = 0
                    ct = 0
                    for i in range (len(wl)):
                        if ct >= len(b_pix):
                            break
                        if i == b_pix[ct]:
                            ct=ct+1
                            count.append(count0)
                            count0 =0
                        else:
                            count0 = count0 + spec[i]
                    x = np.arange(0, pixSize*(len(wl)), pixSize) + pixSize/2.
        
                    for i in range(len(b_pix)-1):               
                        x1 = b_pix[i]*pixSize
                        x3 = (b_pix[i]+1)*pixSize
                        x2 = xedge0[i]                    
                        f= interpolate.interp1d(x, spec , kind ='linear', bounds_error = False)(np.array([x1,x2,x3]))
                        if f[1]+f[0] ==0:
                             f[1] =  f[0] = 1e-30
                        if f[2]+f[1] ==0:
                             f[2] =  f[1] = 1e-30                                                           
    #                    print f
                        A1 = abs((x2-x1)*(f[1]+f[0])/2. )
                        A2 = abs((x3-x2)*(f[2]+f[1])/2. )
                                  
                        S = spec[b_pix[i]]                 
                        S1 = S*A1/(A1+A2)
                        S2 = S*A2/(A1+A2)
    #                    print S, S1, S2                   
                        count[i] = count[i] +S1
                        count[i+1] = count[i+1] +S2   

                    wavcen_list0 = wavcen_list   
                    
                    # now bin for subpixel bins if they exist
                    if cond2 ==1:
                        ct=0                 
                        idx = np.argwhere(wavcen0<np.min(wavcen_list))                   
                        xedge1 = xedge0[idx].T[0]
                        wavcen1 = wavcen0[idx].T[0]
                        xedgepix = xedge1/pixSize    
                                                                           
                        # first add the right side of the first divided pixel to the final count bin     
                        fracRight = xedgepix[ct]-int(xedgepix[ct])
                        SRight = spec[int(xedgepix[ct])]*fracRight
                        count[-1] = count[-1]+SRight 
                        
                        for j in range(len(wavcen1)-1):
                  
                        
                            if int(xedgepix[ct+1]) > int(xedgepix[ct]):
                                
                                if int(xedgepix[ct+1]) == int(xedgepix[ct]) +1:
    #                                print "span"
                                    fracLeft = 1-(xedgepix[ct]-int(xedgepix[ct]))
                                    SLeft =  spec[int(xedgepix[ct])]*fracLeft
                                    fracRight = xedgepix[ct+1]-int(xedgepix[ct+1])
                                    SRight = spec[int(xedgepix[ct +1])]*fracRight
                                    S= SLeft+ SRight
                                    count.append(S)                             
                                    ct= ct+1                                                                    
                                else:
                                    qq= int(xedgepix[ct])
                                    temp=0
                                    fracLeft = 1- (xedgepix[ct]-int(xedgepix[ct]))
    #                                xxxx
                                    SLeft =  spec[qq]*fracLeft   
                                    temp +=SLeft
                                    for i in range(1000):
                                        qq = qq+1
                                        S = spec[qq]
    #                                    print "qq", qq
                                        temp +=S
                                        if xedgepix[ct+1] < qq+2:
                                            fracRight = (xedgepix[ct+1]-int(xedgepix[ct+1]))
                                            SRight = spec[qq+1]*fracRight
    #                                        print SRight
                                            temp +=SRight
                                            count.append(temp)
                                            ct=ct+1
                                            break
                                                                                                          
                            else:
                                frac = xedgepix[ct+1]-xedgepix[ct]
                                S = frac*spec[int(xedgepix[ct])]
                                count.append(S)
                                ct=ct+1      
                                              
                    
                    
                    wavcen_list0 = wavcen0[:-1] 
    
                    
                    exosim_n_plot('binned spectrum', self.opt.diagnostics,
                                xdata = wavcen_list0, ydata = count, marker = 'bo') 
                    
           
                    if exp ==0:
                        count_array = count
                    else:
                        count_array = np.vstack((count_array, count))
#================ new code ==============================================================
# This deals with bins which are < 1 pixel

#=============   

            elif cond == 1:   
                for exp in range(self.spectra.shape[0]): 
                                          
                    spec  = self.spectra[exp]  # pick a 1 D spectrum
    
                    ct=0
                    count=[]
                    xedge1 = xedge0
                    wavcen1 = wavcen0 
                    xedgepix = xedge1/pixSize 
                    
                    for j in range(len(wavcen1)-1):
              
                        #selects if next bin edge is NOT in the same pixel
                        if int(xedgepix[ct+1]) > int(xedgepix[ct]):
                            
                            #selects if next bin edge is in the NEXT pixel
                            if int(xedgepix[ct+1]) == int(xedgepix[ct]) +1:
                                #signal from the left pixel
                                fracLeft = 1-(xedgepix[ct]-int(xedgepix[ct]))
                                SLeft =  spec[int(xedgepix[ct])]*fracLeft
                                #signal from the right pixel
                                fracRight = xedgepix[ct+1]-int(xedgepix[ct+1])
                                SRight = spec[int(xedgepix[ct +1])]*fracRight
                                # add these together
                                S= SLeft+ SRight
                                count.append(S)                             
                                ct= ct+1
                            
                            #selects if next bin edge is NOT in the NEXT pixel                                    
                            else:
                                qq= int(xedgepix[ct])
                                temp=0
                                #signal from the left pixel
                                fracLeft = 1- (xedgepix[ct]-int(xedgepix[ct]))
                                SLeft =  spec[qq]*fracLeft 
                                # add this to a cumulative count
                                temp +=SLeft
                                # move to the next pixel
                                for i in range(1000):
                                    qq = qq+1
                                    S = spec[qq]
                                    # add whole pixel count to cumulative
                                    temp +=S
                                    # check if next pixel has the bin edge
                                    if xedgepix[ct+1] < qq+2:
                                         # add the right pixel fraction to the count               
                                        fracRight = (xedgepix[ct+1]-int(xedgepix[ct+1]))
                                        SRight = spec[qq+1]*fracRight
                                        # final count for bin
                                        temp +=SRight
                                        
                                        count.append(temp)
                                        ct=ct+1
                                        break
                                                                                                      
                        else:
                            #selects if next bin edge is in SAME pixel
                            # find fraction of pixel in the bin 
                            frac = xedgepix[ct+1]-xedgepix[ct]
                            # add count
                            S = frac*spec[int(xedgepix[ct])]                           
                            count.append(S)
                            ct=ct+1
                    
                    wavcen_list0 = wavcen0[1:]                 
                    
                    exosim_n_plot('binned spectrum', self.opt.diagnostics,
                                 xdata = wavcen_list0, ydata=count, marker = 'bo-')
                   
                        
           
                    if exp ==0:
                        count_array = count
                    else:
                        count_array = np.vstack((count_array, count))
   
            self.binnedLC = count_array
            self.binnedWav = wavcen_list0
            
     
            exosim_n_msg('R-power obtained %s'%(np.mean(self.binnedWav/np.gradient(self.binnedWav))),self.opt.diagnostics)
           
            
         
        elif self.opt.pipeline.pipeline_binning.val  == 'fixed-bin' or self.opt.pipeline.pipeline_binning.val  == 'fixed' or self.opt.pipeline.pipeline_binning.val  == 'fixed_bin':
            exosim_n_msg ('binning spectra into fixed-bins of size %s pixel columns'%(bin_size), self.opt.diagnostics)
            
            offs=0          
#============================use only temp for noise budget to match spectral ==================================================
            if self.opt.channel.name == 'NIRSpec_G395M':
                 offs = 5
            if self.opt.channel.name == 'MIRI_LRS':
                 offs = 3
                   
#==============================================================================

            spec = np.add.reduceat(self.spectra, np.arange(int(offs),self.spectra.shape[1])[::int(bin_size)], axis = 1)
            wav = np.add.reduceat(wl, np.arange(offs,len(wl))[::int(bin_size)])  / bin_size            
#            spec = np.add.reduceat(self.spectra, np.arange(self.spectra.shape[1])[::bin_size], axis = 1)
#            wav = np.add.reduceat(wl, np.arange(len(wl))[::bin_size]) / bin_size
            if wav[-1] < wav [-2]:
                wav = wav[0:-2]
                spec = spec[:,0:-2]
                            
            self.binnedLC = spec
            self.binnedWav = wav

    #==============================================================================
    # Binning of the LDCs   
    #==============================================================================


    def binGamma(self):
        
        wl = self.opt.cr_wl.value
        R = self.opt.pipeline.pipeline_R.val
        wavrange_lo = self.opt.channel.pipeline_params.wavrange_lo.val
        wavrange_hi = self.opt.channel.pipeline_params.wavrange_hi.val
        wavrange = [wavrange_lo, wavrange_hi]
        x_wav_osr = self.opt.x_wav_osr
        x_pix_osr = self.opt.x_pix_osr
        pixSize = (self.opt.channel.detector_pixel.pixel_size.val).to(u.um).value
        bin_size = self.opt.pipeline.pipeline_bin_size.val
        useWeightedAv =1 #better than using simple average
#        useWeightedAv =0 
        
        if self.opt.pipeline.pipeline_binning.val == 'R-bin':
            exosim_n_msg ('binning LDCs into R-bins...', self.opt.diagnostics)
            for i in range (len(wl)):
                if wl[i]>0:
                    idx0 = i
                    break
            for i in range (len(wl)-1,0,-1):
                if wl[i]>0:
                    idx1 = i
                    break
            wl0 = wl[idx0:idx1]
            if wl0[-1] < wl0[0]:
                w0 = wl0[-1] 
            else:
                w0 = wl0[0]
                
            dw = w0/(R-0.5)       
            bin_sizes=[dw]
            for i in range(1000):
                dw2 = (1+1/(R-0.5))*dw
                bin_sizes.append(dw2)
                dw = dw2
                if np.sum(bin_sizes) > wavrange[1]-w0:
                    break
            bin_sizes = np.array(bin_sizes)    
            
            wavcen = w0+np.cumsum(bin_sizes)    
            wavedge1 = wavcen-bin_sizes/2.
            wavedge2 =  wavcen+bin_sizes/2.
            wavedge = np.hstack((wavedge1[0],((wavedge1[1:]+wavedge2[:-1])/2.), wavedge1[-1]))
     
            wl_osr = x_wav_osr
            x_osr =  np.arange(pixSize/3./2., (pixSize/3.)*(len(x_wav_osr)), pixSize/3.)
            xedge = interpolate.interp1d(wl_osr,x_osr, kind ='linear', bounds_error = False)(wavedge)            
            x = np.arange(pixSize, pixSize*(len(wl)+1), pixSize)   
  
            
            if wl0[-1]<wl0[0]:          
                xedge = xedge[::-1]  
                wavcen = wavcen[::-1]
            xedge = xedge[1:]   
            idx  =  np.argwhere(np.isnan(xedge))
            xedge0 = np.delete(xedge, idx)           
            wavcen0 = np.delete(wavcen, idx) 
        
            self.gamma = self.opt.ldc[1:]
            for jj in range(self.gamma.shape[0]):  
                                                
                spec  = self.gamma[jj]
                spec_flux  = np.sum(self.opt.fp_signal[1::3,1::3].value, axis=0) #used as weights for each coeff per pixel col.
                

#================ new code ==============================================================
# This deals with bins which are < 1 pixel
#=============          
                ct=0
                count=[]
                xedge1 = xedge0
                wavcen1 = wavcen0 
                xedgepix = xedge1/pixSize
                
                for j in range(len(wavcen1)-1):
          
                    
                    if int(xedgepix[ct+1]) > int(xedgepix[ct]):
                        
                        if int(xedgepix[ct+1]) == int(xedgepix[ct]) +1:
                            
                            if useWeightedAv !=1: 
                                fracLeft = 1-(xedgepix[ct]-int(xedgepix[ct]))
                                SLeft =  spec[int(xedgepix[ct])]*fracLeft                          
                                fracRight = xedgepix[ct+1]-int(xedgepix[ct+1])
                                SRight = spec[int(xedgepix[ct +1])]*fracRight
                                S= SLeft+ SRight
                                binsize = xedgepix[ct+1] - xedgepix[ct]
                                S= S/binsize
                            elif useWeightedAv ==1:  #weights each coeff by the flux in that pixel column
                                # Weighted average
                                fracLeft = 1-(xedgepix[ct]-int(xedgepix[ct]))
                                fracRight = xedgepix[ct+1]-int(xedgepix[ct+1])
                                SLeft_flux = spec_flux[int(xedgepix[ct])]*fracLeft
                                SRight_flux = spec_flux[int(xedgepix[ct])]*fracRight
                                SLeft =  spec[int(xedgepix[ct])]
                                SRight = spec[int(xedgepix[ct +1])]
                                S = (SLeft*SLeft_flux + SRight*SRight_flux)/ (SLeft_flux + SRight_flux)    

                            # add S                 
                            count.append(S)                             
                            ct= ct+1                           
                            
                                                            
                        else:
                            if useWeightedAv !=1: 
                                #simple average
                                qq= int(xedgepix[ct])
                                temp=0
                                binsize_temp = 0
                                fracLeft = 1- (xedgepix[ct]-int(xedgepix[ct]))                                
                                SLeft =  spec[qq]*fracLeft   
                                temp +=SLeft
                                binsize_temp += fracLeft
                                for i in range(1000):
                                    qq = qq+1
                                    S = spec[qq]
    #                              
                                    temp +=S
                                    binsize_temp += 1
                                    if xedgepix[ct+1] < qq+2:
                                                       
                                        fracRight = (xedgepix[ct+1]-int(xedgepix[ct+1]))
                                        SRight = spec[qq+1]*fracRight
    #                                        
                                        temp +=SRight
                                        binsize_temp += fracRight                      
                                        count.append(temp/binsize_temp)
                                        ct=ct+1
                                        break
                            elif useWeightedAv ==1:                                                                                                   
                                #weighted average : make two lists S for gamma*flux(i.e. weight) and flux_list for flux
                                flux_list =[] #weights
                                S_list = []
                                qq= int(xedgepix[ct])
                                fracLeft = 1- (xedgepix[ct]-int(xedgepix[ct]))
                                SLeft =  spec[qq]                             
                                SLeft_flux =  spec_flux[qq]*fracLeft
                                S_list.append(SLeft*SLeft_flux) #S*weight
                                flux_list.append(SLeft_flux)
                                for i in range(1000):
                                    qq = qq+1
                                    S = spec[qq]
                                    S_flux = spec_flux[qq]
                                    S_list.append(S*S_flux) #S*weight
                                    flux_list.append(S_flux)
                                    
                                    if xedgepix[ct+1] < qq+2:                                                  
                                        fracRight = (xedgepix[ct+1]-int(xedgepix[ct+1]))
                                        SRight = spec[qq+1]
                                        SRight_flux = spec_flux[qq+1]*fracRight
                                        
                                        S_list.append(SRight*SRight_flux) #S*weight
                                        flux_list.append(SRight_flux)
                                        weighted_av = np.sum(S_list)/np.sum(flux_list)                     
                                        count.append(weighted_av)
                                        ct=ct+1
                                        break
                                                                                                  
                    else:
#                        frac = xedgepix[ct+1]-xedgepix[ct]
                        S = spec[int(xedgepix[ct])]  # if bin size <1 pixel and within the pixel, no weighting needed                        
                        count.append(S)
                        ct=ct+1                              
                        
                exosim_n_plot('binned Gamma', self.opt.diagnostics,
                             xdata=self.opt.ldc[0], ydata=self.opt.ldc[jj+1], marker = 'ro-')
                exosim_n_plot('binned Gamma', self.opt.diagnostics,
                             xdata=self.binnedWav, ydata=count, marker ='bo-')         
    
  
           
                if jj ==0:
                    count_array = count
                else:
                    count_array = np.vstack((count_array, count))
   
            self.binnedGamma = count_array
            
         
        elif self.opt.pipeline.pipeline_binning.val  == 'fixed-bin':
            bin_size = int(bin_size)
            exosim_n_msg ('binning gamma into fixed-bins of size %s pixel columns...'%(bin_size), self.opt.diagnostics)
            self.gamma = self.opt.ldc[1:]
            for jj in range(self.gamma.shape[0]):            
                spec  = self.gamma[jj]           
                spec = np.add.reduceat(spec, np.arange(len(spec))[::bin_size]) / bin_size
                wav = np.add.reduceat(wl, np.arange(len(wl))[::bin_size]) / bin_size
                if wav[-1] < wav [-2]:
                    wav = wav[0:-2]
                    spec = spec[0:-2]                    
                if jj ==0:
                    count_array = spec
                else:
                    count_array = np.vstack((count_array, spec))
                            
            self.binnedGamma = count_array




#==============================================================================
# processing of OOT simulations                    
#==============================================================================
class processOOT():
          
     def __init__(self, LC, LC_signal, binnedWav, opt):
        self.LC = LC
        self.LC_signal = LC_signal
        self.binnedWav = binnedWav
        self.opt = opt      
        self.obtainSNR() 
          
        if self.opt.pipeline.useAllen.val ==1:
            for i in [0]:
                self.obtainAllen()
                                
                time_target = self.opt.T14.value
#                time_target = 3600               
                # idx = 0
                # if self.ootAllen[0].max() > time_target:
                #     idx = np.argwhere(self.ootAllen[0] > time_target)[0].item()
                #     idx = idx-50 # how many indicies from the end for fit
                #     if idx < 0:
                #         idx = 0
                # t1 = self.ootAllen[0][idx]
                               
                #  starting point
                idx = int(len(self.ootAllen[0])*0.25)
                t1 = self.ootAllen[0][idx]                                
                                                
                # x_fit = np.log10(self.ootAllen[0][idx:])
                # y_fit = np.log10(self.ootAllen[3][idx:])      
                # init = [-2.2,-0.5]
                # res_lsq = optimize.leastsq(self.linear2, init, args=(x_fit, y_fit))
                # c =  res_lsq[0][0]
                # m = res_lsq[0][1]
                # fit_time = np.arange(t1,time_target,10)
                # fit_y= 10**(c)*(fit_time)**(m)

                
                no= []
                for ii in range(self.ootAllen[2].shape[1]):                 
                    x_fit = np.log10(self.ootAllen[0][idx:])
                    y_fit = np.log10(self.ootAllen[2][:,ii][idx:])             
                    init = [-2.2,-0.5]
                    if y_fit[0]<1:
                        m_est = (y_fit[-1]-y_fit[0]) / (x_fit[-1]-x_fit[0])                        
                        init = [y_fit[0],-0.1]   
                    res_lsq = optimize.leastsq(self.linear2, init, args=(x_fit, y_fit))
                    c =  res_lsq[0][0]
                    m = res_lsq[0][1]
                    fit_time = np.arange(t1,time_target,10)
                    fit_y= 10**(c)*(time_target)**(m)*1e6
                    fit_y_= 10**(c)*(fit_time)**(m)
    
                    no.append(fit_y)
                    
                    if self.opt.diagnostics == 1:
                        plt.figure('Allan plots')
                        plt.loglog(self.ootAllen[0], self.ootAllen[2][:,ii] , 'o')
                        plt.plot(fit_time, fit_y_ , '--')
                        plt.ylabel('noise (ppm)')
                        plt.xlabel('binned time')
                        plt.grid(True)
                        
                    
                self.noiseAt1hr = np.array(no)
           
     def obtainSNR(self):
        self.ootSignal = self.LC_signal.mean(axis=0)
        self.ootNoise = self.LC.std(axis=0)
    
     def linear2(self, init, x,y):
        c = init[0]
        m = init[1]
        err = y**2 - (m*x+c)**2
        return err          
     def obtainAllen(self):
         
        timestep = self.opt.exposure_time.value
        binMax = int(self.LC.shape[0]/20  )
        for binSize in range (1, binMax):
            idx = np.arange(0, self.LC.shape[0],binSize)
            sig = (np.add.reduceat(self.LC,idx,axis =0)/ binSize)[:-1]
#            noise = sig.std(axis=0) / sig.mean(axis=0)
            noise = sig.std(axis=0) / self.LC_signal.mean(axis=0)
            if binSize ==1 :
#                sig_stack = sig.mean(axis=0)
                sig_stack = self.LC_signal.mean(axis=0)
                no_stack = noise
            else:
#                sig_stack = np.vstack((sig_stack,sig.mean(axis=0)))
                sig_stack = np.vstack((sig_stack,self.LC_signal.mean(axis=0)))            
                no_stack = np.vstack((no_stack,noise))
        no_median = np.median(no_stack, axis =1)              
        binTimes = np.arange(1, binMax)* timestep
        self.ootAllen = [binTimes, sig_stack, no_stack, no_median]



# =============================================================================
# Photometric                   
# =============================================================================
        
class extractPhot_decorr():             
       
     def __init__(self, data, metadata):
        self.data = data
        self.metadata = metadata
     
        self.applyAperture()  
        
#        self.photGamma = np.array([self.metadata['LDC'][1][0], self.metadata['LDC'][2][0]])
#        self.photGamma = self.photGamma.reshape(len(self.photGamma),1)

        self.photGamma = self.metadata['LDC'][1:,0]
        self.photGamma = self.photGamma.reshape(len(self.photGamma),1)

 
     def applyAperture(self):
        
        wl = self.metadata['WL'][0]
        F_ = self.metadata['F_y']
        ApFactor = self.metadata['ApFactor']
        diff = self.metadata['diff']
        pixSize = self.metadata['pixSize']
        
        
        F = F_*ApFactor     
        print ("aperture of radius = ", ApFactor, "x F lambda, where F = ", F_, "and lambda =", wl)
        
        signal = []
        subymaxlist = []
        subxmaxlist = []    
    
        nWLs =         self.data.shape[2]

    
        ProgMax = 100    # number of dots in progress bar
        if nWLs<ProgMax:   ProgMax = nWLs   # if less than 20 points in scan, shorten bar
        print ("|" +    ProgMax*"-"    + "|     Applyinng photometric aperture")
        sys.stdout.write('|'); sys.stdout.flush();  # print start of progress bar
        nProg = 0   # fraction of progress    

        if diff == 1: 
            print ("using fixed position due to diffuse noise.........")
        if diff ==0:
            for i in range(self.data.shape[2]):
                inImage = self.data[...,i]   
#                print inImage.max()/inImage.std(), inImage.max()/inImage.mean()
    #             if the max signal is less than 10 x noise in any image, it will probably not be dectable thus treat like diffuse.
                if inImage.max()/inImage.std() <10.0 or inImage.max()/inImage.mean() <10.0:
                    diff = 1                   
            if diff ==1:
                print ("using fixed position due to low SNR..........")
        

        for i in range(self.data.shape[2]):
    
            if ( i >= nProg*nWLs/ProgMax ):
                    '''Print dot at some fraction of the loop.'''
                    sys.stdout.write('*'); sys.stdout.flush();  
                    nProg = nProg+1
            if ( i >= nWLs-1 ):
                    '''If done, write the end of the progress bar'''
                    sys.stdout.write('|     done  \n'); sys.stdout.flush();     
        
    
            inImage = self.data[...,i]  
                  
            if diff == 0:
                
                def twoD_Gaussian(X, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
                    x = X[0]
                    y = X[1]
                    xo = float(xo)
                    yo = float(yo)    
                    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
                    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
                    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
                    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                                            + c*((y-yo)**2)))
                    return g.ravel()
                
                # Create x and y indices
    
                x = np.linspace(0, inImage.shape[1]-1, inImage.shape[1])
                y = np.linspace(0, inImage.shape[0]-1, inImage.shape[0])
                x, y = np.meshgrid(x, y)
                
                data_ravel = inImage.ravel()
                Y,X =  np.unravel_index(inImage.argmax(),inImage.shape)
                
                initial_guess = (inImage.max(),X,Y,10,10,0,0)
                popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), data_ravel, p0=initial_guess)
                
                data_fitted = twoD_Gaussian((x, y), *popt)
                data_fitted = data_fitted.reshape(inImage.shape[0], inImage.shape[1])
                
                Y,X =  np.unravel_index(data_fitted.argmax(),data_fitted.shape)
                 
                fact = 10. # 30 is too slow           
                y = np.arange(0,inImage.shape[0],1)
                y2 = np.arange(0,inImage.shape[0], 1/fact)
                x = np.arange(0,inImage.shape[1],1)
                x2 = np.arange(0,inImage.shape[0], 1/fact)
                
                f = interpolate.interp2d(x, y, data_fitted, kind='cubic')
                data_fitted_new = f(x2,y2)
                
                Y,X =  np.unravel_index(data_fitted_new.argmax(),data_fitted_new.shape)
                X /= fact
                Y /= fact
                
                subymax = Y
                subxmax = X
                     
            
            elif diff == 1:
                subymax = inImage.shape[0]/2
                subxmax = inImage.shape[1]/2
                 
            # apply circular aperture centred at ymax, xmax      
            radius = wl * F / pixSize  # F = F * position of 1st or 2nd minimum
            p = CircAp((subxmax,subymax), radius)
            Ap_Count = ApPhot(inImage,p)[0][3]
    
    #        print "radius of aperture in pixels", radius                
    ##                
    #        apertures = p
    #        plt.figure('aperture')
    #        plt.imshow(inImage, cmap='gray_r', origin='lower', interpolation = 'None')
    #        apertures.plot(color='blue', lw=1.5, alpha=0.5)
    #        xx
    #    
            signal.append(Ap_Count)
            
        self.photLC = np.array(signal).reshape(len(signal),1)
        self.photWav = wl

class extractPhot_fixed():             
       
     def __init__(self, data, metadata):
        self.data = data
        self.metadata = metadata
     
        self.applyAperture()        
 
     def applyAperture(self):
        
        wl = self.metadata['WL'][0]
        F_ = self.metadata['F_y']
        ApFactor = self.metadata['ApFactor']
        diff = self.metadata['diff']
        pixSize = self.metadata['pixSize']
              
        F = F_*ApFactor     
        print ("aperture of radius = ", ApFactor, "x F lambda, where F = ", F_, "and lambda =", wl)
        
        signal = []
        subymaxlist = []
        subxmaxlist = []    
    
        for i in range(self.data.shape[2]):
    
            inImage = self.data[...,i]  
                  
            def twoD_Gaussian(X, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
                x = X[0]
                y = X[1]
                xo = float(xo)
                yo = float(yo)    
                a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
                b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
                c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
                g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                                        + c*((y-yo)**2)))
                return g.ravel()
            
            # Create x and y indices

            x = np.linspace(0, inImage.shape[1]-1, inImage.shape[1])
            y = np.linspace(0, inImage.shape[0]-1, inImage.shape[0])
            x, y = np.meshgrid(x, y)
            
            data_ravel = inImage.ravel()
            Y,X =  np.unravel_index(inImage.argmax(),inImage.shape)
            
            initial_guess = (inImage.max(),X,Y,10,10,0,0)
            popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), data_ravel, p0=initial_guess)
            
            data_fitted = twoD_Gaussian((x, y), *popt)
            data_fitted = data_fitted.reshape(inImage.shape[0], inImage.shape[1])
            
            Y,X =  np.unravel_index(data_fitted.argmax(),data_fitted.shape)
             
            fact = 10. # 30 is too slow           
            y = np.arange(0,inImage.shape[0],1)
            y2 = np.arange(0,inImage.shape[0], 1/fact)
            x = np.arange(0,inImage.shape[1],1)
            x2 = np.arange(0,inImage.shape[0], 1/fact)
            
            f = interpolate.interp2d(x, y, data_fitted, kind='cubic')
            data_fitted_new = f(x2,y2)
            
            Y,X =  np.unravel_index(data_fitted_new.argmax(),data_fitted_new.shape)
            X /= fact
            Y /= fact
            
            subymax = Y
            subxmax = X
            
            subymaxlist.append(subymax)
            subxmaxlist.append(subxmax)
 
        subxmax = np.mean(subxmaxlist)
        subymax = np.mean(subymaxlist)
        
#        print np.std(subxmaxlist)
#        print np.mean(subxmaxlist)
#        
#        xxx

        for i in range(self.data.shape[2]):
            inImage = self.data[...,i]  
            # apply circular aperture centred at ymax, xmax      
            radius = wl * F / pixSize  # F = F * position of 1st or 2nd minimum
            p = CircAp((subxmax,subymax), radius)
            Ap_Count = ApPhot(inImage,p)[0][3]
    
#            print "radius of aperture in pixels", radius                
#    #                
#            apertures = p
#            plt.figure('aperture')
#            plt.imshow(inImage, cmap='gray_r', origin='lower', interpolation = 'None')
#            apertures.plot(color='blue', lw=1.5, alpha=0.5)
#            xx
#    #    
            signal.append(Ap_Count)
        
            
        self.photLC = np.array(signal).reshape(len(signal),1)
        self.photWav = wl

# =============================================================================
# Slit loss
# =============================================================================
        
class findSlitloss(): 

    def __init__(self, data, metadata):
        self.data = data
        self.metadata = metadata
        
        self.slitloss()
        
    def slitloss(self):
        # collapse in 1 D
        OneD = self.data.sum(axis=0)
#        plt.figure('test--')
#        plt.plot(OneD[:,0], 'rx')
        
        # interpolate to x100
        InImage = np.zeros((OneD.shape[0]*100, OneD.shape[1]))
        
        x = np.linspace(0,OneD.shape[0],OneD.shape[0]*100)
        xp = np.linspace(0,OneD.shape[0],OneD.shape[0])
        
        max_list =[]
        for i in range(OneD.shape[1]):
            fp  = OneD[:,i]
            InImage[:,i] = np.interp(x,xp,fp)
            max_list.append(np.argmax(InImage[:,i]))
#        plt.figure('test--') 
#        plt.plot(x, InImage[:,0], 'g--')
#        plt.plot(xp, OneD[:,0], 'rx')
        
        self.OneD_PSFs_unclipped = copy.deepcopy(InImage)
        
#        a =  InImage[:,0].sum()
        if self.metadata['opt'].psf_only_use_spectrum == 1 :             
            maxarg = int(np.mean(max_list))
        else:
            maxarg = self.metadata['opt'].cen_slit
        slit_width = np.int(np.round(self.metadata['slit_width']*100))
#        print 'slit width x 100', slit_width
        
        InImage[:maxarg- slit_width/2]= 0
        InImage[(maxarg- slit_width/2) + slit_width:]= 0
        
#        print ((maxarg- slit_width/2) + slit_width) -  (maxarg- slit_width/2)
#        plt.plot(x, InImage[:,0], 'b-')
#        plt.xlabel('x pixel')
#        plt.ylabel('count (e-)')
#        b= InImage[:,0].sum()
#        print "total signal before and after slit clip" , a,b
#        print "percent change in signal", 100*(b-a)/a
        
        self.OneD_PSFs_clipped = InImage # timeline of PSF profiles clipped by slit
        
        
          
    