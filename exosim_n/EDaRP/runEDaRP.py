# -*- coding: utf-8 -*-
"""
Created on Tue May 15 10:23:25 2018

@author: Subhajit Sarkar

Phase A ExoSim pipeline
"""
 

import numpy as np
from exosim_n.EDaRP import calibration, jitterdecorr, binning, detrend
from exosim_n.lib.exolib import exosim_msg, exosim_plot
from astropy import units as u  
 

class pipeline_stage_1():
       
    def __init__(self, opt):
        
        
        self.opt = opt
         
        self.exp_end_time_grid = self.opt.ndr_end_time[self.opt.effective_multiaccum-1::self.opt.effective_multiaccum]
        
        self.opt.exp_end_time_grid = self.exp_end_time_grid
        self.ndr_end_time_grid = self.opt.ndr_end_time
        self.opt.data_raw = opt.data*1
        
        self.ApFactor = opt.pipeline.pipeline_ap_factor.val

        if opt.background.EnableSource.val == 1:
            self.opt.diff =0
        else:
            self.opt.diff =1    
 
        self.loadData()  
        
        opt.init_pix = np.zeros(opt.data[...,0].shape) # placeholder - replace with initial bad pixel map 
        
        self.dqInit()  
        
        self.satFlag() 
        
        # if opt.use_nonlin ==1: 
        #         self.lincorr()
        
        if opt.background.EnableDC.val == 1: 
            if opt.diff==0:
                self.subDark()             
                
        if opt.noise.ApplyPRNU.val ==1: 
            exosim_msg ("mean and standard deviation of qe grid %s %s"%(opt.qe_grid.mean(), opt.qe_grid.std()), opt.diagnostics)
               
            self.flatField()
             
        if opt.background.EnableZodi.val ==1 or opt.background.EnableEmission.val ==1: 
            if opt.diff==0:                  
                self.subBackground()
        exosim_plot('sample NDR image', opt.diagnostics,
                     image =True, image_data=self.opt.data[...,1])
       
        self.doUTR() # lose astropy units at this step
      
        exosim_plot('sample exposure image', opt.diagnostics,
                     image =True, image_data=self.opt.data[...,1]) 
        
        
        # currently this applies zero values to saturated pixels
        self.badCorr()    # exp image obtained here     
   
        if opt.noise.EnableSpatialJitter.val  ==1 or opt.noise.EnableSpectralJitter.val  ==1 :
            exosim_msg ("Decorrelating pointing jitter...", opt.diagnostics)
            self.jitterDecorr()
            
        if opt.pipeline.pipeline_apply_mask.val==1:
            if opt.pipeline.pipeline_auto_ap.val == 1: 
                self.autoSizeAp()            
            
        self.extractSpec()
 
 
 
    def loadData(self):
        self.opt = calibration.loadData(self.opt)
        
    def dqInit(self):   
        self.dq_array = calibration.dqInit(self.opt.data, self.opt)        
        
    def satFlag(self):   
        self.dq_array = calibration.satFlag(self.opt.data, self.opt)  
        
    def subZero(self):   
        self.data = calibration.subZero(self.data, self.metadata)

    def subDark(self):
        self.opt = calibration.subDark(self.opt)

    def flatField(self):
        self.opt = calibration.flatField(self.opt)
   
    def subBackground(self):
        self.opt = calibration.subBackground(self.opt)  
        
    def doUTR(self):
        self.opt.data = calibration.doUTR(self.opt.data, self.opt) 
        
        if self.opt.pipeline.useSignal.val == 1:
            self.opt.data_signal_only = calibration.doUTR(self.opt.data_signal_only, self.opt)
            
    def badCorr(self):   
        self.opt = calibration.badCorr(self.opt.data, self.opt) 
                       
    def jitterDecorr(self):  
         
        self.jitterCode = jitterdecorr.jitterCode(self.opt.data, self.opt)
        self.opt.data = self.jitterCode.decorrData
          
    def driftCorr(self):      
        self.driftCode = jitterdecorr.jitterCode(self.data, self.metadata)
        self.data = self.jitterCode.decorrData        
        
    def fluxconvert(self):   
        self.data = calibration.fluxconvert(self.data, self.metadata)
        if self.opt.useSignal == 1:
            self.signal = calibration.fluxconvert(self.signal, self.metadata)

    def autoSizeAp(self):     
        F = self.opt.channel.wfno.val
        wl_max = self.opt.channel.pipeline_params.wavrange_hi.val 
        wl_min = self.opt.channel.pipeline_params.wavrange_lo.val
        wl = self.opt.cr_wl.value   
        if self.opt.timeline.apply_lc.val ==0 :       
            x=100
            if self.opt.data.shape[2] <x: # has to use the noisy data
                x= self.opt.data.shape[2]
            sample = self.opt.data[...,0:x]
        
        elif self.opt.timeline.apply_lc.val ==1 :     # must use OOT portions else the transit will contribute falsely to the stand dev
            if self.opt.pipeline.split ==0: # can only use the pre-transit values     
                total_obs_time =  self.opt.exposure_time*self.opt.n_exp
                pre_transit_time= self.opt.observation.obs_frac_t14_pre_transit.val * self.opt.T14   
                post_transit_time= self.opt.observation.obs_frac_t14_post_transit.val * self.opt.T14   
                pre_transit_exp = self.opt.n_exp*pre_transit_time/total_obs_time
                post_transit_exp = self.opt.n_exp*post_transit_time/total_obs_time
                x0 = int(pre_transit_exp*0.75)
                if x0 >20:
                    x0=20
                x1 = int(post_transit_exp*0.75)
                if x1 >20:
                    x1=20 
                sample1 =  self.opt.data[...,0:x0]
                sample2 =  self.opt.data[...,-x1:]
                sample = np.dstack((sample1,sample2))   
             
            if self.opt.pipeline.split ==1: # can only use the pre-transit values
                total_obs_time =  self.opt.exposure_time*self.opt.n_exp
                pre_transit_time= self.opt.observation.obs_frac_t14_pre_transit.val * self.opt.T14   
                pre_transit_exp = pre_transit_time/self.opt.exposure_time
                if self.opt.n_exp < pre_transit_exp:
                    x0 = int(self.opt.n_exp*0.75)
                else:
                    x0 = int(pre_transit_exp*0.75)
                if x0 > 40:
                    x0=40
                sample =  self.opt.data[...,0:x0]
                
            
        exosim_msg("SAMPLE SHAPE %s %s %s"%(sample.shape), self.opt.diagnostics)
        pix_size = (self.opt.channel.detector_pixel.pixel_size.val).to(u.um).value
        y_width = self.opt.data.shape[0] * pix_size
        testApFactorMax = int(int(y_width/2.) / (F*wl_max))
        if (testApFactorMax*F*wl_max/pix_size) + 1 >= self.opt.data.shape[0]/2:
            testApFactorMax = int( ((self.opt.data.shape[0]/2)-1)*pix_size/(F*wl_max) )

        ApList = np.arange(1.0,testApFactorMax+1,1.0)     
        exosim_msg ("maximum test ap factor %s"%(testApFactorMax), self.opt.diagnostics)
        for i in ApList:
            testApFactor = 1.0*i  
            exosim_msg("test ApFactor %s"%(testApFactor), self.opt.diagnostics)
            self.extractSample = binning.extractSpec(sample, self.opt, self.opt.diff, testApFactor, 2) 
 
            self.extractSample.applyMask_extract_1D()
            spectra = self.extractSample.spectra
            self.extractSample.binSpectra()   
            binnedLC = self.extractSample.binnedLC   
            wl =  self.extractSample.binnedWav
            SNR = binnedLC.mean(axis =0)/binnedLC .std(axis =0)
             
            if i == ApList[0] :
                SNR_stack = SNR
            else:
                SNR_stack = np.vstack((SNR_stack, SNR))
            
            exosim_plot('test aperture SNR', self.opt.diagnostics, 
                        xdata=wl,ydata=SNR, label = testApFactor)   
        
        idx = np.argwhere( (wl>=wl_min) & (wl<= wl_max))
        SNR_stack = SNR_stack[:,idx][...,0]
        
        best=[]
        for i in range(SNR_stack.shape[1]):
            aa = SNR_stack[:,i]
#            exosim_msg i, np.argmax(aa), ApList[np.argmax(aa)]
            best.append(ApList[np.argmax(aa)])          
        wl=  wl[idx].T[0]

        exosim_plot('highest SNR aperture vs wavelength', self.opt.diagnostics,
                     xdata=wl, ydata=best, marker='bo-')      
        
        AvBest = np.round(np.mean(best),0)
        exosim_msg ("average best aperture factor %s"%(AvBest), self.opt.diagnostics)
        self.opt.AvBest = AvBest
        self.ApFactor = AvBest
        self.opt.pipeline.pipeline_ap_factor.val = self.ApFactor
 
    def extractSpec(self):
        
        exosim_msg ("extracting 1 D spectra....", self.opt.diagnostics)
        #==============================================================================
        # 1) initialise class objects a) noisy data, b) noiseless data if used, c) extraction of n_pix per bin
        #==============================================================================
        
        self.extractSpec = binning.extractSpec(self.opt.data, self.opt, self.opt.diff, self.ApFactor,1)
        # 1 = final ap : means code knows this is the ap factor to use on noisy data and plots
        # final ap  0 = signal only or n_pix evals, so no figures for aperture on image plotted.
        if self.opt.pipeline.useSignal.val == 1:        
            self.extractSpec_signal = binning.extractSpec(self.opt.data_signal_only, self.opt, 0, self.ApFactor, 0) #diff always 0 for signal
        
        # use a set of 3 images with 1 value for all pixels to find out no of pixels per bin
        self.extractSpec_nPix = binning.extractSpec(self.opt.data[...,0:3]*0+1, self.opt, 1, self.ApFactor, 0)
         
        #==============================================================================
        # 2) apply mask and extract 1D spectrum if apply mask selected      
        #==============================================================================
            
        if self.opt.pipeline.pipeline_apply_mask.val==1:             
            
            # a) noisy data   
            exosim_msg ("applying mask and extracting 1D spectrum from noisy data", self.opt.diagnostics)
            self.extractSpec.applyMask_extract_1D()
           
            # b) noiseless data if selected            
            if self.opt.pipeline.useSignal.val == 1:
                exosim_msg ("applying mask and extracting 1D spectrum from signal only data", self.opt.diagnostics)
                self.extractSpec_signal.applyMask_extract_1D()
           
            # c) sample data to find n_pix per bin
            exosim_msg ("applying mask and extracting 1D spectrum to find n_pix per bin" , self.opt.diagnostics)                
            exosim_msg ("extracting 1 D spectra for pixel in bin test",  self.opt.diagnostics)
            self.extractSpec_nPix.applyMask_extract_1D() 
            
            # self.low_SNR_flag = self.extractSpec.low_SNR_flag


        #==============================================================================
        # 3) extract 1D spectrum only if no mask selected 
        #==============================================================================         
        else:
            
            # a) noisy data 
            exosim_msg ("NOT applying mask and extracting 1D spectrum from noisy data", self.opt.diagnostics)
            self.extractSpec.extract1DSpectra()
            
            # b) noiseless data if selected
            if self.opt.pipeline.useSignal.val == 1:  
                exosim_msg ("NOT applying mask and extracting 1D spectrum from signal only data", self.opt.diagnostics)
                self.extractSpec_signal.extract1DSpectra()            
            
            # c) sample data to find n_pix per bin  
            exosim_msg ("NOT applying mask and extracting 1D spectrum to find n_pix per bin", self.opt.diagnostics)
            self.extractSpec_nPix.extract1DSpectra()          
         
            

        self.nPix_1 =  self.extractSpec_nPix.spectra[0] 
        exosim_plot('n_pix per pixel column', self.opt.diagnostics, 
                     xdata=self.opt.cr_wl.value, ydata=self.nPix_1, marker='bo')

        self.OneDSpectra = self.extractSpec.spectra
        
        
        if self.opt.pipeline.pipeline_aperture_corr ==1:
          
            if self.opt.useSignal == 1: # using noiseless signal to get ap corr
                self.extractSpec_no_mask = binning.extractSpec(self.signal, self.metadata)               
                spec_mean = self.extractSpec_signal.spectra.mean(axis=0)
            else: #uses noisy signal to get ap corr - can use - as good and more realistic
                self.extractSpec_no_mask = binning.extractSpec(self.data, self.metadata) 
                spec_mean = self.extractSpec.spectra.mean(axis=0)
            
            self.extractSpec_no_mask.extract1DSpectra()   
            self.OneDSpectra_no_mask = self.extractSpec_no_mask.spectra
            
            spec_mean_no_mask = self.OneDSpectra_no_mask.mean(axis=0)
            ap_corr = spec_mean_no_mask / spec_mean
            
            # plt.figure(11111)
            # plt.plot(self.opt.x_wav_osr[1::3], spec_mean)
            # plt.plot(self.opt.x_wav_osr[1::3], spec_mean_no_mask)
            
            self.extractSpec.spectra  *= ap_corr
            if self.opt.useSignal == 1:  
                self.extractSpec_signal.spectra *= ap_corr
            
            # if self.opt.diagnostics ==1:
            #     plt.figure('check ap corr')
            #     plt.plot(self.opt.x_wav_osr[1::3], self.extractSpec.spectra.mean(axis=0), 'x')
                
            self.OneDSpectra = self.extractSpec.spectra
                  
    
        #==============================================================================
        # 4) Now bin into spectral bins
        #==============================================================================
        # a) noisy data 
        exosim_msg ("binning 1D spectra into spectral bins... from noisy data", self.opt.diagnostics)
        self.extractSpec.binSpectra()    

        # b) noiseless data if selected     
        if self.opt.pipeline.useSignal.val == 1: 
            exosim_msg ("binning 1D spectra into spectral bins... from signal only data", self.opt.diagnostics)
            self.extractSpec_signal.binSpectra() 
            
        # c) sample data to find n_pix per bin    
        exosim_msg ("binning 1D spectra into spectral bins... to find n_pix per bin", self.opt.diagnostics)
        self.extractSpec_nPix.binSpectra() 

        #==============================================================================
        # 5) Define objects from binning process
        #==============================================================================              
        self.binnedLC =  self.extractSpec.binnedLC  
        self.binnedWav =  self.extractSpec.binnedWav 
        
        if self.opt.pipeline.useSignal.val == 1:           
            self.binnedLC_signal =  self.extractSpec_signal.binnedLC    
            self.binnedWav_signal =  self.extractSpec_signal.binnedWav  
        else:
            self.binnedLC_signal =  self.binnedLC 
            self.binnedWav_signal = self.binnedWav        

        if self.opt.timeline.apply_lc.val ==1 :
            self.extractSpec.binGamma()  
            self.binnedGamma =  self.extractSpec.binnedGamma          
         
        self.nPix_2 =  self.extractSpec_nPix.binnedLC[0]
        
        exosim_plot('n_pix per bin', self.opt.diagnostics,
                     xdata=self.binnedWav, ydata =self.nPix_2, marker= 'ro')
  
    

    def extractSpecHotPix(self):    
        self.metadata['diff'] = 1
        bad_pix = np.where(self.dq_array==0,0,1) # makes all pixels that are hot or saturated etc = 1
        self.extractSpecHotPix = binning.extractSpec(bad_pix, self.metadata)           
#        self.extractSpecHotPix.applyMask()               
#        self.extractSpecHotPix.extract1DSpectra()      
        self.extractSpecHotPix.applyMask_extract_1D()                 
        self.extractSpecHotPix.binSpectra()              
        self.binnedHotPix =  self.extractSpecHotPix.binnedLC  
        
        self.extractSpecAllPix = binning.extractSpec(self.dq_array*0+1, self.metadata)           
#        self.extractSpecAllPix.applyMask()               
#        self.extractSpecAllPix.extract1DSpectra()      
        self.extractSpecAllPix.applyMask_extract_1D()                 
        self.extractSpecAllPix.binSpectra()              
        self.binnedAllPix =  self.extractSpecAllPix.binnedLC              
              
        if self.opt.ch == 'AIRS CH1':
            self.binnedHotPix = self.binnedHotPix[:,:-1]
            self.binnedAllPix = self.binnedAllPix[:,:-1] 
        if self.opt.ch == 'AIRS CH0':
            self.binnedHotPix = self.binnedHotPix[:,3:-12]
            self.binnedAllPix = self.binnedAllPix[:,3:-12]
        if self.opt.ch == 'NIR Spec':
            self.binnedHotPix = self.binnedHotPix[:,3:]
            self.binnedAllPix = self.binnedAllPix[:,3:]

         
    def extractPhot_decorr(self):
        self.extractPhot_decorr = binning.extractPhot_decorr(self.data, self.metadata) 
        if self.opt.useSignal == 1:  
            self.extractPhot_decorr_signal = binning.extractPhot_decorr(self.signal, self.metadata)         
        
        self.binnedLC =  self.extractPhot_decorr.photLC    
        if self.opt.useSignal == 1: 
            self.binnedLC_signal =  self.extractPhot_decorr_signal.photLC  
        else:
            self.binnedLC_signal =  self.binnedLC 
   
        self.low_SNR_flag = self.extractPhot_decorr.low_SNR_flag
       
#==============================================================================
#         Aperture correction
#==============================================================================
           
        if self.opt.aperture_corr ==1:   
            s= []   
            if self.opt.useSignal == 1: # using noiseless signal to get ap corr
                data0 = self.signal # maskless signal
                s_mask = self.binnedLC_signal
               
            else:# using noisy signal to get ap corr
                data0 = self.data # maskless signal 
                s_mask = self.binnedLC
              
                
            for i in range(data0.shape[2]):
                
                s.append(data0[...,i].sum())
                                 
            s = np.array(s)
       
            # plt.figure('ap corr test')
            # plt.plot(s, 'rx')
            # plt.plot(s_mask, 'b+')
            # plt.ylim (np.min(s_mask)- np.min(s_mask)*0.1 ,np.max(s)+np.max(s)*0.1)
            
            ap_corr = s.mean() / s_mask.mean()
            # print ("aperture correction factor", ap_corr)
            
            # plt.figure('ap corr test 2')
            # plt.plot(s.mean(), 'ro')
            # plt.plot(s_mask.mean(), 'bo')
            # plt.plot(s_mask.mean()*ap_corr, 'gx')
            
     
            self.binnedLC *= ap_corr
            if self.opt.useSignal == 1:  
                self.binnedLC_signal *= ap_corr
 

#==============================================================================
#         
#==============================================================================

        self.binnedWav =  self.extractPhot_decorr.photWav
            
        self.binnedGamma =  self.extractPhot_decorr.photGamma
        self.OneDSpectra = self.binnedLC 
        
        d_x_wav_osr = np.abs(np.gradient(self.opt.x_wav_osr))
        self.binnedWavEdges = [self.opt.x_wav_osr[0]-d_x_wav_osr[0] , self.opt.x_wav_osr[-1]+ d_x_wav_osr[-1]]
        self.binnedWavBinSizes = (self.opt.x_wav_osr[-1]+ d_x_wav_osr[-1]) - (self.opt.x_wav_osr[0]-d_x_wav_osr[0])
        
        
    def extractPhot_fixed(self):
        self.extractPhot_fixed = binning.extractPhot_fixed(self.data, self.metadata) 
        self.binnedLC =  self.extractPhot_fixed.photLC    
        self.binnedWav =  self.extractPhot_fixed.photWav   
        

    def extractPhot_no_ap(self):
        self.extractPhot_no_ap = binning.extractSpec(self.data, self.metadata) 
        self.extractPhot_no_ap.extract1DSpectra()
        self.binnedLC =   self.extractPhot_no_ap.spectra.sum(axis=1)
        
        self.binnedWav =  self.metadata['WL'][0]
        self.OneDSpectra = self.binnedLC
        

    def findSlitloss(self):
        self.Slitloss = binning.findSlitloss(self.data, self.metadata) 
        self.OneD_PSFs_unclipped = self.Slitloss.OneD_PSFs_unclipped
        self.OneD_PSFs_clipped = self.Slitloss.OneD_PSFs_clipped

        

class pipeline_stage_2():
       
    def __init__(self, opt):
        
        self.opt  = opt
        self.binnedWav = opt.pipeline_stage_1.binnedWav
        self.binnedLC = opt.pipeline_stage_1.binnedLC
        self.binnedLC_signal = opt.pipeline_stage_1.binnedLC_signal
        if self.opt.timeline.apply_lc.val ==1 :
            self.binnedGamma = opt.pipeline_stage_1.binnedGamma
        self.nPix_2 =  opt.pipeline_stage_1.nPix_2
        self.nPix_1 =  opt.pipeline_stage_1.nPix_1
        self.ApFactor=  opt.pipeline_stage_1.ApFactor
        self.exp_end_time_grid = opt.ndr_end_time[opt.effective_multiaccum-1::opt.effective_multiaccum]
        self.ndr_end_time_grid = opt.ndr_end_time
        self.data_raw = opt.data_raw
        
        if self.opt.simulation.sim_use_systematic_model.val ==1 :
            self.detrend()  
        if self.opt.timeline.apply_lc.val ==1 :
            self.fitLC()            
        elif self.opt.timeline.apply_lc.val ==0 :
            self.ootSNR()
   
    def detrend(self):
        self.detrend = detrend.Detrend(self.binnedLC, self.binnedWav, self.opt)   
      
        
    def fitLC(self):
        self.processLC = binning.processLC(self.binnedLC, self.binnedWav, self.binnedGamma, self.opt)   
        self.transitDepths = self.processLC.transitDepths
        self.model_gamma1 = self.processLC.model_gamma1
        self.model_gamma2 = self.processLC.model_gamma2
        self.model_f = self.processLC.model_f
        self.binnedWav = self.processLC.binnedWav

    def ootSNR(self):
        if self.opt.pipeline.useSignal.val == 1:
            self.processOOT = binning.processOOT(self.binnedLC, self.binnedLC_signal,self.binnedWav, self.opt)   
        else:
            self.processOOT = binning.processOOT(self.binnedLC, self.binnedLC, self.binnedWav, self.opt)   
        
        self.ootSignal = self.processOOT.ootSignal     
        self.ootNoise = self.processOOT.ootNoise
             
        if self.opt.pipeline.useAllen.val ==1:
            self.ootAllen = self.processOOT.ootAllen
            self.noiseAt1hr = self.processOOT.noiseAt1hr
          
