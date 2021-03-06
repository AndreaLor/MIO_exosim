'''
exosim_n
2.0
Recipe 2 :
Monte Carlo full transit simulation returning transit depths and noise on transit dept per spectral bin

Each realisation run the stage 2 exosim_n routine with a new noise randomization
The if the pipeline_auto_ap=1, the aperture mask size may vary with each realisation
This should not matter since final measurement is the fractional transit depth

'''

import numpy as np
import time, os, pickle
from datetime import datetime
from exosim_n.modules import astroscene, telescope, channel, backgrounds
from exosim_n.modules import detector, timeline, light_curve, systematics, noise
from exosim_n.pipeline.run_pipeline import pipeline_stage_1, pipeline_stage_2
from exosim_n.lib.exosim_n_lib import exosim_n_msg, exosim_n_plot, write_record
from exosim_n.classes.options import Options

import matplotlib.pyplot as plt
class recipe_2(object):
    def __init__(self, opt):
        DEBUG= Options.DEBUG
        if DEBUG:
            oldDiagnostic=opt.diagnostics
            opt.diagnostics=Options.showPipeline1
        output_directory = opt.common.output_directory.val
        filename=""
        
        self.results_dict ={}
        self.results_dict['simulation_mode'] = opt.simulation.sim_mode.val
        self.results_dict['simulation_realisations'] = opt.simulation.sim_realisations.val
        self.results_dict['ch'] =  opt.observation.obs_channel.val 

        opt.pipeline.useSignal.val=0

        opt.pipeline.split  = 0
        opt.noise.ApplyRandomPRNU.val=1

        opt.timeline.apply_lc.val = 1
        opt.timeline.useLDC.val = 1
        opt.pipeline.useAllen.val =0
        opt.pipeline.fit_gamma.val  =0 #keep zero for uncert on p

        noise_type = int(opt.noise.sim_noise_source.val)
   
        nb_dict = {'rn'            :[1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                'sn'               :[1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
                'spat'             :[1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0],                                   
                'spec'             :[1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0],                                      
                'emm_switch'       :[1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1],                           
                'zodi_switch'      :[1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1],       
                'dc_switch'        :[1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1],
                'source_switch'    :[1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1],
                'diff'             :[0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0],
                'jitter_switch'    :[1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0],
                'noise_tag': [ 'All noise','All photon noise','Source photon noise','Dark current noise',
                        'Zodi noise','Emission noise','Read noise','Spatial jitter noise',
                        'Spectral jitter noise','Combined jitter noise','No noise - no background','No noise - all background'],  
                        'color': ['0.5','b', 'b','k','orange','pink', 'y','g','purple','r', '0.8','c']
              }         

        opt.noise.EnableReadoutNoise.val = nb_dict['rn'][noise_type]
        opt.noise.EnableShotNoise.val = nb_dict['sn'][noise_type]
        opt.noise.EnableSpatialJitter.val= nb_dict['spat'][noise_type]
        opt.noise.EnableSpectralJitter.val= nb_dict['spec'][noise_type]   
        opt.background.EnableEmission.val = nb_dict['emm_switch'][noise_type]
        opt.background.EnableZodi.val = nb_dict['zodi_switch'][noise_type]    
        opt.background.EnableDC.val  =  nb_dict['dc_switch'][noise_type]
        opt.background.EnableSource.val  = nb_dict['source_switch'][noise_type]
        opt.diff = nb_dict['diff'][noise_type]      
        opt.noise_tag = nb_dict['noise_tag'][noise_type]
        opt.color = nb_dict['color'][noise_type]

        start = 0 
        end = int(start + opt.no_real)
        
        if (opt.no_real-start) > 1:
            exosim_n_msg ("Monte Carlo selected", 1) 
                           
        opt = self.run_exosim_nA(opt)
        
                                   
        if opt.observation_feasibility ==0:      
           exosim_n_msg ("Observation not feasible...", opt.diagnostics) 
           self.feasibility = 0
        else:
           self.feasibility = 1
           n_ndr0 = opt.n_ndr*1
           
           ndr_end_frame_number0 = opt.ndr_end_frame_number*1
           frames_per_ndr0 = opt.frames_per_ndr*1
           duration_per_ndr0 = opt.duration_per_ndr*1
           n_exp0 = opt.n_exp
           lc0 = opt.lc_original*1   # important this happens at this stage      
                        
           if n_ndr0 > 10000:

               opt.pipeline.split = 1
               if opt.diagnostics ==1 :
                   exosim_n_msg ('number of NDRs > 10000: using split protocol', opt.diagnostics)
           else:
               opt.pipeline.split = 0 
               
           # #delete    
           # opt.pipeline.split = 1 
           
           for j in range(start, end):
               
               if (opt.no_real-start) > 1:
                   exosim_n_msg ("", 1)    
                   exosim_n_msg ("============= REALIZATION %s ============="%(j), 1)
                   exosim_n_msg (opt.lab, 1)                   
                   exosim_n_msg ("", 1) 
               
               pp = time.time()
               
               opt = self.run_exosim_nA1(opt)  # set QE grid for this realization
               exosim_n_msg ("QE variations set", 1) 
               exosim_n_msg ("Number of exposures %s"%(n_exp0), 1) 
               
               
               
# =============================================================================
#  # split simulation into chunks to permit computation - makes no difference to final results    
# =============================================================================
               if opt.pipeline.split ==1:
                  
                   exosim_n_msg('Splitting data series into chunks', opt.diagnostics)
                   # uses same QE grid and jitter timeline but otherwise randomoses noise
                   ndrs_per_round = opt.effective_multiaccum*int(5000/opt.multiaccum)  
                   # ndrs_per_round = opt.effective_multiaccum*int(500/opt.multiaccum)  
  
                   total_chunks = len(np.arange(0, n_ndr0, ndrs_per_round))
                   
                   idx = np.arange(0, n_ndr0, ndrs_per_round) # list of starting ndrs
                   
                   for i in range(len(idx)):
                       exosim_n_msg('=== Realisation %s Chunk %s / %s====='%(j, i+1, total_chunks), opt.diagnostics)
                       
                       if idx[i] == idx[-1]:
                           opt.n_ndr = n_ndr0 - idx[i]
                           opt.lc_original = lc0[:,idx[i]:]
                           opt.ndr_end_frame_number = ndr_end_frame_number0[idx[i]:]
                           opt.frames_per_ndr=  frames_per_ndr0[idx[i]:]
                           opt.duration_per_ndr = duration_per_ndr0[idx[i]:]  
                          
                       else:
                           opt.n_ndr = idx[i+1]- idx[i]
                           opt.lc_original = lc0[:,idx[i]:idx[i+1]]
                           print ('idx start......', idx[i])
                           
                           opt.ndr_end_frame_number = ndr_end_frame_number0[idx[i]: idx[i+1]]
                           opt.frames_per_ndr=  frames_per_ndr0[idx[i]: idx[i+1]]
                           opt.duration_per_ndr = duration_per_ndr0[idx[i]: idx[i+1]]
                  
                       opt.n_exp = int(opt.n_ndr/opt.effective_multiaccum)
                           
                       if i == 0:
                           opt.pipeline.pipeline_auto_ap.val= 1
                           opt.use_external_jitter = 0
                           
                           opt = self.run_exosim_nB(opt)
                           opt  = self.run_pipeline_stage_1(opt)  
                           
                           opt.pipeline.pipeline_ap_factor.val= opt.AvBest 
                           if (opt.noise.EnableSpatialJitter.val  ==1 or opt.noise.EnableSpectralJitter.val  ==1 or opt.noise.EnableAll.val == 1) and opt.noise.DisableAll.val != 1:
                               opt.input_yaw_jitter, opt.input_pitch_jitter, opt._input_frame_osf = opt.yaw_jitter, opt.pitch_jitter, opt.frame_osf                     
                                                      
                       else:
                           opt.pipeline.pipeline_auto_ap.val = 0
                           opt.use_external_jitter = 1 # uses the jitter timeline from the first realization

                           opt = self.run_exosim_nB(opt)
                           opt  = self.run_pipeline_stage_1(opt)
                       
                                  
                       exosim_n_msg('Aperture used %s'%(opt.pipeline.pipeline_ap_factor.val), opt.diagnostics)
                       binnedLC = opt.pipeline_stage_1.binnedLC
                       data = opt.pipeline_stage_1.opt.data_raw
                                                                  
                       if i ==0:
                           data_stack = data
                           binnedLC_stack = binnedLC    
                       else:
                           data_stack = np.dstack((data_stack,data))
                           binnedLC_stack = np.vstack((binnedLC_stack,binnedLC))                  

                   aa = data_stack.sum(axis=0)
                   bb = aa.sum(axis=0)
                   exosim_n_plot('test_from_sim', opt.diagnostics,
                            ydata=bb[opt.effective_multiaccum::opt.effective_multiaccum] )
                   aa = binnedLC_stack.sum(axis=1)
                  
                   exosim_n_plot('test_from_pipeline', opt.diagnostics,
                                        ydata=aa)                            

                   opt.n_ndr  = n_ndr0             
                   opt.ndr_end_frame_number  = ndr_end_frame_number0  
                   opt.frames_per_ndr  = frames_per_ndr0 
                   opt.duration_per_ndr = duration_per_ndr0
                   opt.n_exp = n_exp0                                
                           
               elif opt.pipeline.split ==0:

                   opt  = self.run_exosim_nB(opt)
                   if j==start:  # first realization sets the ap, then the other use the same one
                       opt.pipeline.pipeline_auto_ap.val= 1
                   else:
                       opt.pipeline.pipeline_auto_ap.val= 0
                   opt  = self.run_pipeline_stage_1(opt)
                   if j==start:  # first realization sets the ap, then the other use the same one
                       opt.pipeline.pipeline_ap_factor.val= opt.AvBest

        
                   binnedLC_stack  = opt.pipeline_stage_1.binnedLC                   
                   exosim_n_plot('testvvv', opt.diagnostics,
                                ydata=binnedLC_stack.sum(axis=1) )  
                   
               
               if DEBUG:
                         if opt.diagnostics:
                                 plt.show()
                         opt.diagnostics=oldDiagnostic 
               opt.pipeline_stage_1.binnedLC = binnedLC_stack     
               opt = self.run_pipeline_stage_2(opt)
               pipeline = opt.pipeline_stage_2
       
               p = pipeline.transitDepths
               if j == start:
                   p_stack = p
               else:
                   p_stack = np.vstack((p_stack,p))
                                      
               exosim_n_msg ("time to complete realization %s %s"%(j, time.time()-pp ) ,opt.diagnostics)
        

               self.results_dict['wl'] = pipeline.binnedWav   
               self.results_dict['input_spec'] = opt.cr
               self.results_dict['input_spec_wl'] = opt.cr_wl
              
               if j==start:  # if only one realisation slightly different format
                   self.results_dict['p_stack'] = np.array(p)
                   self.results_dict['p_std']= np.zeros(len(p))
                   self.results_dict['p_mean'] = np.array(p)             
               else:
                   self.results_dict['p_stack'] = np.vstack((self.results_dict['p_stack'], p))
                   self.results_dict['p_std'] = self.results_dict['p_stack'].std(axis=0)  
                   self.results_dict['p_mean'] = self.results_dict['p_stack'].mean(axis=0)

                      
               time_tag = (datetime.now().strftime('%Y_%m_%d_%H%M_%S'))
                   
               self.results_dict['time_tag'] =  time_tag
               self.results_dict['bad_map'] = opt.bad_map
               self.results_dict['example_exposure_image'] = opt.exp_image
               self.results_dict['pixel wavelengths'] = opt.x_wav_osr[1::3].value
                                 
    
               if j != start:
                   os.remove(filename)  # delete previous temp file
     
               filename = '%s/Full_transit_%s_TEMP.pickle'%(output_directory, opt.lab)
               with open(filename, 'wb') as handle:
                   pickle.dump(self.results_dict , handle, protocol=pickle.HIGHEST_PROTOCOL)
                   
                
           os.remove(filename)  # delete previous temp file
           # write final file
           filename = '%s/Full_transit_%s_%s.pickle'%(output_directory, opt.lab, time_tag)
           with open(filename, 'wb') as handle:
                  pickle.dump(self.results_dict , handle, protocol=pickle.HIGHEST_PROTOCOL)
        
           exosim_n_msg('Results in %s'%(filename), 1)
           self.filename = 'Full_transit_%s_%s.pickle'%(opt.lab, time_tag)
               
           write_record(opt, output_directory, self.filename, opt.params_file_path)
            
    def run_exosim_nA(self, opt):
      exosim_n_msg('Astroscene', 1)
      astroscene.run(opt)
      exosim_n_msg('Telescope', 1)
      telescope.run(opt)
      exosim_n_msg('Channel', 1)
      channel.run(opt)
      exosim_n_msg('Backgrounds', 1)
      backgrounds.run(opt) 
      exosim_n_msg('Detector', 1)
      detector.run(opt)
      if opt.observation_feasibility ==1: # if detector does not saturate continue
          exosim_n_msg('Timeline', 1)
          timeline.run(opt)
          exosim_n_msg('Light curve', 1)
          light_curve.run(opt)     
          return opt       
      else: # if detector saturates end sim      
          return opt

    def run_exosim_nA1(self, opt):
      exosim_n_msg('Systematics', 1)
      systematics.run(opt)               
      return opt
        
    def run_exosim_nB(self, opt):
      exosim_n_msg('Noise', 1)
      noise.run(opt)                 
      return opt
           
    def run_pipeline_stage_1(self, opt):
      exosim_n_msg('Pipeline stage 1', 1)
      opt.pipeline_stage_1 = pipeline_stage_1(opt)   
      return opt  
             
    def run_pipeline_stage_2(self, opt):    
      exosim_n_msg('Pipeline stage 2', 1)
      opt.pipeline_stage_2 = pipeline_stage_2(opt)             
      return opt 
    
