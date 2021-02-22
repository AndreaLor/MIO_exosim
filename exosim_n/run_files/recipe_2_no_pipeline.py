'''
exosim_n
2.0
Recipe 2 :
Full transit simulation

This recipe only generates raw data to be used in a subsequent pipeline.

'''

import numpy as np
import time, os, pickle
from datetime import datetime
from exosim_n.modules import astroscene, telescope, channel, backgrounds, output
from exosim_n.modules import detector, timeline, light_curve, systematics, noise
from exosim_n.lib.exosim_n_lib import exosim_n_msg, exosim_n_plot
from exosim_n.lib.output_lib import write_record_no_pipeline


class recipe_2_no_pipeline(object):
    def __init__(self, opt):
        
        output_directory = opt.common.output_directory.val
        filename=""
        
        self.results_dict ={}
        self.results_dict['simulation_mode'] = opt.simulation.sim_mode.val
        self.results_dict['simulation_realisations'] = opt.simulation.sim_realisations.val
        self.results_dict['ch'] =  opt.observation.obs_channel.val 

        opt.pipeline.useSignal.val=0

        opt.noise.ApplyRandomPRNU.val=1

        opt.timeline.apply_lc.val = 1
        opt.timeline.useLDC.val = 1

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
         
        opt.observation_feasibility = 1   
                          
        if opt.observation_feasibility ==0:      
           exosim_n_msg ("Observation not feasible...", opt.diagnostics) 
           self.feasibility = 0
        else:
           self.feasibility = 1
           
        for j in range(start, end):
            
            if (opt.no_real-start) > 1:
                exosim_n_msg ("", 1)    
                exosim_n_msg ("============= REALIZATION %s ============="%(j), 1)
                exosim_n_msg (opt.lab, 1)                   
                exosim_n_msg ("", 1) 
            
            pp = time.time()
            
            opt = self.run_exosim_nA1(opt)  # set QE grid for this realization
            exosim_n_msg ("QE variations set", 1) 
            exosim_n_msg ("Number of exposures %s"%(opt.n_exp), 1) 

            opt = self.run_exosim_nB(opt)
            
            filename = output.run(opt)
        
            write_record_no_pipeline(opt, output_directory, filename, opt.params_file_path)
            exosim_n_msg('File saved as %s/%s'%(output_directory, filename), 1)

         
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
    
