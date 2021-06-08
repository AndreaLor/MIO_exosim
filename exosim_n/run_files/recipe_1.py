'''
exosim_n 
2.0
Recipe 1 - OOT simulation returning stellar signal and noise per spectral bin with Allan analysis

'''
from exosim_n.modules import astroscene, telescope, channel, backgrounds
from exosim_n.modules import detector, timeline, light_curve, systematics, noise, output 
from exosim_n.pipeline.run_pipeline import pipeline_stage_1, pipeline_stage_2
from exosim_n.lib.exosim_n_lib import exosim_n_msg, exosim_n_plot, write_record
from astropy import units as u
from datetime import datetime
import pickle, os
import numpy as np
from exosim_n.classes.options import Options
class recipe_1(object):
    def __init__(self, opt):
        DEBUG=Options.DEBUG
        output_directory = opt.common.output_directory.val
        filename =""
        
        self.results_dict ={}
        self.results_dict['simulation_mode'] = opt.simulation.sim_mode.val
        self.results_dict['simulation_realisations'] = opt.simulation.sim_realisations.val
        self.results_dict['ch'] =  opt.observation.obs_channel.val 
        self.noise_dict ={}    
   
        opt.pipeline.useSignal.val=1
        opt.pipeline.split  = 0
        opt.noise.ApplyRandomPRNU.val=1
                      
        opt.timeline.apply_lc.val = 0
        opt.timeline.useLDC.val = 0
        opt.pipeline.useAllen.val =1
        opt.timeline.use_T14.val=0
        opt.timeline.obs_time.val = 0*u.hr # 0 means n_exp overides obs_time
        opt.timeline.n_exp.val = 1000 # uses 1000 exposures 
        # opt.timeline.n_exp.val = 50
     
        noise_type = int(opt.noise.sim_noise_source.val)
   
               
        nb_dict = {'rn'           :[1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                   'sn'           :[1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
                   'spat'         :[1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0],                                   
                   'spec'         :[1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0],                                      
                   'emm_switch'   :[1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1],                    
                   'zodi_switch'  :[1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1],    
                   'dc_switch'    :[1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                   'source_switch':[1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
                   'diff'         :[0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
                   'jitter_switch':[1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0],
                   'noise_tag': [ 'All noise','All photon noise','Source photon noise','Dark current noise',
                        'Zodi noise','Emission noise','Read noise','Spatial jitter noise',
                        'Spectral jitter noise','Combined jitter noise','No noise - no background','No noise - all background', 'Fano noise'],  
                        'color': ['0.5','b', 'b','k','orange','pink', 'y','g','purple','r', '0.8', 'c']
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
                      
        self.noise_dict[nb_dict['noise_tag'][noise_type]] ={}
        
        exosim_n_msg ("Noise type: %s"%(nb_dict['noise_tag'][noise_type]), 1) 
  
        start = 0 
        end = int(start + opt.no_real)

        for j in range(start, end):
            
            
            if (opt.no_real-start) > 1:
                   exosim_n_msg ("", 1)    
                   exosim_n_msg ("============= REALIZATION %s ============="%(j), 1)
                   exosim_n_msg (opt.lab, 1)                   
                   exosim_n_msg ("", 1)  
            
            opt = self.run_exosim_nA(opt)
            opt = self.run_exosim_nA1(opt)
            
            if opt.observation_feasibility ==0:      
                exosim_n_msg ("Observation not feasible...", opt.diagnostics) 
                self.feasibility = 0
            else:
                self.feasibility = 1
                              
                opt = self.run_exosim_nB(opt)
                
                if opt.simulation.sim_output_type.val == 1:  
                    opt = self.run_pipeline_stage_1(opt)       
                    opt = self.run_pipeline_stage_2(opt)              
                    self.pipeline = opt.pipeline_stage_2
                    
                    
                    self.noise_dict[opt.noise_tag]['wl'] = self.pipeline.binnedWav
                    self.results_dict['input_spec'] = opt.cr
                    self.results_dict['input_spec_wl'] = opt.cr_wl            
                    
                             
                    if j == start:                    
                        self.noise_dict[opt.noise_tag]['signal_std_stack'] = self.pipeline.ootNoise
                        self.noise_dict[opt.noise_tag]['signal_mean_stack'] = self.pipeline.ootSignal
                        if opt.pipeline.useAllen.val == 1:
                            self.noise_dict[opt.noise_tag]['fracNoT14_stack'] = self.pipeline.noiseAt1hr
                            
                        self.noise_dict[opt.noise_tag]['signal_std_mean'] = self.pipeline.ootNoise
                        self.noise_dict[opt.noise_tag]['signal_mean_mean'] = self.pipeline.ootSignal
                        if opt.pipeline.useAllen.val == 1:
                            self.noise_dict[opt.noise_tag]['fracNoT14_mean'] = self.pipeline.noiseAt1hr
                            
                        self.noise_dict[opt.noise_tag]['signal_std_std'] = np.zeros(len(self.pipeline.binnedWav))
                        self.noise_dict[opt.noise_tag]['signal_mean_std'] = np.zeros(len(self.pipeline.binnedWav))
                        if opt.pipeline.useAllen.val == 1:
                            self.noise_dict[opt.noise_tag]['fracNoT14_std'] = np.zeros(len(self.pipeline.binnedWav))
                       

                    else:
                        self.noise_dict[opt.noise_tag]['signal_std_stack'] = np.vstack((self.noise_dict[opt.noise_tag]['signal_std_stack'], self.pipeline.ootNoise))
                        self.noise_dict[opt.noise_tag]['signal_mean_stack'] = np.vstack((self.noise_dict[opt.noise_tag]['signal_mean_stack'], self.pipeline.ootSignal))
                        if opt.pipeline.useAllen.val == 1:
                            self.noise_dict[opt.noise_tag]['fracNoT14_stack'] = np.vstack((self.noise_dict[opt.noise_tag]['fracNoT14_stack'], self.pipeline.noiseAt1hr))
    
                        self.noise_dict[opt.noise_tag]['signal_std_mean'] = self.noise_dict[opt.noise_tag]['signal_std_stack'].mean(axis=0)
                        self.noise_dict[opt.noise_tag]['signal_mean_mean'] = self.noise_dict[opt.noise_tag]['signal_mean_stack'].mean(axis=0)
                        if opt.pipeline.useAllen.val == 1:
                            self.noise_dict[opt.noise_tag]['fracNoT14_mean'] = self.noise_dict[opt.noise_tag]['fracNoT14_stack'].mean(axis=0)
                            
                        self.noise_dict[opt.noise_tag]['signal_std_std'] = self.noise_dict[opt.noise_tag]['signal_std_stack'].std(axis=0)
                        self.noise_dict[opt.noise_tag]['signal_mean_std'] = self.noise_dict[opt.noise_tag]['signal_mean_stack'].std(axis=0)
                        if opt.pipeline.useAllen.val == 1:
                            self.noise_dict[opt.noise_tag]['fracNoT14_std'] = self.noise_dict[opt.noise_tag]['fracNoT14_stack'].std(axis=0)
    
                    self.noise_dict[opt.noise_tag]['bad_map'] = opt.bad_map
                    self.noise_dict[opt.noise_tag]['example_exposure_image'] = opt.exp_image
                    self.noise_dict[opt.noise_tag]['pixel wavelengths'] = opt.x_wav_osr[1::3].value

                    self.results_dict['noise_dic'] = self.noise_dict
    
                    time_tag = (datetime.now().strftime('%Y_%m_%d_%H%M_%S'))
                    self.results_dict['time_tag'] =  time_tag
             
                    if j == start:
                        filename = '%s/OOT_SNR_%s_%s_TEMP.pickle'%(output_directory, opt.lab, time_tag)
                        self.filename = 'OOT_SNR_%s_%s_TEMP.pickle'%(opt.lab, time_tag)
                        with open(filename, 'wb') as handle:
                            pickle.dump(self.results_dict , handle, protocol=pickle.HIGHEST_PROTOCOL)
                        write_record(opt, output_directory, self.filename, opt.params_file_path)
                        
                    if j != start:
                        os.remove(filename)  # delete previous temp file  
                        txt_file = '%s.txt'%(filename)
                        os.remove(txt_file)         
                        filename = '%s/OOT_SNR_%s_%s_TEMP.pickle'%(output_directory, opt.lab, time_tag)
                        self.filename = 'OOT_SNR_%s_%s_TEMP.pickle'%(opt.lab, time_tag)
                        with open(filename, 'wb') as handle:
                            pickle.dump(self.results_dict , handle, protocol=pickle.HIGHEST_PROTOCOL)    
                        write_record(opt, output_directory, self.filename, opt.params_file_path)
                       
                    if j == end-1:
                        
                        os.remove(filename)  # delete previous temp file
                        txt_file = '%s.txt'%(filename)
                        os.remove(txt_file)             
                        filename = '%s/OOT_SNR_%s_%s.pickle'%(output_directory, opt.lab, time_tag)
                        filename = '%s/FIXED_NAME.pickle'%(output_directory) 
                        
                        with open(filename, 'wb') as handle:
                            pickle.dump(self.results_dict , handle, protocol=pickle.HIGHEST_PROTOCOL)
                        
                        exosim_n_msg('Results in %s'%(filename), 1)
                        #self.filename = 'OOT_SNR_%s_%s.pickle'%(opt.lab, time_tag)
                        self.filename = 'FIXED_NAME.pickle'
                        write_record(opt, output_directory, self.filename, opt.params_file_path)
                                
                     
                elif opt.simulation.sim_output_type.val == 2:
                     
                    output.run(opt)
                    
                             
                            
    def run_exosim_nA(self, opt):
      exosim_n_msg('Exosystem', 1)
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
  

  
    
  
    
  
