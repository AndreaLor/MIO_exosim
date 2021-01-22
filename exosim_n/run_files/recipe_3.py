'''
exosim_n

Recipe 3 - OOT simulation returning a noise budget
 
'''
import numpy as np
import time, os, pickle
from datetime import datetime
from exosim_n.modules import astroscene, telescope, channel, backgrounds
from exosim_n.modules import detector, timeline, light_curve, systematics, noise
from exosim_n.pipeline.run_pipeline import pipeline_stage_1, pipeline_stage_2
from exosim_n.lib.exosim_n_lib import exosim_n_msg, exosim_n_plot, write_record
from astropy import units as u

 
class recipe_3(object):
    def __init__(self, opt):
        
        output_directory = opt.common.output_directory.val
        filename=""
        
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
        opt.pipeline.useAllen.val =0
        
        opt.pipeline.pipeline_auto_ap.val = 0 # for noise budget keep this to a fixed value (i.e. choose 0) so same for all sims
         
        opt.timeline.obs_time.val = 0*u.hr
        opt.timeline.n_exp.val = 1000.0
         
        noise_list = [0,2,3,4,5,6,7,8,9]
        
 
        start = 0 
        end = int(start + opt.no_real)       
                
        nb_dict = {'rn'           :[1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                   'sn'           :[1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
                   'spat'         :[1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0],                                   
                   'spec'         :[1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0],                                      
                   'emm_switch'   :[1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1],                           
                   'zodi_switch'  :[1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1],       
                   'dc_switch'    :[1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1],
                   'source_switch':[1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1],
                   'diff'         :[0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0],
                   'jitter_switch':[1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0],
                   'noise_tag': [ 'All noise','All photon noise','Source photon noise','Dark current noise',
                        'Zodi noise','Emission noise','Read noise','Spatial jitter noise',
                        'Spectral jitter noise','Combined jitter noise','No noise - no background','No noise - all background'],  
                        'color': ['0.5','b', 'b','k','orange','pink', 'y','g','purple','r', '0.8','c']
              } 
        
        for i in noise_list:
            self.noise_dict[nb_dict['noise_tag'][i]] ={}
            
        for j in range(start,end):
                
            seed = np.random.randint(100000000)
            for i in noise_list:
                                      
                opt.noise.EnableReadoutNoise.val = nb_dict['rn'][i]
                opt.noise.EnableShotNoise.val = nb_dict['sn'][i]
                opt.noise.EnableSpatialJitter.val= nb_dict['spat'][i]
                opt.noise.EnableSpectralJitter.val= nb_dict['spec'][i]   
                opt.background.EnableEmission.val = nb_dict['emm_switch'][i]
                opt.background.EnableZodi.val = nb_dict['zodi_switch'][i]    
                opt.background.EnableDC.val  =  nb_dict['dc_switch'][i]
                opt.background.EnableSource.val  = nb_dict['source_switch'][i]
                opt.diff = nb_dict['diff'][i]      
                opt.noise_tag = nb_dict['noise_tag'][i]
                opt.color = nb_dict['color'][i] 
                
                exosim_n_msg('==========================================', 1)
                exosim_n_msg('Noise source:%s'%(opt.noise_tag), 1)
                
                if (opt.no_real-start) > 1:      
                   print ("============= REALIZATION %s ============="%(j))
                   print (opt.lab)               
                
                np.random.seed(seed)
                opt = self.run_exosim_n_A(opt)
                opt = self.run_exosim_n_A1(opt) 
                
                opt.observation_feasibility =1  
                    
                if opt.observation_feasibility ==0:      
                    exosim_n_msg ("Observation not feasible...", opt.diagnostics)
                    self.feasibility = 0    
                else:
                    self.feasibility = 1    
                    opt = self.run_exosim_n_B(opt)
                    
                    if opt.simulation.sim_output_type.val == 1:  
                    
                        opt = self.run_pipeline_stage_1(opt) 
                        opt = self.run_pipeline_stage_2(opt) 
                        
                        self.pipeline = opt.pipeline_stage_2
                        self.opt = opt
                        
                        self.noise_dict[opt.noise_tag]['wl'] = self.pipeline.binnedWav
                         
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

                            self.noise_dict[opt.noise_tag]['bad_map'] = opt.bad_map
                            self.noise_dict[opt.noise_tag]['example_exposure_image'] = opt.exp_image
                            self.noise_dict[opt.noise_tag]['pixel wavelengths'] = opt.x_wav_osr[1::3].value
                                             

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
                        
                    elif opt.simulation.sim_output_type.val == 2:
                        
                        filename = output.run(opt)
                        write_record(opt, output_directory, filename, opt.params_file_path)

            if self.feasibility ==1:      

                # dump pickle file at end of each cycle of noise sims
                if opt.simulation.sim_output_type.val == 1:         
                    time_tag = (datetime.now().strftime('%Y_%m_%d_%H%M_%S'))
                    self.results_dict['time_tag'] =  time_tag
             
                    # if j != start:
                    #     os.remove(filename)  # delete previous temp file
                    # filename = '%s/Noise_budget_%s_TEMP.pickle'%(output_directory, opt.lab)
                    # with open(filename, 'wb') as handle:
                    #     pickle.dump(self.results_dict , handle, protocol=pickle.HIGHEST_PROTOCOL)
                       
                    if j == end-1:
                        # os.remove(filename)  # delete previous temp file
                        filename = '%s/Noise_budget_%s_%s.pickle'%(output_directory, opt.lab, time_tag)
                        with open(filename, 'wb') as handle:
                            pickle.dump(self.results_dict , handle, protocol=pickle.HIGHEST_PROTOCOL)
         
                            exosim_n_msg('Results in %s'%(filename), 1)
                            self.filename = 'Noise_budget_%s_%s.pickle'%(opt.lab, time_tag)

                            write_record(opt, output_directory, self.filename, opt.params_file_path)
    def run_exosim_n_A(self, opt):
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
      exosim_n_msg('Timeline', 1)
      timeline.run(opt)
      exosim_n_msg('Light curve', 1)
      light_curve.run(opt)        
      return opt

    def run_exosim_n_A1(self, opt):
      exosim_n_msg('Systematics', 1)
      systematics.run(opt)               
      return opt     
      
    def run_exosim_n_B(self, opt):
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
  

