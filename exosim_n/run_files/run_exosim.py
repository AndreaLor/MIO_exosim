from exosim_n.classes.options import Options
from exosim_n.classes.params import Params

from exosim_n.generate.gen_planet_xml_file import make_planet_xml_file
from exosim_n.run_files.recipe_1 import recipe_1
from exosim_n.run_files.recipe_2 import recipe_2
from exosim_n.run_files.recipe_3 import recipe_3
from exosim_n.run_files.recipe_4 import recipe_4
from exosim_n.run_files  import results

import numpy           as     np
import sys, time, os
import exosim_n 
from exosim_n.lib.exosim_n_lib import exosim_n_msg
import matplotlib.pyplot as plt
import exosim_n.pipeline as pipeline
import gc
import os
plt.style.use('classic')

def run(params_file):
  
  exosim_n_msg('exosim_n is running!\n', 1) 
  exosim_n_msg('User-defined input parameter file: %s\n  '%(params_file), 1)  
  exosim_n_path =  os.path.dirname((os.path.dirname(exosim_n.__file__)))
  exosim_n_msg('Reading instrument file ... \n', 1)
  
  opt = Options(filename='').opt

  exosim_n_msg('Apply user-defined adjustments to default... \n', 1)
  paths_file = '%s/exosim_n/input_files/exosim_n_paths.txt'%(exosim_n_path)
  params_file = '%s/exosim_n/input_files/%s'%(exosim_n_path, params_file)
  
  #first find the system from the input text file and then find the XML ICF file for system
  params_to_opt = Params(opt, params_file, 3)  
  system = (params_to_opt.params['obs_system'])
  ICF = '%s/exosim_n/systems/%s/%s.XML'%(exosim_n_path, system, system)

  # now use ICF
  opt = Options(ICF).opt    
  opt.exosim_n_path = exosim_n_path
  opt.__path__ = '%s/%s'%(exosim_n_path, 'exosim_n')
     
  #read in path information set by user
  params_to_opt = Params(opt, paths_file, 0)  
  paths = params_to_opt.params
  opt = params_to_opt.opt # adjust default values to user defined ones
  
  # read in input parameters for this simulation
  params_to_opt = Params(opt, params_file, 1)  
  input_params = params_to_opt.params
  opt = params_to_opt.opt # adjust default values to user defined ones   
  
  #read in path information set by user
  params_to_opt = Params(opt, paths_file, 0)  
  paths = params_to_opt.params
  opt = params_to_opt.opt # adjust default values to user defined ones
  if opt.common.output_directory.val == '__path__/output':
      opt.common.output_directory.val = '%s/output'%(exosim_n_path)
 
  opt.no_real = opt.simulation.sim_realisations.val
  opt.diagnostics = opt.simulation.sim_diagnostics.val
  opt.input_params = input_params  
  opt.system = opt.observation.obs_system.val
  opt.params_file_path = params_file 

  # select channel
  for ch in opt.channel:
      if ch.name == opt.observation.obs_channel.val:  
          opt.channel = ch

  if opt.channel.is_spec.val == 'True':
  
      opt.channel.is_spec.val = True
  elif opt.channel.is_spec.val == 'False':
      opt.channel.is_spec.val = False


  pl = opt.exosystem_params.planet_name.val
  make_planet_xml_file(opt, pl) 
 
  ex_file = '%s/exosim_n/exosystems/%s.xml'%(exosim_n_path, pl)
  opt3 = Options(ex_file).opt
  for key in opt3.__dict__:
     setattr(opt, str(key), opt3.__dict__[key])

    # if opt.exosystem_params.planet_use_database.val == 0:
    #     pl =  input_params['user_defined_planet_name']
  
  opt.sim_mode = opt.simulation.sim_mode.val
  opt.lab = '%s_%s'%(opt.observation.obs_channel.val,pl)
  
  
  if opt.sim_mode == 1:
          recipe  = recipe_1(opt)       
  if opt.sim_mode == 2:
          recipe  = recipe_2(opt)
  if opt.sim_mode == 3:
          recipe  = recipe_3(opt) 
  if opt.sim_mode == 4:
          recipe  = recipe_4(opt)      
          
  results_file = recipe.filename
  results.run(results_file)

if __name__ == "__main__":     
    
     run('exosim_n_input_params_ex1.txt')
 
    
