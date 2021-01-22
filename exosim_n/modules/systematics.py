"""
exosim_n

Systematics module

"""

from exosim_n.lib import exosim_n_lib, systematics_lib
from exosim_n.lib.exosim_n_lib import exosim_n_msg
import copy

def run(opt):
    
    opt = systematics_lib.gen_prnu_grid(opt)
    opt = systematics_lib.gen_systematic_grid(opt)
    
    opt.qe_original = copy.deepcopy(opt.qe)
    opt.qe_uncert_original = copy.deepcopy(opt.qe_uncert)
    opt.syst_grid_original = copy.deepcopy(opt.syst_grid)
 
    return opt

      
    
