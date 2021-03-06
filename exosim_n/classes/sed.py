import numpy as np
import sys 
from ..lib import exosim_n_lib


class Sed(object):
  """
    A class container for SED-like objects
    
    Attributes
    ----------
    sed : array
	  the SED array
    wl  : array
	  the wavelength array
  """
  
  def __init__(self, wl=None, sed=None):
    self.sed        = sed
    self.wl         = wl
  
  
  def rebin(self, wl):
    """
      Rebins the SED to the given wavelength
      
      Parameters
      __________
	wl: array
	    wavelength
    """
    self.wl, self.sed = exosim_n_lib.rebin(wl, self.wl, self.sed)
    
    