"""
exosim_n

Planet class

"""

import numpy as np
from exosim_n.classes import sed
from astropy import units as u
from astropy import constants as const
from exosim_n.lib.exosim_n_lib import exosim_n_msg, planck
import os, sys
from pytransit import QuadraticModel 


class Planet(object):
  """
    Instantiate a Planet class 
  """  
 
  def __init__(self, opt, exosystem):
    """
    Initialize a new Planet class 
    
    Parameters
    __________
    contrast_path		string
				path name of file containg the planet-star contrast spectrum
    """
    
    self.planet = exosystem
    self.opt =opt  
    
        
    if opt.exosystem_params.planet_spectrum_model.val == 'simple':
        self.simple_model() 
        exosim_n_msg ("simple planet spectrum chosen with fixed CR of %s"%(self.sed.sed[0]), 1)
    elif opt.exosystem_params.planet_spectrum_model.val =='file':
        self.file_model()  
        exosim_n_msg ("filed planet spectrum chosen", 1)
    else:
        exosim_n_msg('Error1 planet class: no compatible entry for planet_spectrum_model', 1)
        sys.exit()

    
            
  def simple_model(self):    #need to complete    
    if self.opt.observation.obs_type.val == 1:
        wl = np.arange(0.3,30,0.1)*u.um
        cr = np.array([( (self.planet.R).to(u.m)/(self.planet.star.R).to(u.m))**2]*len(wl))*u.dimensionless_unscaled
        self.cr = self.sed = sed.Sed(wl,cr)      
    elif self.opt.observation.obs_type.val == 2:
        wl = np.arange(0.3,30,0.1)*u.um
        star_flux =  np.pi*planck(wl, self.planet.star.T)*(  (self.planet.star.R).to(u.m) / (self.planet.star.d).to(u.m))**2
        planet_flux_em =  np.pi*planck(wl, self.planet.T)*(  (self.planet.R).to(u.m) / (self.planet.star.d).to(u.m))**2
        planet_flux_ref =  star_flux * (self.planet.albedo / 4.) *(  (self.planet.R).to(u.m) / (self.planet.a).to(u.m))**2
        planet_flux = planet_flux_em + planet_flux_ref 
        cr = planet_flux / star_flux
        self.cr = self.sed = sed.Sed(wl,cr)   
        # import matplotlib.pyplot as plt
        # plt.plot(wl, planet_flux)
        # plt.plot(wl, star_flux)
        # xxxx
    else:
        exosim_n_msg("Error2  planet class: no compatible entry for obs_type", 1)
        sys.exit() 


  def file_model(self):
      filename = self.opt.exosystem_params.planet_spectrum_file.val   
      if '/' in  filename:
          pass
      else:
          filename = '%s/exosim_n/data/planet_spectra/%s'%(self.opt.exosim_n_path, filename)
      
      exosim_n_msg('planet spectrum from file: %s'%(filename),1)
      try:
          aa = np.loadtxt(filename)
      except IOError:
          exosim_n_msg("Error6  planet class: No spectrum file found",  1) 
      if aa[:,0][0]< 1e-4: 
          wl=  aa[:,0]*1e6*u.um 
      else:
          wl=  aa[:,0]*u.um 
      cr = aa[:,1]*u.dimensionless_unscaled
      # import matplotlib.pyplot as plt
      # plt.figure(999999)
      # plt.plot(wl, cr)
      # xxxx
      self.cr = self.sed = sed.Sed(wl,cr)   
      

  def calc_T14(self, inc, a, per, Rp, Rs):
    b = np.cos(inc)*a/Rs
    rat = Rp/Rs
    self.t14 = per/np.pi* Rs/a * np.sqrt((1+rat)**2 - b**2)
    return self.t14

         
  def calc_gravity(self, Mp, Rp):
    Mp = 1*u.M_earth
    Rp = 1*u.R_earth
    
    Mp = Mp.to(u.kg)
    Rp = Rp.to(u.m)
    G = const.G
    g  = G*Mp/(Rp**2)
    return g 
    

      
  def get_t14(self, inc, a, period, planet_radius, star_radius):
    """ t14
    Calculates the transit time 
    
    Parameters
    __________
    inc:			scalar
				Planet oprbital inclination [rad]
    a:				scalar
				Semimajor axis [meters]
    period:			scalar
				Orbital period [seconds]
    planet_radius	:	scalar
				Planet radius [meters]
    star_radius	:		scalar
				Star radius [meters]
    
    Returns
    __________
    transit duration : float
	Returns the transit duration [seconds]
    
    Notes
    _____
    Seager, S., & Mallen-Ornelas, G. 2003, ApJ, 585, 1038
      
    
    """
    impact_parameter = np.cos(inc)*a/star_radius
    dtmp = 1+planet_radius/star_radius
    if impact_parameter < dtmp:
      self.t14 = period/np.pi * \
	    star_radius/a * \
	    np.sqrt(dtmp**2 - impact_parameter**2)
    else:
      #print "WARNING: planet not transiting"
      self.t14 = np.nan
    return self.t14


  def get_orbital_phase(self, t14, period, N=1000):
    f = t14/period
    self.phi = np.linspace(-f.item(), f.item(), N)*f.units
    return self.phi
    
  def eccentric(self, phi, inc, ecc, omega, a, period, star_radius):
    """ eccentric
    Implements an eccentric orbit and calculates the projectedseparation (z)
    
    Parameters
    __________
    phi : 	array
      orbital phase 
    inc :	float
      orbital inclination [radiants]
    ecc : 	float
      orbital eccentricity [radiants]
    omega :	float
      argument of periastron [radiants]
    a	:	float
      semimajor axis [meters]
    star_radius : 	float
      star radius [meters]
    
    Returns
    _______
    z	: 	array
      orbital separations
    
    Notes
    _____
    This implementation written by Ingo Waldman 2012
    """
    theta = 2.0 * np.pi * phi 
    aR    = a / star_radius
    
    if ecc < 1e-5:
      # calculating z for a circular orbit
      self.z = aR * np.sqrt(1-((np.cos(theta))**2.*(np.sin(inc))**2.))
      return self.z
    
    # calculating z for eccentric orbits
    n = len(theta)
    E = np.zeros(n)
    ecc2 = np.sqrt((1.+ecc)/(1.-ecc))
    fref = (np.pi / 2.) - omega #setting reference point for the true anomaly
    
    Eref = 2. * np.arctan(1./ecc2 * np.tan(fref/2.)) #reference point for eccentric anomaly
    if Eref < (-np.pi/2.):
      Eref = Eref + 2. * np.pi
      
    Mref = Eref - (ecc * np.sin(Eref)) #reference point for mean anomaly
   
    for i in range(n): 
      # calculating eccentric anomaly using Newton-Rhapson method
      Mtmp = theta[i] + Mref #calculating mean anomaly for phase point
      Etmp = Mtmp
      for j in range(10):
        Etmp = Etmp + ((Mtmp + ecc*np.sin(Etmp) - Etmp) / (1.-ecc*np.cos(Etmp)))
      E[i] = Etmp
    
      # calculating true anomaly
      f = 2.*np.arctan(ecc2*np.tan(E/2.))
      # calculating distance from true anomaly as fraction
      r_frac = (1.-ecc**2.)/(1. + ecc*np.cos(f))
      # computing z
      self.z = 1. - (np.sin(inc)**2.*np.sin(f+omega)**2.)
      self.z = aR*r_frac*np.sqrt(self.z)

    return self.z
    


