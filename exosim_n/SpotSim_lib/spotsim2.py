# -*- coding: utf-8 -*-
"""
Created on Sat Dec 24 12:50:40 2016

@author: c1341133
"""
import numpy as np
import matplotlib.pyplot as plt
from photutils import CircularAperture as CircAp
from photutils import aperture_photometry as ApPhot
#import scipy.constants as cnst
from scipy.misc import imrotate
from scipy.ndimage.interpolation import shift
#from scipy import interpolate
#from scipy import signal
#import numpy.random as rnd
import time
#import multiprocessing as mp
import pytransit
from matplotlib.patches import Ellipse
from scipy import ndimage
import os


#==============================================================================
# Primary Transit Class:
#==============================================================================



class SpotSim():
  def __init__(self, Tstar, Tspot, Tfac, Q, ImHalfSizePx, RstarPx, SUR, opt):  
            
   
    print "Using set files?", opt.SpotSim_params['useFileSpots']
    self.filling = opt.SpotSim_params['filling']    
        
    self.Tstar = Tstar
    self.RstarPx = RstarPx
#    self.ImHalfSizePx = pix*3/2
    self.ImHalfSizePx = ImHalfSizePx
    self.Tspot = Tspot
    self.Tfac = Tfac
    self.b = opt.b
    self.useSpotfile = opt.SpotSim_params['useFileSpots']  
    self.usePhoenix = opt.SpotSim_params['usePhoenix']
    self.realization = opt.realization
#    self.spotfile = opt.spotfile    
#    self.facfile = opt.facfile  
    self.Q =  Q
    self.path = opt.__path__ 
    self.addFacs = opt.SpotSim_params['addFacs']    
    self.SpatDist = opt.SpotSim_params['spatDist']    
    self.SUR = SUR
    
    
  def planck(self, w1, T):
    a = 1.191042768e8    #*pq.um**5 *pq.W/ pq.m**2 /pq.sr/pq.um
    b = 14387.7516       #*1*pq.um * 1*pq.K
    x = b/(w1*T)
    bb = a/w1**5 / (np.exp(x) - 1.0)
    return bb
  

  def create_spotted_disc(self): 
    
    Xspot =[]; Yspot = [] ; Dspot =[] ; Dfac =[]
    for i in range (len(self.spot_positions)):
        Xspot.append(self.spot_positions[i]['x']*self.RstarPx)
        Yspot.append(self.spot_positions[i]['y']*self.RstarPx)
        Dspot.append(self.spot_positions[i]['r']*self.RstarPx*2)
        Dfac.append(self.spot_positions[i]['r']*self.RstarPx*2*np.sqrt(self.Q + 1))
    Xspot = np.array(Xspot); Yspot = np.array(Yspot)
    Dspot = np.array(Dspot); Dfac = np.array(Dfac)       
    r = np.sqrt(Xspot**2 + Yspot**2) # radial distance of spot from centre in pixels
    a = 1.0 * Dspot ; b = a*np.sin(np.arccos(r/self.RstarPx))  # change in width of ellipse due to projection    
    a_fac = 1.0 * Dfac ; b_fac = a_fac*np.sin(np.arccos(r/self.RstarPx))  # change in width of ellipse due to projection
    
    angles=[]
    for i in range(len(self.spot_positions)):  
        Xspot0 = Xspot[i]  # find rotation angle theta and prevent bugs due to 0 vals
        Yspot0 = Yspot[i]
        if Xspot[i] == 0:
            Xspot0 = 1e-3
        if Yspot[i] ==0:
            Yspot0 = 1e-3
        if Xspot[i] == 0 and Yspot[i] == 0:
            theta = 0
        else:
            theta = -np.arctan(Xspot0*1.0/Yspot0)*360./(2*np.pi)
        angles.append(theta)        
#    self.dpi_set = 118
    self.dpi_set=62.3
    print "DPI SET", self.dpi_set
    fig  = plt.figure('starXX', figsize=(10,10), dpi = self.dpi_set)
    ax = plt.subplot(111, aspect='equal')
    e = Ellipse((0,0), self.RstarPx*2-1.5, self.RstarPx*2-1.5,0) # -1.5 needed as a correction factor
    e.set_clip_box(ax.bbox); e.set_alpha(0.1); ax.add_artist(e)        
    plt.xticks([]); plt.yticks([])
    plt.xlim(-self.ImHalfSizePx, self.ImHalfSizePx); plt.ylim(-self.ImHalfSizePx, self.ImHalfSizePx)
    plt.tight_layout(); plt.axis('off')
    fig.canvas.draw()
    data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    data = np.array(data)[...,1]
    data = np.where(data==data.max(), 0.0, 1.0)[11:-11, 11:-11]
    plt.close()    
    self.star_mask = data
#    plt.figure(888)
#    plt.imshow(data)       
    print data[:,data.shape[1]/2].sum(), data.shape
    if data.shape[0] != self.ImHalfSizePx*2 or data[:,data.shape[1]/2].sum() != self.RstarPx*2:
        print "ERROR IN SIZING - look at code!"
    
    if self.useSpotfile == 0:

        fig  = plt.figure('spotsXX', figsize=(10,10), dpi = self.dpi_set) 
        ax = plt.subplot(111, aspect='equal')
        for i in range(len(self.spot_positions)):
            e = Ellipse((Xspot[i], Yspot[i]), a[i]-1.5, b[i]-1.5, angles[i])     
            e.set_clip_box(ax.bbox); e.set_alpha(0.1); ax.add_artist(e)        
        plt.xticks([]); plt.yticks([])
        plt.xlim(-self.ImHalfSizePx, self.ImHalfSizePx); plt.ylim(-self.ImHalfSizePx, self.ImHalfSizePx)
        plt.tight_layout(); plt.axis('off')
        fig.canvas.draw()
        data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
        data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
        data = np.array(data)[...,1]
        data = np.where(data==data.max(), 0, 1.0)[11:-11, 11:-11]
        plt.close()    
        self.spots_mask = data
#        plt.figure(8889)
#        plt.imshow(data)    
            
        fig  = plt.figure('facsXX', figsize=(10,10), dpi = self.dpi_set) 
        ax = plt.subplot(111, aspect='equal')
        for i in range(len(self.spot_positions)):
            e = Ellipse((Xspot[i], Yspot[i]), a_fac[i]-1.5, b_fac[i]-1.5, angles[i])     
            e.set_clip_box(ax.bbox); e.set_alpha(0.1); ax.add_artist(e)        
        plt.xticks([]); plt.yticks([])
        plt.xlim(-self.ImHalfSizePx, self.ImHalfSizePx); plt.ylim(-self.ImHalfSizePx, self.ImHalfSizePx)
        plt.tight_layout(); plt.axis('off')
        fig.canvas.draw()
        data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
        data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
        data = np.array(data)[...,1]
        data = np.where(data==data.max(), 0.0, 1.0)[11:-11, 11:-11]
        plt.close()
        self.facs_mask = np.where(self.spots_mask+data==1.0,1.0,0.0)
#        plt.figure(8890)
#        plt.imshow(self.facs_mask)
        
    elif self.useSpotfile == 1:      
        self.spots_mask =  self.spotfile[...,self.realization]  
        self.facs_mask =self.facfile[...,self.realization]

    return np.sum(self.spots_mask)/np.sum(self.star_mask)


  def limb_darkening(self, I0, c1, c2, Rimg, pix):
      x = np.linspace(-1, 1, 2*Rimg)*Rimg/pix
      y = np.linspace(-1, 1, 2*Rimg)*Rimg/pix
      xx, yy = np.meshgrid(x, y)
      arg = np.sqrt(xx**2 + yy**2)
      u = np.sqrt(1 - arg**2)
      arg_l = I0 * ( 1 - c1*(1-u) - c2*(1-u)**2)
      arg_l = np.where(np.isnan(arg_l), 0, arg_l)
      return arg_l
      
      
  def extra_spots_flash(self, facs_only):
   
    if facs_only == 0:   
        FF_in = 100* np.sum(self.spots_mask)/ np.sum(self.star_mask)
        n = int(0.6*(len(self.Area_spots) * (self.filling - FF_in)/FF_in))
        print "Flash addition of spots/facs"
        spotX = self.spots_mask*0; facX = self.facs_mask*0
    if facs_only == 1:
        Q_in = np.sum(self.facs_mask)/ np.sum(self.spots_mask)
        n = int(0.6*(len(self.Area_spots) * (self.Q - Q_in)/Q_in))
        
        
        print "Flash addition of facs"
        facX = self.facs_mask*0   
    if self.SpatDist == 'Equatorial':
        if self.filling <= 70.0 and self.filling>=35 :
            n  = int(n*4.0)
            if facs_only == 1:
                n = int(n*2)
        if self.filling >= 15.0 and self.filling<=30 :
            n  = int(n*1.6)
            if facs_only == 1:
                n = int(n*2) 


    if self.SpatDist == 'Longitudinal':
        if self.filling <= 17.0 and self.filling>=7.5 :
            if facs_only == 1:
                n = int(n*4)
        if self.filling >= 15.0 and self.filling<=29 :
            n  = int(n*4.0)
#            if facs_only == 1:
#                n = int(n*3.5)    
        if self.filling <= 70.0 and self.filling>=30 :
            n  = (n*4.0)
            if facs_only == 1:
                n = (n*1.2)       
        if self.filling == 5.0:
            if facs_only == 1:
                n = (n*3.5)       

        if self.filling == 10.0:
            if facs_only == 1:
                n = (n*1.5)   
    
        n = int(n)                
                
            
    for i in range(50): 
        print i
        if n == 9:
            break
            print "error in rescaling extra spots or facs" 
            xxxxxxxx
            
        elif n<0:
            n ==1       

                        
        Spot_areas = np.random.choice(self.A, n, p=self.A_weights)*self.SUR*1e-6
        Area_hemi = 2*np.pi*self.RstarPx**2 #Area of hemisphere in pixels
        Spot_areas *=Area_hemi # Area of whole spots in pixels
  
        if self.SpatDist == 'Uniform':
            angles1 = np.random.uniform(0,np.pi, n)  # random angular position in longintude
            angles2 = np.random.uniform(0,np.pi, n)  
        elif self.SpatDist == 'Equatorial':
            angles1 = np.random.uniform(0,np.pi, n)  # random angular position in longintude
            angles2 = np.random.normal(np.pi/2.,1*self.band_factor, n)  # equatorial qith 10 degree sd                
        elif self.SpatDist == 'Longitudinal':
            
            angles1_0=[]
            angles2_0=[]
            x_Spot_areas =[]
            for q in range(self.groups):
                cen = self.cen1 + q*(2*np.pi)/self.groups
                print "CEN", cen, q       
                angles1 = np.random.normal(cen, 1*self.band_factor, n)
                angles1 = np.where(angles1 <= 2*np.pi, angles1, angles1-2*np.pi)
                angles1 = np.where(angles1 >=0, angles1, 2*np.pi-angles1)
                angles2 = np.random.normal(np.pi/2.,1*self.band_factor, n) # equatorial qith 10 degree sd        
                xx_Spot_areas =  Area_hemi* (np.random.choice(self.A, n , p=self.A_weights)*self.SUR*1e-6)
        
                idx = np.argwhere(angles1<=np.pi)
                angles1_0 = angles1_0 + angles1[idx].tolist()
                angles2_0 = angles2_0 + angles2[idx].tolist()
                x_Spot_areas = x_Spot_areas  + xx_Spot_areas[idx].tolist()
            angles1 = np.array(angles1_0)
            angles2 = np.array(angles2_0)
            Spot_areas = np.array(x_Spot_areas)
 
        R_spot = ((Spot_areas/np.pi)**0.5) # array of star spot radii relative to star        
        R_fac = R_spot*np.sqrt(self.Q + 1)
        if facs_only == 1:
            R_fac = R_spot*np.sqrt(self.Q)
        Dfac = R_fac*2 ;  Dspot = R_spot*2
         
        Xspot = np.cos(angles1)*self.RstarPx; Yspot =  np.cos(angles2)*self.RstarPx
        r = np.sqrt(Xspot**2 + Yspot**2) # radial distance of spot from centre in pixels
        a = 1.0 * Dspot ; b = a*np.sin(np.arccos(r/self.RstarPx))  # change in width of ellipse due to projection   
        a_fac = 1.0 * Dfac ; b_fac = a_fac*np.sin(np.arccos(r/self.RstarPx))  # change in width of ellipse due to projection
        if facs_only == 1:
            a = a_fac
            b = b_fac
           
        angles=[]
        for i in range(len(angles1)):  
            
            Xspot0 = Xspot[i]  # find rotation angle theta and prevent bugs due to 0 vals
            Yspot0 = Yspot[i]
            if Xspot[i] == 0.0:
                Xspot0 = 1.0e-1
            if Yspot[i] ==0.0:
                Yspot0 = 1.0e-1
            if Xspot[i] == 0.0 and Yspot[i] == 0.0:
                theta = 0.0
            else:
                theta = -np.arctan(Xspot0*1.0/Yspot0)*360./(2*np.pi)
            angles.append(theta)    

        fig  = plt.figure('spotsXX', figsize=(10,10), dpi = self.dpi_set) 
        ax = plt.subplot(111, aspect='equal')
        for i in range(len(angles1)):
            e = Ellipse((Xspot[i], Yspot[i]), a[i]-1.5, b[i]-1.5, angles[i])     
            e.set_clip_box(ax.bbox); e.set_alpha(0.1); ax.add_artist(e)        
        plt.xticks([]); plt.yticks([])
        plt.xlim(-self.ImHalfSizePx, self.ImHalfSizePx); plt.ylim(-self.ImHalfSizePx, self.ImHalfSizePx)
        plt.tight_layout(); plt.axis('off')
        fig.canvas.draw()
        data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
        data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
        data = np.array(data)[...,1]
        data = np.where(data==data.max(), 0.0, 1.0)[11:-11, 11:-11]
        plt.close()
        if facs_only == 0:
            spotX = data + self.spots_mask
            spotX = np.where(spotX==0.0,0.0,1.0)
            FF_out = np.sum(spotX)/ np.sum(self.star_mask)
        
            if FF_out <= self.filling/100.:
                
                fig  = plt.figure('facsXX', figsize=(10,10), dpi = self.dpi_set) 
                ax = plt.subplot(111, aspect='equal')
                for i in range(len(angles1)):
                    e = Ellipse((Xspot[i], Yspot[i]), a_fac[i]-1.5, b_fac[i]-1.5, angles[i])     
                    e.set_clip_box(ax.bbox); e.set_alpha(0.1); ax.add_artist(e)        
                plt.xticks([]); plt.yticks([])
                plt.xlim(-self.ImHalfSizePx, self.ImHalfSizePx); plt.ylim(-self.ImHalfSizePx, self.ImHalfSizePx)
                plt.tight_layout(); plt.axis('off')
                fig.canvas.draw()
                data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
                data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
                data = np.array(data)[...,1]
                data = np.where(data==data.max(), 0.0, 1.0)[11:-11, 11:-11]
                plt.close()
                facX = data + self.facs_mask + self.spots_mask
                facX = np.where(facX==0.0,0.0,1.0)
                spot_facX = facX + spotX
                self.facs_mask = np.where(spot_facX==1.0,1.0,0)  
                self.spots_mask = spotX
                
                break
            else:
                print "OVERSHOT FF ->>> repeating"
                n = int(n-n*0.1)
            
        elif facs_only == 1:

            facX = data + self.facs_mask + self.spots_mask
            facX = np.where(facX==0.0,0.0,1.0)            
            spot_facX = facX + self.spots_mask
            facX = np.where(spot_facX==1.0,1.0,0.0)
            Q_out = np.sum(facX)/ np.sum(self.spots_mask)
            if self.SpatDist == 'Equatorial':
                if Q_out <= self.Q and Q_out > 0.7*self.Q:
                    self.facs_mask =  facX
                    break
                else:
                    print "OVERSHOT Q ->>> repeating"
                    n = int(n-n*0.1)
            else:
                if Q_out <= self.Q:
                    self.facs_mask =  facX
                    break
                else:
                    print "OVERSHOT Q ->>> repeating"
                    n = int(n-n*0.1)
            
    return self.spots_mask     
      

 
  def extra_spots_top_up(self, facs_only):
      
    if facs_only == 0:  
        all_spots = self.spots_mask ; all_spots0 = all_spots*0
        all_facs = self.facs_mask;  all_facs0 = all_facs*0
        print "Adding extra spots to exact filling factor..."
    elif facs_only ==1:
        print "Adding extra facs to exact Q factor..."
        all_pix = self.spots_mask + self.facs_mask + self.star_mask
        all_pix0 = np.where(all_pix==1.0,1.0,0.0)
        idx = np.argwhere(all_pix0 ==1)
        x=[];y=[]
        for i in range (len(idx)):
            x.append(idx[i][1])
            y.append(idx[i][0])
        x=np.array(x); y= np.array(y)
        x=np.array(x)-300; y= np.array(y)-300
        x_y_list =idx        
        r_list= np.sqrt(x**2.0 + y**2.0)/self.RstarPx    
        r0 = np.linspace(0,.99,600)
        dr = r0[1]-r0[0]
        w0 =( dr * 1/(np.sqrt(1-r0**2)) )  /  np.sum( dr* 1/(np.sqrt(1-r0**2)) )          
        weights=[]
        for i in range (len(r_list)):
            idx = np.argwhere(r0<=r_list[i])[-1]
            weights.append(w0[idx])  
        weights = np.array(weights) ; weights = weights/np.sum(weights)
        r_list = np.reshape(r_list, len(r_list)) ; weights = np.reshape(weights, len(weights))  
        
        if self.SpatDist == 'Equatorial':
#            gauss_weight = np.exp(-(y**2/(2*(self.RstarPx*np.tan(self.band_factor*1))**2)))           
            gauss_weight = np.exp(-( (np.arcsin(y/self.RstarPx))**2 / (1*(self.band_factor*1)**2)   ) )
#            plt.figure('Gaussian weighting')
#            plt.plot(y, gauss_weight)
            weights *= gauss_weight ; weights = weights/np.sum(weights)  


        if self.SpatDist == 'Longitudinal':
        
            alpha = (2*np.pi/ self.groups)            
            cen = self.cen1 + np.arange(0,np.pi*2, alpha)            
            cen = np.where(cen > 2*np.pi, cen-2*np.pi, cen)
            cen = np.where(cen < 0, 2*np.pi-cen, cen)
            cen = cen[np.argwhere(cen<np.pi)]
            idx = np.argwhere(cen >= np.pi/2.)
            cen = -(np.pi/2. - cen)
            x_cen = - np.sin(cen)*self.RstarPx
            
#            fact = 1
#            if self.filling >=25.0 and facs_only == 1:
#                fact = 1
#            if self.filling >=50.0:
#                fact = 4
                       
            print x_cen
            long_weight = np.zeros(len(x))
                    
            for i in range(len(cen)):

                cen0 = abs(cen[i]).item()                      
                r = np.arange(-1000,1000,1)
                weight  = np.exp(- (np.arcsin ( (x - x_cen[i])  / (self.RstarPx*np.cos(cen0))  )  )**2 / ((self.band_factor*1)**2)    )
                weight[np.isnan(weight)] = 0
                long_weight += weight
            gauss_weight = np.exp(-( (np.arcsin(y/self.RstarPx))**2 / ((self.band_factor*1)**2)   ) )

            weights *= long_weight*gauss_weight ; weights = weights/np.sum(weights) 
         

        all_spots = self.spots_mask
        all_spots0 = all_spots*0
        all_facs = self.facs_mask
        all_facs0 = all_facs*0    


    ct = 0
    for i in range(1000):
        
        Spot_area = np.random.choice(self.A, 1, p=self.A_weights)[0].item()*self.SUR*1e-6
        Area_hemi = 2*np.pi*self.RstarPx**2 #Area of hemisphere in pixels
        Spot_area *=Area_hemi # Area of whole spots in pixels
        R_spot = ((Spot_area/np.pi)**0.5)    
               
        if facs_only == 1:
            R_fac = R_spot*np.sqrt(self.Q)
#            if self.filling <=2.0:
#               R_fac = R_spot            
            idx0  = np.arange(0, len(x_y_list), 1)     
            idx  = np.random.choice(idx0, 1, p=weights)
            idx =  int(idx.item())
            Xspot = x_y_list[idx][1]-300.0; Yspot = x_y_list[idx][0]-300.0

        elif facs_only == 0:
            R_fac = R_spot*np.sqrt(self.Q + 1) 
            if self.SpatDist == 'Uniform':
                angle1 = np.random.uniform(0,np.pi)  # random angular position in longintude
                angle2 = np.random.uniform(0,np.pi)
            elif self.SpatDist == 'Equatorial':
                angle1 = np.random.uniform(0,np.pi)  # random angular position in longintude
                angle2 = np.random.normal(np.pi/2.,1*self.band_factor)  # equatorial qith 10 degree sd                     
            elif self.SpatDist == 'Longitudinal':
                    if ct >= self.groups: ct = 0
                    cen = self.cen1 + ct*(2*np.pi)/self.groups
                    ct +=1
                    angle1 = np.random.normal(cen, 1*self.band_factor)
                    if angle1 > 2*np.pi: angle1 = angle1-2*np.pi
                    if angle1 < 0: angle1 = 2*np.pi-angle1
                    angle2 = np.random.normal(np.pi/2.,1*self.band_factor) # equatorial qith 10 degree sd                    
                    if angle1 > np.pi:                    
                        angle1 = angle1  - np.round((np.pi)/(2*np.pi/ self.groups)) * (2*np.pi/ self.groups)
                    if angle1 < 0:
                        angle1 += (2*np.pi/ self.groups)
                    if angle1 <0 or angle1 >np.pi:
                        print "ERROR IN ANGLE!!!!!!!!!!!", angle1
            Xspot = np.round(np.cos(angle1)*self.RstarPx); Yspot = np.round(np.cos(angle2)*self.RstarPx)

     
        r = np.sqrt(Xspot**2.0 + Yspot**2.0) # radial distance of spot from centre in pixels
        a = 1.0   
        b = a*np.sin(np.arccos(r/self.RstarPx)) # change in width of ellipse due to projection
          
        Xspot0 = Xspot  # find rotation angle theta and prevent bugs due to 0 vals
        Yspot0 = Yspot
        if Xspot == 0:
            Xspot0 = 1e-3
        if Yspot ==0:
            Yspot0 = 1e-3
        if Xspot == 0 and Yspot == 0:
            theta = 0
        else:
            theta = np.arctan(Xspot0*1.0/Yspot0)*360./(2*np.pi)
           
        y_s = np.linspace(-self.ImHalfSizePx+0.5, self.ImHalfSizePx-0.5, 2*self.ImHalfSizePx) #set up ellipse
        x_s = np.linspace(-self.ImHalfSizePx+0.5, self.ImHalfSizePx-0.5, 2*self.ImHalfSizePx)
        xx_s, yy_s = np.meshgrid(x_s, y_s)
        arg_s = np.sqrt((xx_s/a)**2 + (yy_s/b)**2) 
        
        if facs_only == 0:
            spot = np.where(arg_s> R_spot, 0.0, 1.0)
            fac = np.where(arg_s> R_fac, 0.0, 1.0)
    
            spot0 = imrotate(spot, theta,interp='bilinear') # rotate spot and facula 
            fac0 = imrotate(fac, theta,interp='bilinear') # rotate spot and facula      
            
            spotMask = np.where(spot0==0, spot0, spot.max())# account for artifact from rotating
            facMask = np.where(fac0==0, fac0, fac.max())# account for artifact from rotating  
            spot_mask = np.ones( (2*int(self.ImHalfSizePx)+1, 2*int(self.ImHalfSizePx)+1))
            fac_mask = np.ones( (2*int(self.ImHalfSizePx)+1, 2*int(self.ImHalfSizePx)+1))              
            spot_mask *= shift(np.where(spotMask>0, -1, spotMask)+1, 
                                      (Yspot, Xspot), cval = 1)
            fac_mask *= shift(np.where(facMask>0, -1, facMask)+1, 
                                      (Yspot, Xspot), cval = 1)                                                                  
            spot_mask = np.round(spot_mask)
            fac_mask = np.round(fac_mask)        
            spot_mask = np.where(spot_mask==0.0,1.0,0.0)
            fac_mask = np.where(fac_mask==0.0,1.0,0.0)            
            all_spots0 +=  spot_mask
            all_facs0 +=  fac_mask
            all_spots0 = np.where(all_spots0==0,0.0,1.0)  
            all_facs0 = np.where(all_facs0==0,0.0,1.0)                   
            all_spots += spot_mask
            all_spots = np.where(all_spots==0,0.0,1.0)
            
            FF = np.sum(all_spots)/ np.sum(self.star_mask)
#            print FF, np.sum(all_spots), np.sum(self.star_mask)
            
            if FF >= self.filling/100:
                facs = self.spots_mask+self.facs_mask +  all_facs0
                facs = np.where(facs==0.0,0.0,1.0)   
                spots =   self.spots_mask + all_spots0
                spots = np.where(spots==0.0,0.0,1.0)  
                spot_fac = spots + facs    
                self.facs_mask = np.where(spot_fac==1.0,1.0,0.0)
                self.spots_mask = spots  
                
                plt.figure('facs addedl')
                plt.imshow(all_facs0)
            
                plt.figure('spots added')
                plt.imshow(all_spots0)  
#                
                break
        if facs_only == 1:
            fac = np.where(arg_s> R_fac, 0.0, 1.0)
            fac0 = imrotate(fac, theta,interp='bilinear') # rotate spot and facula      
            facMask = np.where(fac0==0, fac0, fac.max())# account for artifact from rotating  
            fac_mask = np.ones( (2*int(self.ImHalfSizePx)+1, 2*int(self.ImHalfSizePx)+1))              
            fac_mask *= shift(np.where(facMask>0, -1, facMask)+1, (Yspot, Xspot), cval = 1)                                                                  
            fac_mask = np.round(fac_mask)        
            fac_mask = np.where(fac_mask==0,1.0,0.0)            
            all_facs0 += fac_mask
  
            facX = all_facs0 + self.facs_mask + self.spots_mask
            facX = np.where(facX==0.0,0.0,1.0)            
            spot_facX = facX + self.spots_mask
            facX = np.where(spot_facX==1.0,1.0,0.0)
            QX = np.sum(facX)/np.sum(self.spots_mask) 
#            print QX
            if QX>= self.Q:
                self.facs_mask = facX
#                plt.figure('facs added2')
#                plt.imshow(all_facs0)
                break  
    
    return
    
  def size_dist_new(self, c_s, SUR,size_type):
      
      if self.useSpotfile == 0: 

            if size_type == 'LogN':  
              print "LOG-NORMAL SIZE DISTRIBUTION"  
                   
              Std = 4.0
              A_mean = 29.979*(c_s*0.01/SUR)+0.4769
              
        
              dA = 5.0 
              A_init = np.arange(3.0, 5000, dA)
              Grad_max = 1.0
              
              Grad = Grad_max * np.exp(-(np.log(A_init)-np.log(A_mean))**2/(2*np.log(Std))) 
              
              idx = np.argwhere(Grad < 0.001*Grad.max())[0].item()
              A_max =  A_init[idx] # A_max is area when the grad function falls to < 0.1% of its maximum
              
              self.A = np.arange(3.0, A_max, dA)
              Grad = Grad_max * np.exp(-(np.log(self.A)-np.log(A_mean))**2/(2*np.log(Std))) 
              self.A_weights = Grad/np.sum(Grad)
                             
              Total_area = 0
              Area_list = []
              for i in range(10000000):
                Spot_area = np.random.choice(self.A, 1, p=self.A_weights)[0].item()*SUR*1e-6
                Total_area += Spot_area
                if Total_area/(c_s/100.) >= 1.0:
                    break
                else:
                    Area_list.append(Spot_area)
                    
                    
            if size_type == 'Uniform':  
                print "UNIFORM SIZE DISTRIBUTION"  
                                      
             
                theta_deg = 2.0 # radius of spot in degrees
                theta_rad = theta_deg*2*np.pi/360.
                r =  np.sin(theta_rad)
                Area_spot = r**2 / 2                 # check this...
                N_spots = np.round((c_s/100.)/Area_spot)
                Area_list = np.array( int(N_spots)*[Area_spot] )
                self.A = Area_list / (self.SUR*1e-6)
                self.A_weights = np.ones(len(Area_list))/ np.sum(np.ones(len(Area_list)))                
   
            self.Area_spots = Area_list                    
                                       
    
      Area_list = np.array(Area_list) 


      
      self.Area_spots = Area_list    
        
                
  def size_dist(self,Filling0, SUR, size_type):
      
      if self.useSpotfile == 0: 

            if self.SpatDist == 'Uniform':
                corr = 1.0
                if Filling0 >= 0.1: corr = 0.65
                if Filling0 >= 1.0: corr = 0.85
                if Filling0 >= 5.0: corr = 0.94
                if Filling0 >= 10.0: corr = 0.98              
                if Filling0 >= 25.0: corr = 1.2
                if Filling0 >= 50.0: corr = 1.5 
            elif self.SpatDist == 'Equatorial' :
                corr = 1.0
                if Filling0 >= 0.1: corr = 0.4
                if Filling0 >= 1.0: corr = 0.4
                if Filling0 >= 5.0: corr = 0.4
                if Filling0 >= 10.0: corr = 0.52              
                if Filling0 >= 25.0: corr = 0.70
                if Filling0 >= 50.0: corr = 1.1
            elif self.SpatDist == 'Longitudinal':
                corr = 1.0
                if Filling0 >= 0.1: corr = 0.4
                if Filling0 >= 1.0: corr = 0.4
                if Filling0 >= 5.0: corr = 0.6
                if Filling0 >= 10.0: corr = 0.7           
                if Filling0 >= 25.0: corr = 0.8
                if Filling0 >= 50.0: corr = 3.5   
                          
            
            coverage = np.sqrt(2)*self.filling
            
            if self.SpatDist == 'Longitudinal':
                corr =corr * 0.5   
          
            if size_type == 'LogN':  
                print "LOG-NORMAL SIZE DISTRIBUTION" 
     
                A_mean = 0.57
                coverage = np.sqrt(2.0)*self.filling
                                
                print "Input filling factor %", self.filling
                print "Calculated goal coverage %", coverage
                
                # these are from table 6 in Solanki and Unruh
                hem_cov = np.array([0.0003,0.003,0.014,0.07,0.28,0.62])
                sigma = np.array([3.75,4.0,5.04,10.4,30.3,62.0])
                #        gradmax = np.array([4.5,38.5,106,106,41.2,20.8])
                
                r = 1 #using linear fit but r=3 gives closer values to the above - but linear adequate re accuracy
                z = np.polyfit(hem_cov, sigma, r)
                
#                x = np.linspace(0,100,100)
#                y = 0
#                for i in range (0,r+1):
#                    y = y + z[i]*(x/100.)**(r-i)
#                Std = y
                y = 0
                for i in range (0,r+1):
                    y = y + z[i]*(coverage/100.)**(r-i)
                Std = y
                
                print ""
                print "For log-normal distribution using relationship from Solanki and Unruh:"
                print "for coverage of ", coverage, "%, calculated Std is ", Std
                print "keeping A_mean fixed to:", A_mean
                print ""                
                
                dA = 5.0 
                A = np.arange(3.0, 5000, dA)
                Grad_max = 1.0
                Grad = Grad_max * np.exp(-(np.log(A)-np.log(A_mean))**2/(2*np.log(Std))) 
                idx = np.argwhere(Grad < 0.001*Grad.max())[0].item()
                A_max =  A[idx] # A_max is area when the grad function falls to < 0.1% of its maximum
                self.A = np.arange(3.0, A_max, dA)
                Grad = Grad_max * np.exp(-(np.log(self.A)-np.log(A_mean))**2/(2*np.log(Std))) 
                self.A_weights = Grad/np.sum(Grad)
                     
                Total_area = 0
                Area_list = []
                for i in range(10000000):
                    Spot_area = np.random.choice(self.A, 1, p=self.A_weights)[0].item()*SUR*1e-6
                    Total_area += Spot_area
                    if Total_area/(coverage*corr/100.) >= 1.0:
                        break
                    else:
                        Area_list.append(Spot_area)
                
                Area_list = np.array(Area_list);  

 
            if size_type == 'Uniform':  
                   
                print "Input filling factor %", self.filling
                print "Calculated goal coverage %", coverage
             
                theta_deg = 2.0 # radius of spot in degrees
                theta_rad = theta_deg*2*np.pi/360.
                r =  np.sin(theta_rad)
                Area_spot = r**2 / 2                
                N_spots = np.round((coverage*corr/100.)/Area_spot)
                Area_list = np.array( int(N_spots)*[Area_spot] )
                self.A = Area_list / (self.SUR*1e-6)
                self.A_weights = np.ones(len(Area_list))/ np.sum(np.ones(len(Area_list)))                
   
            self.Area_spots = Area_list

#            print "XXXXXXXXXXXXXXXXXXXXXXXX"
#            print Area_list
#            print Area_list.shape
#            
#
#            plt.figure('spot hisotgram')
#            plt.hist(Area_list/(SUR*1e-6), bins = 10)

      

      else:
         self.Area_spots = []
                        
      return self.Area_spots      
        

  def spatial_dist(self, RstarPx, SpatDist): 
      
      if self.useSpotfile == 0:
        # Area of whole spots (umbra+penumbra) in 10^-6 of the hemisphere
        Area_hemi = 2*np.pi*RstarPx**2 #Area of hemisphere in pixels
        Area_spots = self.Area_spots
        Area_spots *=Area_hemi # Area of whole spots in pixels
        r_spots = ((Area_spots/np.pi)**0.5)/RstarPx # array of star spot radii relative to star
        
     
        if SpatDist == 'Uniform':
            print "spatial distribution is UNIFORM"
            angles1 = np.random.uniform(0,np.pi,len(r_spots))  # random angular position in longintude
            angles2 = np.random.uniform(0,np.pi,len(r_spots)) # random angular position in latitutude
        elif SpatDist == 'Equatorial':
            print "spatial distribution is EQUATORIAL"            
#            self.band_factor = -0.0003*self.filling**2 + 0.0241*self.filling + 0.7876
#            self.band_factor  = 0.02*self.filling + 0.8
            self.band_factor  = 1.2*(0.02*self.filling + 0.8)
            
            print 1*self.band_factor
            
            self.band_factor = 0.24870941840919195/1
            
#            self.band_factor = ((6.6/360.) *( 2*np.pi) +   ((self.Q*self.filling/100+self.filling/100)*np.pi/6.)) 
 
            self.band_factor = ((7/360.) *( 2*np.pi) +   ((self.filling/100)*np.pi/6.)) 

   
            print "Band standard deviation", self.band_factor *10, "degrees"
            angles1 = np.random.uniform(0,np.pi,len(r_spots))  # random angular position in longintude
            angles2 = np.random.normal(np.pi/2.,1*self.band_factor, len(r_spots)) # equatorial qith 10 degree sd
    #        angles2 = np.random.normal(np.pi/4.,0.1, len(r_spots)) # random angular position in latitutude
        elif SpatDist == 'Longitudinal':
            print "spatial distribution is LONGITUDINAL"     
            self.groups = 4
#            self.band_factor  = 0.02*self.filling/(self.groups/2.) + 0.8
            self.band_factor  = 1.1*(0.02*self.filling + 0.8)
            
            
            self.band_factor = ((7/360.) *( 2*np.pi) +   ((self.filling/100)*np.pi/6.)) 

            
            
            self.cen1 = np.random.uniform(0,np.pi)

            angles1_0=[]
            angles2_0=[]
            r_spots0=[]
            for q in range(self.groups):
                cen = self.cen1 + q*(2*np.pi)/self.groups
                print "CEN", cen, q       
                angles1 = np.random.normal(cen, 1*self.band_factor, len(r_spots))
                angles1 = np.where(angles1 <= 2*np.pi, angles1, angles1-2*np.pi)
                angles1 = np.where(angles1 >=0, angles1, 2*np.pi-angles1)
                angles2 = np.random.normal(np.pi/2.,1*self.band_factor, len(r_spots)) # equatorial qith 10 degree sd        
                idx = np.argwhere(angles1<=np.pi)
                angles1_0 = angles1_0 + angles1[idx].tolist()
                angles2_0 = angles2_0 + angles2[idx].tolist()
                r_spots0 = r_spots0 + r_spots[idx].tolist()  
            r_spots = np.array(r_spots0)
            angles1 = np.array(angles1_0)
            angles2 = np.array(angles2_0)
            
            print "R_SPOTS", len(r_spots)
      
      
        x =  np.cos(angles1)   # translate angles into x position on star disc (relative to star radius)
        y =  np.cos(angles2)  # translate angles into y position on star disc (relative to star radius)
    
    #    spot_positions = [{'x':0, 'y':0, 'r':0.0000001}, {'x':0.8, 'y':0.0, 'r':0.005}]
    #    spot_positions = [{'x':0, 'y':0, 'r':0.000001}]
        spot_positions =[]   # list of spot coordinates and relative radii per spot
        for i in range (len(r_spots)):
            spot_positions.append({'x':x[i], 'y':y[i], 'r':r_spots[i]})
        
        self.spot_positions = spot_positions
      
      else: 
        self.spot_positions =[]

                        
      return self.spot_positions    


  def calcStellarModel(self):
    print "initial filling factor and Q", 100*np.sum(self.spots_mask)/ np.sum(self.star_mask) , np.sum(self.facs_mask)/ np.sum(self.spots_mask)
 
    if self.useSpotfile == 0:   
        
        print 100*np.sum(self.spots_mask)/ np.sum(self.star_mask) , np.sum(self.facs_mask)/ np.sum(self.spots_mask)
  
        self.extra_spots_flash(0)
        
        print 100*np.sum(self.spots_mask)/ np.sum(self.star_mask) , np.sum(self.facs_mask)/ np.sum(self.spots_mask)        
        
        self.extra_spots_top_up(0)
        
        print 100*np.sum(self.spots_mask)/ np.sum(self.star_mask) , np.sum(self.facs_mask)/ np.sum(self.spots_mask)        
        
        self.extra_spots_flash(1)
        print 100*np.sum(self.spots_mask)/ np.sum(self.star_mask) , np.sum(self.facs_mask)/ np.sum(self.spots_mask)        
        
        self.extra_spots_top_up(1)
        
        print 100*np.sum(self.spots_mask)/ np.sum(self.star_mask) , np.sum(self.facs_mask)/ np.sum(self.spots_mask)        
        
        self.FF = 100*np.sum(self.spots_mask)/ np.sum(self.star_mask)
      
        print "final filling factor and Q", 100*np.sum(self.spots_mask)/ np.sum(self.star_mask) , np.sum(self.facs_mask)/ np.sum(self.spots_mask)
        print "required filling factor and Q", self.filling , self.Q

#    plt.figure('spots and facs')
#    plt.imshow((self.facs_mask*20) + (self.spots_mask*10) + self.star_mask) 
#    
#    
#    xxxx
       
    return self.spots_mask, self.facs_mask
    
     
  
  def calcLC(self, wl, cr, c1, c2, xArr):

    RplanetPx = self.RstarPx * cr**.5
     
    self.calcStellarModel() 
    star0 = self.star_mask
    spots = self.spots_mask
    facs = self.facs_mask
    FF = np.sum(self.spots_mask)/ np.sum(self.star_mask)
    spots *=star0 #trims the overhangs
    facs *=star0  
    
    Tspot = self.Tspot
    Tfac  = self.Tfac
    Tstar = self.Tstar
    
 
    if self.addFacs == 1:    
        star = star0 - spots - facs
        a = (star/star.max())*Tstar+ (spots/spots.max())*Tspot + (facs/facs.max())*Tfac
        a = (star/star.max())*60+ (spots/spots.max())*10 + (facs/facs.max())*80

    elif self.addFacs ==0 :
        star = star0 - spots
        a = (star/star.max())*Tstar+ (spots/spots.max())*Tspot
        a = (star/star.max())*60+ (spots/spots.max())*10 


    # option to use if only one limb darkening coeff used for all wavelengths or zero
    limb = self.limb_darkening(I0 = 1, c1 = c1[len(c1)/2], c2 = c2[len(c1)/2], Rimg = self.ImHalfSizePx, pix = self.RstarPx)
    limb = 2-limb #inverts to get dimming on image at edge
    star_dem  = a*limb
    
    plt.figure(23)   
    plt.imshow(star_dem,vmin= 0, vmax=100,origin='lower', interpolation = 'none')
  
    
    lc_star = np.zeros_like(xArr)
    lc_spots = np.zeros_like(xArr)
    lc_facs =  np.zeros_like(xArr)
    lc_star0 = np.zeros_like(xArr)
         
    print 'start'
    a = time.clock()
    for i in range(len(xArr)):
      xPos = xArr[i]
      yPos = self.ImHalfSizePx-0.5 + self.b*self.RstarPx
      P_pos = (xPos, yPos)
      P_mask = CircAp(P_pos, RplanetPx)
       
      planet_star = ApPhot(star, P_mask)[0][3]
      lc_star[i] = np.sum(star) - planet_star
      
      planet_spots = ApPhot(spots, P_mask)[0][3]
      lc_spots[i] = np.sum(spots) - planet_spots
           
      planet_star0 = ApPhot(star0, P_mask)[0][3]
      lc_star0[i] = np.sum(star0) - planet_star0
      
      planet_facs = ApPhot(facs, P_mask)[0][3]
      lc_facs[i] = np.sum(facs) - planet_facs
      
      apertures = P_mask
 
      plt.figure(23)
      apertures.plot(color='0.25', lw=1.5, alpha=0.3)  
     
    lc_array = np.ones((len(wl), len(xArr)))
    lc_star = lc_star*lc_array
    lc_spots = lc_spots*lc_array
    lc_star0 = lc_star0 * lc_array
    lc_facs = lc_facs * lc_array
    
    if self.usePhoenix ==True:
        pass
    else: 
        bb_star = np.pi*self.planck(wl, self.Tstar)
        Tspot = self.Tspot
        Tfac  = self.Tfac
        bb_sp = np.pi*self.planck(wl, Tspot)
        bb_fa = np.pi*self.planck(wl, Tfac)  
            
        lc_star *= bb_star.reshape(len(bb_star),1)
        lc_spots *= bb_sp.reshape(len(bb_sp),1)
        lc_facs *= bb_fa.reshape(len(bb_fa),1)
        lc_star0 *= bb_star.reshape(len(bb_star),1)
           
        sim_ideal = star0*bb_star[0]
        sim = star*bb_star[0] + spots*bb_sp[0] +  facs*bb_fa[0]

    if self.addFacs == 1:
        lc = lc_star+lc_spots+lc_facs
    elif self.addFacs == 0: 
        lc = lc_star+lc_spots

    lc_i = lc_star0

    b = time.clock()
    print "time", b-a
  
    
    # option to add limb darkening here for diff wavlengths (c1 and c2 now have to be arrays)
    m = pytransit.MandelAgol(eclipse=False)
    y4z = self.b*self.RstarPx
    x4z = xArr-(self.ImHalfSizePx-0.5)
    z = ((x4z**2+y4z**2)**0.5)/(self.RstarPx)
   
    for i in range(lc.shape[0]):
        lc_nolimb = m(z, np.sqrt(cr), [0,0]) 
        lc_limb = m(z, np.sqrt(cr), [c1[i],c2[i]])

        lc[i] *= lc_limb/lc_nolimb
        lc_i[i] *= lc_limb/lc_nolimb        
                
#         don't normalise here otherwise OOT change in count does not appear       
#        lc[i] /= lc[i][0]
#        lc_i[i] /= lc_i[i][0]

    return wl, lc, lc_i, sim, sim_ideal, FF
    
