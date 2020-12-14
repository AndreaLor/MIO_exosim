# -*- coding: utf-8 -*-
"""
Created on Sat Dec 24 12:50:40 2016

@author: c1341133
"""
import numpy as np
import matplotlib.pyplot as plt
from photutils import CircularAperture as CircAp
from photutils import aperture_photometry as ApPhot
import scipy.constants as cnst
from scipy.misc import imrotate
from scipy.ndimage.interpolation import shift
from scipy import interpolate
from scipy import signal
import numpy.random as rnd
import time
import multiprocessing as mp



#==============================================================================
# Primary Transit Class:
#==============================================================================



class SpotSim():
  def __init__(self, spotsList, spotParams, Tstar, b, Frame_pix, pix, useSpotfile, realization):
    """
    PTC: Primary Transit Class
    
    :INPUTS:
    spotsList  -- List containing dicrtiories of spot radius/position (deg)
    spotPrams  -- Dictionary containing:
      CONTspot   -- spot-to-photosphere temperature contrast /kelvin      (scalar)
      CONTfac    -- facula-to-photosphere temperature contrast /kelvin    (scalar)
      Q          -- facula-to-spot area ratio                             (scalar)
    Tstar      -- stellar photosphere temperature /kelvin               (scalar)
    a          -- orbital semi-major axis /metres                       (scalar)
    b          -- planetary impact parameter /stellar radii             (scalar)
    pix        -- number of pixels per star radius                      (integer)
    
    :NOTES:
    v1.0. Built  for implementation in ExoSim stellar variability module (Sarkar 
    & Pascale, 2015) and for independent investigation into stellar variability 
    as part of PX4310 dissertation 'TWINKLE: a British space mission to explore
    faraway worlds'. Luke Johnson (C1216542). Last edit: 17 April 2016.
    
   """
    
    self.Tstar = Tstar
    self.RstarPx = pix
#    self.ImHalfSizePx = pix*3/2
    self.ImHalfSizePx = Frame_pix
    self.spotsList = spotsList
    self.spotParams = spotParams
    self.b = b
    self.useSpotfile = useSpotfile
    self.realization = realization
    
    
    self.spotsLayerMask = None
    self.spotMasksList = None
    self.initSpotsMask()
   
    
    
  def planck(self, w1, T):
    a = 1.191042768e8    #*pq.um**5 *pq.W/ pq.m**2 /pq.sr/pq.um
    b = 14387.7516       #*1*pq.um * 1*pq.K
    x = b/(w1*T)
    bb = a/w1**5 / (np.exp(x) - 1.0)
    return bb
  
  def calcSpotMask(self, Rspot, Xspot, Yspot, RstarPx, ImHalfSizePx, Q):
    R  = RstarPx # radius of star in pixels
    r = np.sqrt(Xspot**2 + Yspot**2) # radial distance of spot from centre in pixels
    a = 1.0   
    b = a*np.sin(np.arccos(r/R)) # change in width of ellipse due to projection
   
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
   
    y_s = np.linspace(-ImHalfSizePx, ImHalfSizePx, 2*ImHalfSizePx+1) #set up ellipse
    x_s = np.linspace(-ImHalfSizePx, ImHalfSizePx, 2*ImHalfSizePx+1) 
    xx_s, yy_s = np.meshgrid(x_s, y_s)
    arg_s = np.sqrt((xx_s/a)**2 + (yy_s/b)**2)
    R_fac = Rspot*np.sqrt(Q + 1)
    
    spot = np.where(arg_s> Rspot, 0.0, 1.0)
    
    fac = np.where((arg_s> R_fac), 0.0, 2.0)
    spot0 = imrotate(spot, theta,interp='bilinear') # rotate spot and facula 
    fac0 = imrotate(fac, theta,interp='bilinear')               
    spot = np.where(spot0==0, spot0, spot.max())# account for artifact from rotating
    fac = np.where(fac0==0, fac0, fac.max())
    
#    print "spot", spot.max(), Rspot, Xspot, Yspot, RstarPx, ImHalfSizePx, Q, r
    
    return {'spot':spot, 'fac':fac}
  
  
  def star(self, Num_phot, Rimg, pix):
    F_star = Num_phot
    x = np.linspace(-Rimg, Rimg, 2*Rimg+1, dtype=np.float)
    yy, xx = np.meshgrid(x, x)
    arg = np.sqrt(yy**2 + xx**2)
    arg[arg.shape[0]/2, arg.shape[1]/2] = 1
    arg = np.where(arg< pix, arg, 0)
    arg = np.where(arg<= 0, arg, 1)
    argI = arg*F_star
    return argI
  
  def calcSpotsLayerMask(self, spotsList, spotMaskList):
      if self.useSpotfile == False:
        spotsLayerFac  = np.ones( (2*self.ImHalfSizePx+1, 2*self.ImHalfSizePx+1))
        spotsLayerSpot = np.ones( (2*self.ImHalfSizePx+1, 2*self.ImHalfSizePx+1))
        spotsLayerMask  = np.zeros( (2*self.ImHalfSizePx+1, 2*self.ImHalfSizePx+1), dtype = int)
        for i in range(len(spotsList)):
          spotPosPx = {
              'x':np.round(spotsList[i]['x']*self.RstarPx),
              'y':np.round(spotsList[i]['y']*self.RstarPx)
            }
          # Convert Masks to 0 where spot, 1 else and shift. 
          # This can be done properly in the createMasks(), keeping it refactoring purposes now. 
          spotsLayerSpot *= shift( np.where(spotMaskList[i]['spot']>0, -1, spotMaskList[i]['spot'])+1, 
                                  (spotPosPx['y'], spotPosPx['x']), cval = 1)
          spotsLayerFac  *= shift( np.where(spotMaskList[i]['fac']>0, -1, spotMaskList[i]['fac'])+1, 
                                  (spotPosPx['y'], spotPosPx['x']), cval = 1)
          
          # Following rounding very important
          spotsLayerSpot = np.round(spotsLayerSpot)
          spotsLayerFac = np.round(spotsLayerFac)
          
        spotsLayerMask = np.where(spotsLayerFac <1, 1, spotsLayerMask)
        spotsLayerMask = np.where(spotsLayerSpot<1, 2, spotsLayerMask)
        fac = np.where(spotsLayerMask==1,1,0)
        spots = np.where(spotsLayerMask==2,1,0)
      else:
          print "using file spots"      
          a = np.load('/Users/c1341133/Desktop/spot_array.npy')
          b = np.load('/Users/c1341133/Desktop/fac_array.npy')
          spots = a[...,self.realization]
          fac = b[..., self.realization]
  
      return fac, spots

  def limb_darkening(self, I0, c1, c2, Rimg, pix):
      x = np.linspace(-1, 1, 2*Rimg+1)*Rimg/pix
      y = np.linspace(-1, 1, 2*Rimg+1)*Rimg/pix
      xx, yy = np.meshgrid(x, y)
      arg = np.sqrt(xx**2 + yy**2)
      u = np.sqrt(1 - arg**2)
      arg_l = I0 * ( 1 - c1*(1-u) - c2*(1-u)**2)
      arg_l = np.where(np.isnan(arg_l), 0, arg_l)
      return arg_l

  def initSpotsMask(self):
      if self.useSpotfile == False:
          self.spotMasksList = []
          for i in range(len(self.spotsList)):
              self.spotMasksList.append(self.calcSpotMask( self.spotsList[i]['r']*
                    self.RstarPx, self.spotsList[i]['x']*self.RstarPx, 
                    self.spotsList[i]['y']*self.RstarPx,
                    self.RstarPx, self.ImHalfSizePx, self.spotParams['Q']))
          self.spotsLayerMask = self.calcSpotsLayerMask(self.spotsList, self.spotMasksList)
      else:
          self.spotsLayerMask = self.calcSpotsLayerMask(self.spotsList, self.spotMasksList)
          

    

  def calcStellarModel(self, wl, c1, c2, useP):
    if useP == False:
        bb_star = self.planck(wl, self.Tstar)
        Tspot = self.Tstar + self.spotParams['CONTspot']
        Tfac  = self.Tstar + self.spotParams['CONTfac']
    
        bb_sp = self.planck(wl, Tspot)
        bb_fa = self.planck(wl, Tfac)
        
        star = self.star(bb_star, self.ImHalfSizePx, self.RstarPx)
        spots = self.spotsLayerMask[1]*(bb_sp-bb_star)
        
        FF = np.sum((spots/(bb_sp-bb_star)))/np.sum(star/(bb_star))
      
        facs = self.spotsLayerMask[0]*(bb_fa-bb_star)
        
        limb = self.limb_darkening(I0 = 1, c1 = c1, c2 = c2, Rimg = self.ImHalfSizePx, pix = self.RstarPx)
        sim = limb*(star+spots+facs)
        sim_IDEAL = limb*star
        
    elif useP == True:
        star = self.star(1, self.ImHalfSizePx, self.RstarPx)
        spots = self.spotsLayerMask[1]*(100-1)
        facs = self.spotsLayerMask[0]*(1000-1)
    
        sim = limb*(star+spots+facs)
        sim_IDEAL = limb*star
        
        #apply phoenix fluxes after transit
    
    
    return sim, sim_IDEAL, FF
    
  def store_spots(self):
      return self.spotsLayerMask[1],self.spotsLayerMask[0]
  
  def calcLC(self, wl, cr, c1, c2, xArr, useP, output):
#    zArr_sign = np.zeros_like(zArr)
#    zArr_sign[:-1] = np.sign(zArr[1:] - zArr[:-1])
#    zArr_sign[-1] = zArr_sign[-2]

    RplanetPx = self.RstarPx * cr**.5
    
    sim, sim_IDEAL, FF = self.calcStellarModel(wl, c1, c2, useP)    
    star_spotted = np.sum(sim)
    star_IDEAL = np.sum(sim_IDEAL)
          
    lc = np.zeros_like(xArr)
    lc_i = np.zeros_like(xArr)
    a =time.clock()
    print 'start'
    for i in range(len(xArr)):
      xPos = xArr[i]
      yPos = self.ImHalfSizePx + self.b*self.RstarPx
      P_pos = (xPos, yPos)
      P_mask = CircAp(P_pos, RplanetPx)
      planet = ApPhot(sim, P_mask)[0][3]
      planet_I = ApPhot(sim_IDEAL, P_mask)[0][3]
      lc[i] = star_spotted - planet
      lc_i[i] = star_IDEAL - planet_I
        
#        apertures = P_mask
#        plt.figure('aperture')
#        plt.imshow(sim, cmap='gray_r', origin='lower', interpolation = 'None')
#        apertures.plot(color='blue', lw=1.5, alpha=0.5)      
#    xx            
        
    b = time.clock()
    print "time", b-a
        
    lc = np.where(lc==0, np.max(lc), lc)   
    lc_i = np.where(lc_i==0, np.max(lc_i), lc_i)
    
    lc /=np.nanmax(lc)
    lc_i /=np.nanmax(lc_i)
    


#    plt.figure('Spotted Star & Correction Matrix')
#    plt.subplot(211).axhspan(self.Rimg - self.B - self.Rplanet, self.Rimg -
#self.B + self.Rplanet, alpha = 0.5, color='black')
#    plt.ylim(0, 2*self.Rimg)
       
    output.put( (wl, lc, lc_i, sim, sim_IDEAL, FF))

def size_dist(filling, SUR):
      
    A_mean = 0.57
    
    # these are from table 6 in Solanki and Unruh
    ff = np.array([0.0003,0.003,0.014,0.07,0.28,0.62])
    sigma = np.array([3.75,4.0,5.04,10.4,30.3,62.0])
    gradmax = np.array([4.5,38.5,106,106,41.2,20.8])

    r = 3 # degree 7 works well
    z = np.polyfit(ff, sigma, r)

    y = 0
    for i in range (0,r+1):
        y = y + z[i]*(filling/100.)**(r-i)
    Std = y
    
    # array of spot areas, limit max and min areas
    A = np.linspace(0,10000,100000)
    idx = np.int( np.argwhere(A>1.0)[0] )
    #idx2 =np.argwhere(A>60.0)[0]
    idx2 =np.int( np.argwhere(A>5000.0)[0])
    Grad_max = 10
    # generate list of areas
    for i in range(10000000):
        Grad = Grad_max * np.exp(-(np.log(A)-np.log(A_mean))**2/(2*np.log(Std)))  # dn/da
        Grad = Grad[idx:idx2] # limit dn/da to within max and min spot areas
        A0 = A[idx:idx2]   # limit max and min spot areas 
        N = Grad*np.gradient(A0)  # convert dn/da to N per da, where da = np.gradient of A0
        N=np.round(N*100)    # expand the N direction by 100 and obtain whole values of N
        idx0 = np.int( np.argwhere(N<1)[0])  # exclude N < 1
        N = N[0:idx0]   # as above
        A0 = A0[0:idx0]  # as above
        TA = np.sum(N*A0)  # obtain total area (*100)
        if (TA/100)*1e-6*SUR*100 < filling:  # condition - if total area as percent of hemispehere < ff
            Grad_max += 0.1    # increase maximum of dn/da - adjust the distribution maximum
        if (TA/100)*1e-6*SUR*100 >= filling:  # once filling factor achieved, record grad_max and quit
            print "Grad_max found", Grad_max
            print "Std", Std
            print "A mean", A_mean
            print "Filling factor per hemisphere percent", (TA/100)*1e-6*SUR*100 # actual filling factor
            break
    
    
    A=A0  # reset A as between max and min areas - now are bins of width np.gradient A
    
    alist = []   # set up list of spot areas
    for i in range (A.shape[0]):
        i = np.int(i)

        alist = alist+[A[i]]*np.int(N[i])  #list all spots as areas
    
    total_area = (TA/100)   #total area of umbra
    asum = 0  # accumlated sum of areas
    alist2=[]  # list of areas of selected spots
    for i in range (10000000):
        if asum > total_area:
            break
        a = np.int(np.random.uniform(0,len(alist))) # select a spot at random from all spots
        
        asum += alist[a]# add spot area to total sum
        alist2.append(alist[a]) # add spoty to list of spot areas

    alist2.sort()  # sort into ascending order the list of selected spot areas
    
#        plt.figure(222)
#        plt.hist(alist2,bins = len(A))
    
    return alist2
    

def spatial_dist(Area_spots, Star_pix, dist): 
    r_spots = ((Area_spots/np.pi)**0.5)/Star_pix # array of star spot radii relative to star
    if dist == 'uniform':
        angles1 = np.random.uniform(0,np.pi,len(r_spots))  # random angular position in longintude
        angles2 = np.random.uniform(0,np.pi,len(r_spots)) # random angular position in latitutude
    elif dist == 'equatorial':
        angles1 = np.random.uniform(0,2*np.pi,len(r_spots))  # random angular position in longintude
        angles2 = np.random.normal(0,np.pi/2.,len(r_spots)) # random angular position in latitutude

    x = Star_pix*np.cos(angles1)/Star_pix  # translate angles into x position on star disc
    y = Star_pix*np.cos(angles2)/Star_pix # translate angles into y position o
##
#    spotsList = [{'x':0.0, 'y':0.3, 'r':0.1}, {'x':0.8, 'y':0.0, 'r':0.1}]
#
    spotsList =[]   # list of spot coordinates and relative radii per spot
    for i in range (len(r_spots)):
        spotsList.append({'x':x[i], 'y':y[i], 'r':r_spots[i]})
#        
    return spotsList 