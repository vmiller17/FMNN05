# -*- coding: utf-8 -*-
from scipy import *
import numpy as np
	
	
def init_woodpecker():
    
    y = np.array([0,0,0])
    ypp = np.array([0,0,0])
    
    return y,yp
    
def woodpecker(t, y, ypp):
	"""
	Residual function of the 7-bar mechanism in
	Hairer, Vol. II, p. 533 ff, see also formula (7.11)
	written in residual form
	y,yp vector of dim 20, t scalar
	
	y = [z phiS phiB]
	ypp = [zpp phiSpp phiBpp]
	"""
	
	# -*- coding: utf-8 -*-
        # Mechanical and geometric constants of the woodpecker model
        # according to 
        # Ch. Glocker: Dynamik von StarrkÃ¶rpersystemen mit Reibung und StÃ¶ssen
        # PhD thesis TU MÃ¼nchen, July 1995, pp. 162ff
        #
        mS = 3.0e-4 # Mass of sleeve [kg]
        JS = 5.0e-9 # Moment of inertia of the sleeve [kgm]
        mB = 4.5e-3 # Mass of bird [kg]
        masstotal=mS+mB # total mass
        JB = 7.0e-7 # Moment of inertia of bird [kgm]
	r0 = 2.5e-3 # Radius of the bar [m]
        rS = 3.1e-3 # Inner Radius of sleeve [m]
        hS = 5.8e-3 # 1/2 height of sleeve [m]
        lS = 1.0e-2 # verical distance sleeve origin to spring origin [m]
        lG = 1.5e-2 # vertical distance spring origin to bird origin [m]
        hB = 2.0e-2 # y coordinate beak (in bird coordinate system) [m]
        lB = 2.01e-2 # -x coordinate beak (in bird coordinate system) [m]
        cp = 5.6e-3 # rotational spring constant [N/rad]
        g  = 9.81 #  [m/s^2]
        
        
        
        res_1 = -(mS+mB)*zpp - mB*lS*phiSpp - mB*lG*phiBpp - (mS+mB)*g
        res_2 = -(mB*lS)*zpp - (JS+mB*lS**2)*phiSpp - (mB*lS*lG)*phiBpp + cp(phiS-phiB)-mB*lS*g-lamb[1]
        res_3 = -mB*lG*zpp - (mB*lS*lG)*phiSpp - (JB+mB*lG**2)*phiBpp + cp(phiS-phiB)-mB*lS*g-lamb[2]
        res_4 = lamb[1]
        res_5 = lamb[2]
        
	
	return hstack((res_1,res_2,res_3,res_4,res_5))
	
	
	
	
	
	
	
	
	
	
	