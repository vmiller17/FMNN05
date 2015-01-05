# -*- coding: utf-8 -*-
from scipy import *
import numpy as np
from assimulo.solvers import IDA
from assimulo.problem import Implicit_Problem
	
		
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
        
def init_woodpecker():
    
    y = np.array([1,0,0,0,0,0,0,0])
    yp = np.array([0,0,0,0,0,-g,0,0])
    
    return y,yp
    
def woodpecker(t, y, yp, sw):
	"""
	Residual function of the 7-bar mechanism in
	Hairer, Vol. II, p. 533 ff, see also formula (7.11)
	written in residual form
	y,yp vector of dim 20, t scalar
	
	y = [z phiS phiB lamb1 lamb2 zp phiSp phiBp]
	yp = [zp phiSp phiBp lamb1p lamb2p zpp phiSpp phiBpp]
	"""

        
        
        
        # STATE I
        #y = [z phiS phiB lamb1 lamb2 zp phiSp phiBp]
	#yp = [zp phiSp phiBp lamb1p lamb2p zpp phiSpp phiBpp]
	#      0   1     2      3       4    5    6     7  
	if sw[0]:
            res_0 = y[5:7] - yp[0:2] #check index        
            res_1 = -(mS+mB)*yp[5] - mB*lS*yp[6] - mB*lG*yp[7] - (mS+mB)*g
            res_2 = -(mB*lS)*yp[5] - (JS+mB*lS**2)*yp[6] - (mB*lS*lG)*yp[7] + cp*(y[7]-y[6])-mB*lS*g-y[3]
            res_3 = -mB*lG*yp[5] - (mB*lS*lG)*yp[6] - (JB+mB*lG**2)*yp[7] + cp*(y[6]-y[7])-mB*lG*g-y[4]
            res_4 = y[3]
            res_5 = y[4]
	
	
	# STATE II
	#y = [z phiS phiB lamb1 lamb2 zp phiSp phiBp]
	#yp = [zp phiSp phiBp lamb1p lamb2p zpp phiSpp phiBpp]
	#      0   1     2      3       4    5    6     7  
	if sw[1]:
   	    res_0 = y[5:7] - yp[0:2] #check index        
            res_1 = -(mS+mB)*yp[5] - mB*lS*yp[6] - mB*lG*yp[7] - (mS+mB)*g - y[4]
            res_2 = -(mB*lS)*yp[5] - (JS+mB*lS**2)*yp[6] - (mB*lS*lG)*yp[7] + cp*(y[7]-y[6])-mB*lS*g-hS*y[3]-rS*y[4]
            res_3 = -mB*lG*yp[5] - (mB*lS*lG)*yp[6] - (JB+mB*lG**2)*yp[7] + cp*(y[6]-y[7])-mB*lG*g
            res_4 = (rS-r0) + hS*y[1]
            res_5 = yp[0] + rS*yp[1]
	
	
	
	# STATE III
	#y = [z phiS phiB lamb1 lamb2 zp phiSp phiBp]
	#yp = [zp phiSp phiBp lamb1p lamb2p zpp phiSpp phiBpp]
	#      0   1     2      3       4    5    6     7 
	if sw[2]: 
   	    res_0 = y[5:7] - yp[0:2] #check index        
            res_1 = -(mS+mB)*yp[5] - mB*lS*yp[6] - mB*lG*yp[7] - (mS+mB)*g - y[4]
            res_2 = -(mB*lS)*yp[5] - (JS+mB*lS**2)*yp[6] - (mB*lS*lG)*yp[7] + cp*(y[7]-y[6])-mB*lS*g-hS*y[3]-rS*y[4]
            res_3 = -mB*lG*yp[5] - (mB*lS*lG)*yp[6] - (JB+mB*lG**2)*yp[7] + cp*(y[6]-y[7])-mB*lG*g
            res_4 = (rS-r0) - hS*y[1]
            res_5 = yp[0] + rS*yp[1]
            
            
        # STATE IV
        # Switch sign
        if sw[3]:
            yp[1] = - yp[1]
            
	
	
	
        return hstack((res_0,res_1,res_2,res_3,res_4,res_5))
	
def state_events(t,y,yp,sw):
    
	e1 = -hS*y[1] - (rS-r0)
	
	e2 = -hS*y[1] + (rS-r0)
	
	e3 = y[3]
	
	e4 = -hB*y[2] +lS+lG-lB-r0
	
	return np.array([e1, e2 , e3 , e4])
	
def handle_event(solver,event_info):
    
    state_info = event_info[0]
    
    # State I => State II 
    if state_info[0] != 0 and solver.sw[0]:
        solver.sw[0] = not solver.sw[0]
        solver.sw[1] = not solver.sw[1]
        
    # State I => State III
    if state_info[0] != 0 and solver.sw[1]:
        solver.sw[0] = not solver.sw[0]
        solver.sw[2] = not solver.sw[2]
        
    # State II => State I
    if state_info[1] != 0 and solver.sw[1]:
        solver.sw[0] = not solver.sw[0]
        solver.sw[2] = not solver.sw[2]
        
    # State III => State I
    if state_info[2] != 0 and solver.sw[2]:
        solver.sw[0] = not solver.sw[0]
        solver.sw[2] = not solver.sw[2] 
        
    # State III => State IV
    if state_info[3] != 0 and solver.sw[2]:
        solver.sw[2] = not solver.sw[2]
        solver.sw[3] = not solver.sw[3]
        
    # State IV => State III
    if solver.sw[3]:
        solver.sw[3] = not solver.sw[3]
        solver.sw[2] = not solver.sw[2]          
        
        
        
        
y0,yp0 = init_woodpecker()
t0 = 0.0
sw0 = np.array([1,0,0,0])

model = Implicit_Problem(woodpecker,y0,yp0,t0,sw0)

model.state_events = state_events
model.handle_event = handle_event

solver = IDA(model)
solver.simulate(3)





        

        
	
	
	 
	