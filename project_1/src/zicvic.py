import numpy as np
import pylab as P
from assimulo.problem import Explicit_Problem
from assimulo.solvers.sundials import CVode
from BDF2 import BDF_2
from BDF3 import BDF_3
from BDF4 import BDF_4


#Define the rhs 
def getRHS(k=1):
    def f(t,y):
        y1 = y[2]
        y2 = y[3]
        y3 = -y[0]*k*(np.sqrt(y[0]**2 + y[1]**2)-1)/np.sqrt(y[0]**2 + y[1]**2)
        y4 = -y[1]*k*(np.sqrt(y[0]**2 + y[1]**2)-1)/np.sqrt(y[0]**2 + y[1]**2)-1

        return np.array([y1,y2,y3,y4])
    return f
    
    
    
    
#pend_mod=Explicit_Problem(f, y0=np.array([1, 1, 0, 0] ))
#pend_mod.problem_name='Nonlinear Pendulum'

#Define an explicit solver
#exp_sim = CVode(pend_mod) #Create a BDF solver

#t, y = exp_sim.simulate(1)

rhs = getRHS(100)
pend_mod=Explicit_Problem(rhs, y0=np.array([2, 2, 0, 0]))
pend_mod.problem_name='Nonlinear Pendulum'

#Define an explicit solver
solver = CVode(pend_mod)

time = 10
t1,y1 = solver(time)

P.plot(y1[:,0],y1[:,1])
P.grid()
P.show()
P.title('CVode for k = 1, time = 30s')
P.xlabel('X-pos')
P.ylabel('Y-pos')


