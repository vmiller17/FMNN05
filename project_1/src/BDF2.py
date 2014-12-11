from assimulo.explicit_ode import *
from assimulo.ode import *
import numpy as N
import pylab as P
import scipy.linalg as SL
from assimulo.solvers import CVode
from assimulo.support import *

class BDF_2(Explicit_ODE):
    """
    Explicit Euler.
    """
    tol=1.e-8
    maxit=100
    maxsteps=5000
    
    def __init__(self, problem):
        Explicit_ODE.__init__(self, problem) #Calls the base class
        
        #Solver options
        self.options["h"] = 0.01
        
   
    def _set_h(self,h):
            self.options["h"] = float(h)

    def _get_h(self):
        return self.options["h"]
        
    h=property(_get_h,_set_h)
        
    def integrate(self, t, y, tf, opts):
        """
        _integrates (t,y) values until t > tf
        """
        h = self.options["h"]
        h = min(h, abs(tf-t))
        
        #Lists for storing the result
        tres = []
        yres = []
        
        for i in range(self.maxsteps):
            if t >= tf:
                break
            
            if i==0:  # initial step
                t_np1,y_np1 = self.step_EE(t,y, h)
            else:   
                t_np1, y_np1 = self.step_BDF2([t,t_nm1], [y,y_nm1], h)
                
            t_nm1 = t    
            t = t_np1
            
            y_nm1 = y
            y =y_np1
            
            tres.append(t)
            yres.append(y.copy())
        
            h=min(self.h,N.abs(tf-t))
        else:
            raise Explicit_ODE_Exception('Final time not reached within maximum number of steps')
        
        return ID_PY_OK, tres, yres
    
    def step_EE(self, t, y, h):
        """
        This calculates the next step in the integration with explicit Euler.
        """
        
        f = self.problem.rhs
        return t + h, y + h*f(t, y) 
        
    def step_BDF2(self,T,Y, h):
        """
        BDF-2 with Fixed Point Iteration and Zero order predictor
        
        alpha_0*y_np1+alpha_1*y_n+alpha_2*y_nm1=h f(t_np1,y_np1)
        alpha=[3/2,-2,1/2]
        """
        alpha=[3./2.,-2.,1./2]
        f=self.problem.rhs
        
        t_n,t_nm1=T
        y_n,y_nm1=Y
        # predictor
        t_np1=t_n+h
        y_np1_i=y_n   # zero order predictor
        # corrector with fixed point iteration
        for i in range(self.maxit):
            
            y_np1_ip1=(-(alpha[1]*y_n+alpha[2]*y_nm1)+h*f(t_np1,y_np1_i))/alpha[0]
            if SL.norm(y_np1_ip1-y_np1_i) < self.tol:
                return t_np1,y_np1_ip1
            y_np1_i=y_np1_ip1
        else:
            raise Explicit_ODE_Exception('Corrector could not converge within % iterations'%i)
            
##Define the rhs     
#def f(t,y):
#
#    k = 50
#    y1 = y[2]
#    y2 = y[3]
#    y3 = -y[0]*k*(N.sqrt(y[0]**2 + y[1]**2)-1)/N.sqrt(y[0]**2 + y[1]**2)
#    y4 = -y[1]*k*(N.sqrt(y[0]**2 + y[1]**2)-1)/N.sqrt(y[0]**2 + y[1]**2)-1
#
#    return N.array([y1,y2,y3,y4])
#     
#    
#pend_mod=Explicit_Problem(f, y0=N.array([1, 1, 0, 0]))
#pend_mod.problem_name='Nonlinear Pendulum'
#
##Define an explicit solver
#exp_sim = BDF_2(pend_mod) #Create a BDF solver
#exp_sim2 = CVode(pend_mod)
#
#time = 30
#t1,y1 = exp_sim.simulate(time)
#t2,y2 = exp_sim2(time)
#
#P.plot(y1[:,0],y1[:,1],y2[:,0],y2[:,1])
#P.grid()
#P.show()
