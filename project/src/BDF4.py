from assimulo.explicit_ode import *
from assimulo.ode import *
import numpy as N
import pylab as P
import scipy.linalg as SL
from assimulo.solvers import CVode

class BDF_4(Explicit_ODE):
    """
    Explicit Euler.
    """
    tol=1.e-8
    maxit=100
    maxsteps=500
    
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
            
            if i == 0:  # initial step
                t_np1,y_np1 = self.step_EE(t,y,h)
                t_nm1,y_nm1 = t_np1,y_np1
                t_nm2,y_nm2 = t_nm1,y_nm1
            else:   
                t_np1, y_np1 = self.step_BDF4([t,t_nm1,t_nm2,t_nm3],[y,y_nm1,y_nm2,y_nm3],h)
  
           
            t_nm3 = t_nm2
            t_nm2 = t_nm1
            t_nm1 = t
            t = t_np1
            
            y_nm3 = y_nm2
            y_nm2 = y_nm1
            y_nm1 = y
            y = y_np1
            
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
        
    def step_BDF4(self,T,Y, h):
        """
        BDF-4 with Fixed Point Iteration and Zero order predictor
        
        alpha_0*y_np1+alpha_1*y_n+alpha_2*y_nm1+alpha_3*y_nm2+alpha_4*y_nm3=h f(t_np1,y_np1)
        alpha=[25/12,-4,3,-4/3,1/4]
        """
        alpha=[25./12,-4.,3.,-4./3,1./4]
        f=self.problem.rhs
        
        t_n,t_nm1,t_nm2,t_nm3=T
        y_n,y_nm1,y_nm2,y_nm3=Y
        # predictor
        t_np1=t_n+h
        y_np1_i=y_n   # zero order predictor
        # corrector with fixed point iteration
        for i in range(self.maxit):            
            y_np1_ip1=(-(alpha[1]*y_n+alpha[2]*y_nm1+alpha[3]*y_nm2+alpha[4]*y_nm3)+h*f(t_np1,y_np1_i))/alpha[0]
            if SL.norm(y_np1_ip1-y_np1_i) < self.tol:
                return t_np1,y_np1_ip1
            y_np1_i=y_np1_ip1
        else:
            raise Explicit_ODE_Exception('Corrector could not converge within % iterations'%i)
            
            
#Define the rhs
def f(t,y):
    ydot = -y[0]
    return N.array([ydot])
    
#Define an Assimulo problem
exp_mod = Explicit_Problem(f, 4)
exp_mod.problem_name = 'Simple BDF-4 Example'

def pend(t,y):
    #g=9.81    l=0.7134354980239037
    gl=13.7503671
    return N.array([y[1],-gl*N.sin(y[0])])
    
pend_mod=Explicit_Problem(pend, y0=N.array([2.*N.pi,1.]))
pend_mod.problem_name='Nonlinear Pendulum'

#Define an explicit solver
exp_sim = BDF_4(pend_mod) #Create a BDF solver

t, y = exp_sim.simulate(1)

P.plot(t,y)
P.grid()
P.show()