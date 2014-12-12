import numpay as np
import pylab as P
from assimulo.problem import Implicit_Problem
from assimulo.solvers import IDA

y0 = np.array([0 0 0 0 0 0 0])
yp0 = np.array([0 0 0 0 0 0 0])
t0 = 0.0



model = Implicit_Problem(squeezer, y0, yp0, t0)
model.name = 'Squeezer'

solver = IDA(model)


