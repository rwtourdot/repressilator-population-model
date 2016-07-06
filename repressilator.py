#!/usr/bin/env python
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.integrate import odeint

beta = 0.5
alpha0 = 0 #0
alpha = 100
n = 2
# declare model
def mymodel(var,t):
    m = var[:3]
    p = var[3:]
    dm0 = - m[0] + alpha/(1+p[2]**n) + alpha0
    dm1 = - m[1] + alpha/(1+p[0]**n) + alpha0
    dm2 = - m[2] + alpha/(1+p[1]**n) + alpha0
    dp0 = - beta*(p[0] - m [0])
    dp1 = - beta*(p[1] - m [1])
    dp2 = - beta*(p[2] - m [2])
    return [dm0,dm1,dm2,dp0,dp1,dp2]

time = np.linspace(0.0,100.0,1000)
minit = np.array([1.0,4.0,1.0,2.0,1.0,1.0])
y = odeint(mymodel,minit,time)

plt.plot(time,y[:,0],time,y[:,1],time,y[:,2])
plt.xlabel('t')
plt.ylabel('y')
plt.show()
