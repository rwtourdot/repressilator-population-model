#!/usr/bin/env python
import numpy as np
import math
import sys
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import random
import matplotlib.animation as animation

# transition rates
basecelldivide = 0.3
basecelldeath = 0.05

# grid parameters
rows = 40
cols = 40

# repressilator parameters
#beta = 1
alpha0 = 0 #0
alpha = 100
n = 2

# declare model - repressilator
def repressilator(var,t,beta,q):
    m = var[:3]
    p = var[3:]
    dm0 = - m[0] + alpha/(1+p[2]**n) + alpha0
    dm1 = - m[1] + alpha/(1+p[0]**n) + alpha0
    dm2 = - m[2] + alpha/(1+p[1]**n) + alpha0
    dp0 = - beta*(p[0] - m [0])
    dp1 = - beta*(p[1] - m [1])
    dp2 = - beta*(p[2] - m [2])
    return [dm0,dm1,dm2,dp0,dp1,dp2]

# pixel class
class pixel:
    def __init__(self):
        self.numneigh = 0
	self.beta = 1
        self.neighx = []; self.neighy = []
	self.p0 = 0; self.p1 = 0; self.p2 = 0;
	self.m0 = 0; self.m1 = 0; self.m2 = 0;
	self.dp0 = []; self.dp1 = []; self.dp2 = [];
	self.dm0 = []; self.dm1 = []; self.dm2 = [];
	self.state = ['none']
	self.cell = False
    def add_neighbor(self,i,j):
        self.numneigh += 1
        self.neighx.append(i)
        self.neighy.append(j)
    def cell_division(self,delt,grid,i):
	prob = (1-math.exp(-(basecelldivide+self.dp2[i]/(alpha*2))*delt))
	if random.random() < prob:
	   print "cell division!",i
	   newcell = int(random.random()*self.numneigh)
	   ncx = self.neighx[newcell]
	   ncy = self.neighy[newcell]
	   newobj = grid[ncx][ncy]
	   if newobj.cell == False:
	      newobj.cell = True  
	      newobj.state[-1] = self.state[-1]  
	#else:
	#   print "tried to divide - failed",delt 
    def cell_death(self,delt,i):
	prob = (1-math.exp(-(basecelldeath)*delt))   #+self.dp0[i]/(alpha*2)
	if random.random() < prob:
	   print "cell death!",i
	   self.cell = False
	   #self.state = 'none'

############ initialize grid data structure ##################
grid = [[ pixel() for col in range(cols)] for row in range(rows)]

#construct neighbor list
for row in range(rows):
  for col in range(cols):
    if row != 0:
      grid[row][col].add_neighbor(row-1,col)
    if row != (rows-1):
      grid[row][col].add_neighbor(row+1,col)
    if col != 0:
      grid[row][col].add_neighbor(row,col-1)
    if col != (cols-1):
      grid[row][col].add_neighbor(row,col+1)

y = {}
##################### State A #############################
stateA = grid[30][15]
stateA.cell = True
stateA.p0 = 2; stateA.p1 = 1; stateA.p2 = 1;
stateA.m0 = 1; stateA.m1 = 1; stateA.m2 = 1;
stateA.state = ['A']
stateA.beta = 0.5
y[stateA.state[0]] = []

##################### State B #############################
stateB = grid[20][24]
stateB.cell = True
stateB.p0 = 2; stateB.p1 = 3; stateB.p2 = 1;
stateB.m0 = 1; stateB.m1 = 1; stateB.m2 = 1;
stateB.state = ['B']
stateB.beta = 1
y[stateB.state[0]] = []

##################### State C #############################
stateC = grid[10][15]
stateC.cell = True
stateC.p0 = 1; stateC.p1 = 2; stateC.p2 = 3;
stateC.m0 = 1; stateC.m1 = 1; stateC.m2 = 1;
stateC.state = ['C']
stateC.beta = 0.7
y[stateC.state[0]] = []

##############################################################
tfinal = 400
tsample = 8000
intervals = tsample/tfinal

for row in range(rows):
  for col in range(cols):
    obj = grid[row][col]
    if obj.cell == True:
	time = np.linspace(0.0,tfinal,tsample)   #/(row+1)
    	minit = np.array([obj.m0,obj.m1,obj.m2,obj.p0,obj.p1,obj.p2])				#np.array([1.0,1.0,1.0,2.0,1.0,1.0])
        y[obj.state[0]] = odeint(repressilator,minit,time,args=(obj.beta,1))  #[obj.dm0,obj.dm1,obj.dm2,obj.dp0,obj.dp1,obj.dp2]
    	#obj.dm0 = y[:,0];  obj.dm1 = y[:,1];  obj.dm2 = y[:,2]; obj.dp0 = y[:,3];  obj.dp1 = y[:,4];  obj.dp2 = y[:,5]

markov_steps = rows*cols
numcells = rows*cols
delt = float(tfinal)/float(tsample)
for i in range(tsample):
  for row in range(rows):
    for col in range(cols):
      obj = grid[row][col]
      if obj.cell == True:
	ytemp = y[obj.state[i]]
        obj.dm0.append(ytemp[i,0]); obj.dm1.append(ytemp[i,1]); obj.dm2.append(ytemp[i,2]);
        obj.dp0.append(ytemp[i,3]); obj.dp1.append(ytemp[i,4]); obj.dp2.append(ytemp[i,5]);  
	obj.state.append(obj.state[-1])
      elif obj.cell == False: 
	obj.state.append('none')
        obj.dm0.append(0.0); obj.dm1.append(0.0); obj.dm2.append(0.0);
        obj.dp0.append(0.0); obj.dp1.append(0.0); obj.dp2.append(0.0);
  for m in range(markov_steps):	
    rrow = random.randrange(0,rows)
    rcol = random.randrange(0,cols)
    obj = grid[rrow][rcol]
    if obj.cell == True:
	obj.cell_division(delt,grid,i)
	obj.cell_death(delt,i)

####################### animation section
typecolor = {}
typecolor['none'] = 0,0,0
typecolor['A'] = 1,0,0
typecolor['B'] = 0,1,0
typecolor['C'] = 0,0,1

fig = plt.figure(1)
toshow = [[ typecolor[grid[row][col].state[0]] for col in range(cols)] for row in range(rows)]
im = plt.imshow(toshow,interpolation='none',animated=True,vmin = 0,vmax = alpha,cmap = 'summer')  #'coolwarm'

m = 0
def updatefig(*args):
  global m
  toshow = [[ typecolor[grid[row][col].state[m]] for col in range(cols)] for row in range(rows)]
  #toshow2 = [[ grid[row][col].dm0[m] for col in range(cols)] for row in range(rows)]
  im.set_array(toshow)
  m += 1
  return im,

anim = animation.FuncAnimation(fig,updatefig,interval=10,blit = False,frames=(tsample-1),repeat=False)
writer = animation.writers['ffmpeg'](fps=80)
plt.xlabel('x')
plt.ylabel('y')
anim.save('dynamic_image.mp4',writer=writer,dpi=200)



fig2 = plt.figure(2)
toshow2 = [[ grid[row][col].dm0[0] for col in range(cols)] for row in range(rows)]
im2 = plt.imshow(toshow2,interpolation='none',animated=True,vmin = 0,vmax = alpha,cmap = 'summer')
plt.colorbar() 

m = 0
def updatefig2(*args):
  global m
  #toshow = [[ typecolor[grid[row][col].state[m]] for col in range(cols)] for row in range(rows)]
  toshow2 = [[ grid[row][col].dm0[m] for col in range(cols)] for row in range(rows)]
  im2.set_array(toshow2)
  m += 1
  return im2,

anim = animation.FuncAnimation(fig2,updatefig2,interval=10,blit = False,frames=(tsample-1),repeat=False)
writer = animation.writers['ffmpeg'](fps=80)
plt.xlabel('x')
plt.ylabel('y')
anim.save('dynamic_image2.mp4',writer=writer,dpi=200)
#plt.plot(time,obj.dp0,time,obj.dp1,time,obj.dp2)

#plt.show()
