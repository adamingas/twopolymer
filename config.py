import numpy as np
import math

noise = "yes" #yes or no
steps = 1000
runs = 100
constant = [1,5,10] #np.linspace(0,10,21).tolist()
ldensity = [100]
time_step = 0.1
#Particles a b c d with initial positions:




p1 = np.array([0,-0.05,0])
p2 = np.array([0,0.05,0])
def randomiser(l):
	radius = np.random.uniform()*l
	costh = np.random.uniform(-1, 1, size=2)
	fi = np.random.uniform(size=2) * 2 * math.pi
	sinth = np.sqrt(1 - costh ** 2)

	com1 = (p1 + p2) / 2
	radius2 = np.linalg.norm(p2 - p1)

	particle3 = com1 + np.array([radius * math.cos(fi[0]) * sinth[0], radius * math.sin(fi[0]) * sinth[0], radius * costh[0]])
	particle4 = particle3 + np.array([radius2 * math.cos(fi[1]) * sinth[1], radius2 * math.sin(fi[1]) * sinth[1], radius2 * costh[1]])
	return particle3, particle4
ar_ratio =  0.1 #Rmax over radius of sphere
chi = [10**(-3)] #Thermal energy over max potential energy of FENE spring
if ar_ratio == 0:
	hydro = "no"
else:
	hydro = "yes"

#Distribution constants
distr_constant = [1,5,10]
bins = 30
