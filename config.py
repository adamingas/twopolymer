import numpy as np

noise = "yes" #yes or no
steps = 10000
runs = 250
constant = np.linspace(0,10,21).tolist()
time_step = 0.1
#Particles a b c d with initial positions:
ayinitial = -0.05
axinitial = 0.
azinitial = 0.
byinitial = 0.05
bxinitial = 0.
bzinitial = 0.
cyinitial = 0.05
cxinitial = 0.
czinitial = 100000.
dyinitial = -0.05
dxinitial = 0.
dzinitial = 100000.
initial_positions = [axinitial,ayinitial,azinitial,bxinitial,byinitial,bzinitial,cxinitial,cyinitial,czinitial,dxinitial,dyinitial,dzinitial]
ar_ratio =  0.1 #Rmax over radius of sphere
chi = [10**(-3)] #Thermal energy over max potential energy of FENE spring
if ar_ratio == 0:
	hydro = "no"
else:
	hydro = "yes" 
