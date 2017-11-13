# driftsim.py
#
# Omkar H. Ramachandran
# omkar.ramachandran@colorado.edu
#
# Radial drift of dust particles in a protoplanetary disk
#

import numpy as np
import matplotlib.pyplot as plt

# Constants. Everything's in SI
AU = 1.46e11
G = 6.67e-11
rho_s = 2e-8
C_D = 0.48

# Number of timesteps, dt
Nt = 5000
dt = 1e11

# Ratio of gas to dust velocity
beta = 0.996

# Number of dust particles
Npart = 101

# Radius of dust particles. Constant for now. Make distribution later
Rpart = 1e-3
A = np.pi*Rpart**2

# Drag parameter
alpha = 0.5*C_D*A*rho_s*(1-beta**2)
#alpha 1e-12

# Mass of central star
M = 2e30

# Radius of disk. For now, stay outside frost line.
Rmax = 40*AU
Rmin = 10*AU

# Density distribution parameters
rho_mean = 3e3
rho_var = 1e3
N_rho = 101
rho_max = 5e3
rho = np.linspace(0,rho_max,N_rho)
drho = rho_max/N_rho

# Generate distribution. Divide by sum to get the approximation of the continuous
# distribution to sum to 1
P_rho =  np.exp(-(rho-rho_mean)**2/2/rho_var**2)/np.sqrt(2*np.pi*rho_var)*drho
P_rho /= np.sum(P_rho)

# Create list to hold Particle data. 
Part = np.zeros([Nt,Npart,5])

# Populate list. Start everything off at v_theta = v_c and v_r = 0.
# This is a one time computation, so don't worry about optimization
print("Generating particles:")
for i in range(Npart):
    # Generate R between Rmin and Rmax
    Part[0,i,0] = Rmax
    # Generate theta between 0 and 2*pi
    Part[0,i,1] = np.random.uniform()*2*np.pi
    # Enforce velocity initial conditions
    Part[0,i,2] = 0.
    Part[0,i,3] = np.sqrt(G*M/Part[0,i,0]**3)
    Part[0,i,4] = 4/3.*np.pi*Rpart**3*np.random.choice(rho,p=P_rho)
print(Part[0,:,4])

# Enuma Elish!
print("Look up and behold, Enuma Elish!")
for i in range(Nt-1):
    print("Running timestep",i)
    for j in range(Npart):
        # r^{t} and v_th^{t} are needed. Store them seperately
        r = Part[i,j,0]
        if(r < 1.*AU):
            Part[i+1,j,0] = r
            Part[i+1,j,1] = Part[i,j,1]
            Part[i+1,j,2] = 0
            Part[i+1,j,3] = 0
            Part[i+1,j,4] = Part[i,j,4]
        v_th = Part[i,j,3]

        # v_th equation
        Part[i+1,j,3] = Part[i,j,3] - alpha/Part[i,j,4]*Part[i,j,3]**2*dt
#        Part[i+1,j,3] = Part[i,j,3] - alpha/Part[i,j,4]/Part[i,j,0]*Part[i,j,3]*dt
        # th equation
        Part[i+1,j,1] = Part[i,j,1] + Part[i+1,j,3]*dt
        if(Part[i+1,j,1] > 2*np.pi):
            while(Part[i+1,j,1] > 2*np.pi):
                Part[i+1,j,1] -= 2*np.pi

        # r equation
        Part[i+1,j,0] = Part[i,j,0] + Part[i,j,2]
        Part[i+1,j,2] = Part[i,j,2] + (-G*M/r**2 + r*v_th**2)*dt
        # Mass
        Part[i+1,j,4] = Part[i,j,4]
