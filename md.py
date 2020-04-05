"""====================================================================
PREAMBLE
===================================================================="""

# Import useful Libraries
import numpy as np
from random import uniform
from sys import argv

# Define interaction constants for later use
epsilon = 4.91151e-21
sigma = 2.725

#Define mass to determine Average Kinetic Energy and Temprature
mass = 2.9915e-20 # kg

# Creates the size of one side of the box in angstroms
box = float(argv[1]) # angstroms

"""====================================================================
Class
===================================================================="""

# Define a class that will initalize each particle from an input file
class particle:

    """a simple particle"""

    def __init__(self, __array__):
        velocityx = uniform(-1.0,1.0)
        velocityy = uniform(-1.0,1.0)
        velocityz = uniform(-1.0,1.0)

        # Gives every particle a Six dimensional array to work with
        # through out the program
        self.P_vec = np.array([__array__[0], __array__[1],
            __array__[2], velocityx, velocityy, velocityz], float)

"""====================================================================
Initialization
===================================================================="""

# Load in input.xy to initalize the
coordinates = np.loadtxt("input.xy")

# Finds the number of particles in the system
size = coordinates.shape
size = size[0]

# Initialize a particle array that uses every single particle in the
# box as an index. Each index is assigned a six-dimensional vector
# for calculation.
particle_array = []
for p in range(size):
    foo = particle(np.array(coordinates[p][:]))
    particle_array.append(foo.P_vec)

"""====================================================================
Definitions
===================================================================="""

# Define interaction between two forces
def Force(P1,P2):

    def repulse(P1,P2):
        x_dir = P1[0] - P2[0]
        y_dir = P1[1] - P2[1]
        z_dir = P1[2] - P2[2]
        R = np.sqrt((x_dir)**2 + (y_dir)**2 + z_dir**2)
        product_1 = 24*epsilon/R**2
        diff_1 = 2*(sigma/R)**12
        diff_2 = (sigma/R)**6
        product_2 = diff_1 - diff_2
        force = product_1*product_2
        return np.array([x_dir*force, y_dir*force,
            z_dir*force,P1[3], P1[4],P1[5]], float)

    return (repulse(P1,P2))

# Define the Verlet Algorithm
def vverlet(r1,r2,h):
    f_half_step = 0.5 * h * Force(r1,r2)
    vx_half_step = r1[3] + f_half_step[0]
    vy_half_step = r1[4] + f_half_step[1]
    vz_half_step = r1[5] + f_half_step[2]
    r1[0] += h*vx_half_step
    r1[1] += h*vy_half_step
    r1[2] += h*vz_half_step
    Return = np.array([r1[0],r1[1],r1[2],
        vx_half_step,vy_half_step,vz_half_step],float)
    k = h*Force(Return,r2)
    Return[3] += k[0]
    Return[4] += k[1]
    Return[5] += k[2]
    return Return

# Define the interaction between the box and a particle. This makes
# the system N,V,T
def boundary(P):
    if P[0] >= box:
        P[0] = box
        P[3] = -P[3]
    elif P[0] <= -box:
        P[0] = -box
        P[3] = -P[3]
    elif P[1] >= box:
        P[1] = box
        P[4] = -P[4]
    elif P[1] <= -box:
        P[1] = -box
        P[4] = -P[4]
    elif P[2] >= box:
        P[2] = box
        P[5] = -P[5]
    elif P[2] <= -box:
        P[2] = -box
        P[5] = -P[5]
    return P

# Define kinetic energy
def T(P1):
    vx2 = P1[3]**2
    vy2 = P1[4]**2
    vz2 = P1[5]**2
    return 0.5*mass*(vx2 + vy2 + vz2)

# Define Potential
def U(P1,P2):
    x_dir = P1[0] - P2[0]
    y_dir = P1[1] - P2[1]
    z_dir = P1[2] - P2[2]
    R = np.sqrt((x_dir)**2 + (y_dir)**2 + z_dir**2)
    product_1 = 4*epsilon
    diff_1 = (sigma/R)**12
    diff_2 = (sigma/R)**6
    product_2 = diff_1 - diff_2
    potential = product_1*product_2
    return potential

"""====================================================================
Main
===================================================================="""

# Import visualization libraries
import matplotlib.pyplot as mp
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

# Create a timerange of one picosecond
t0, tf = 0, 1
h= (tf-t0)/1000
tp = np.arange(t0,tf,h)

# Create the figure
fig = mp.figure()
ax = fig.add_subplot(111, projection='3d',
        autoscale_on=False, xlim=(-box, box),
        ylim=(-box, box), zlim=(-box,box))

# Create the file to view the energy in the system and the temprature.
__file__ = open("output.dat", "w")

# Initialize the plot
def init():
    for k in range(size):
        coord = particle_array[k]
        ax.scatter(coord[0], coord[1], coord[2], c='b' , marker='o')
        energy = 0.0
        kinetic_sum = 0.0

# Calculate and plot the position of every particle in the system. Then
# Calculate the temprature of the system.
def animate(i):
    energy, kinetic_sum, Temprature = 0.0, 0.0, 0.0
    mp.cla()
    ax.set_xlim(-box,box)
    ax.set_ylim(-box,box)
    ax.set_zlim(-box,box)

    # Calculation loop
    for k in range(size):
        for l in range(size):
            if k != l:

                # Calculate the position and interaction for every
                # particle in the system.
                particle_array[k] = boundary(vverlet(particle_array[k],
                        particle_array[l],h))
                energy += U(particle_array[k], particle_array[l])
        energy += T(particle_array[k])
        kinetic_sum += T(particle_array[k])

        # Plot the Particle's position
        coord = particle_array[k]
        ax.scatter(coord[0],coord[1],coord[2], c='b', marker='o')

    # Calculate the temprature then print the energy and Temprature
    kinetic_sum = kinetic_sum / size
    Temprature = (2/3)*kinetic_sum/1.38e-23
    print(energy, Temprature, file=__file__)

# Load the animation.
anix = animation.FuncAnimation(fig, animate, np.arange(t0,tf,h),
                              interval=h, blit=False, init_func=init)

# Save the animation
anix.save('myAnimation.mp4', writer='ffmpeg', fps=60)



