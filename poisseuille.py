from __future__ import print_function

import espressomd
from espressomd import shapes
from espressomd import lb
from espressomd import lbboundaries


import numpy as np
import argparse

print( " " )
print( "===================================================" )
print( "=            Lattice-Boltzmann Fluid              =" )
print( "===================================================" )
print( " " )

print( "Program Information: \n" )
print( espressomd.code_info.features() )
system = espressomd.System()
system.time_step = 0.01
system.cell_system.skin = 0.2

box_l = 32
system.box_l = [box_l, box_l, box_l]

# the force density on the fluid nodes
force_density=0.001

lbf = lb.LBFluid_GPU(agrid=1, dens=1, visc=1, tau=0.01, ext_force=[force_density, 0, 0])

system.actors.add(lbf)


system.thermostat.set_lb(kT=0)

# create the boundary "shape"
upper_wall=shapes.Wall(normal=[0,1,0], dist=1.5)
lower_wall=shapes.Wall(normal=[0,-1,0], dist=-(box_l-1.5))

# from these shapes, define the LB boundary
upper_bound=lbboundaries.LBBoundary(shape=upper_wall)
lower_bound=lbboundaries.LBBoundary(shape=lower_wall)

system.lbboundaries.add(upper_bound)
system.lbboundaries.add(lower_bound)


# save the center x-velocity in this file
center_output = open("center_velocity.dat", "w")

max_time=1000
for i in range(max_time):
    system.integrator.run(500)
    print("Running simulation, step: {}/{}\r".format(i,max_time))

    vx=(lbf[box_l/2, box_l/2, box_l/2].velocity[0])
    center_output.write("{}\n".format(vx))

# save the x-velocity  of a "line" of fluid
profile_output = open("profile_velocity.dat", "w")
for i in range(box_l):
    vx=(lbf[0, i, 0].velocity[0])
    profile_output.write("{} \t {}\n".format(i,vx))


center_output.close()
profile_output.close()


