from __future__ import print_function

import espressomd
from espressomd import shapes
from espressomd.io.writer import h5md

import numpy as np
import argparse



print( "Program Information: \n" )
print( espressomd.code_info.features() )
system = espressomd.System()
system.time_step = 0.01
system.cell_system.skin = 5.

box_l = 16
system.box_l = [box_l, box_l, box_l]
system.thermostat.set_langevin(kT=1.0, gamma=1.0)

E_ext=0.5
run_ID=0
n_mono=10

#===========================
#  File i/o
#===========================
sim_ID="constraints-{}R".format(run_ID)
h5md_file_name = "./{}.h5".format(sim_ID)
save_h5md = True




#===========================
#  counters
#===========================
type_counter=0
bond_counter=0
inter_counter=0


#===========================
# Particle types
#===========================
type_counter = 0
type_mono = type_counter
type_counter+=1
type_wall = type_counter

system.non_bonded_inter[type_mono, type_mono].lennard_jones.set_params( epsilon=1.0, sigma=1.0, cutoff=2.**(1.0 / 6.0), shift="auto")
system.non_bonded_inter[type_mono, type_wall].lennard_jones.set_params( epsilon=1.0, sigma=1.0, cutoff=2.**(1.0 / 6.0), shift="auto")
    

#===========================
# define the particle constraints (two walls)
#===========================

left_constr=system.constraints.add(particle_type=type_wall, penetrable=1, only_positive=0, shape=shapes.Wall(normal=[1,0,0], dist=1.12))
right_constr=system.constraints.add(particle_type=type_wall, penetrable=1, only_positive=0, shape=shapes.Wall(normal=[-1,0,0], dist=-(box_l-0.5)))


for p in range(0,n_mono):
    system.part.add(id=p, pos=(np.random.random(3) * system.box_l*0.5)+(box_l*0.25,0,0), type=type_mono )

cap=1
system.non_bonded_inter.set_force_cap(cap)
print( "Warming up... ",)
dist=system.analysis.mindist(p1=[type_mono], p2=[type_mono])
for t in range(4000):
    print( "Warming step: {} min_dist={} cap={}\r".format(t, dist, cap))
    system.integrator.run(200)
    dist=system.analysis.mindist(p1=[type_mono], p2=[type_mono])
    cap = cap + 1.
    system.non_bonded_inter.set_force_cap(cap)
    if (dist >= 1.1):
        break
print( "Done.")
system.non_bonded_inter.set_force_cap(10)
system.integrator.run(5000)
print( "Remove capping of LJ-interactions." )
system.non_bonded_inter.set_force_cap(0)
system.integrator.run(5000)
print( "Warmup finished.")

# apply a force
for p in system.part:
    p.ext_force = [E_ext, 0 , 0]

    
    
print("Equilibrate...",)
system.integrator.run(500)
print("** Done equilibration")


if save_h5md:
    h5_file = h5md.H5md(filename=h5md_file_name, write_pos=True, write_vel=False,
                write_force=False, write_species=True, write_mass=False, write_charge=True, write_ordered=True)
    h5_file.write()
    h5_file.flush()



if True:
    print("integrating...")
    while (system.time <1000):
        system.integrator.run(500)

        if save_h5md:
            h5_file.write()
            h5_file.flush()

