from __future__ import print_function

import espressomd
from espressomd import thermostat
from espressomd import interactions
from espressomd.io.writer import h5md
from espressomd import polymer

import numpy as np
import argparse
import os






# System parameters
#############################################################
print(espressomd.code_info.features())

system = espressomd.System()


parser = argparse.ArgumentParser(description='Read simulation parameters')
parser.add_argument('-R', type=int,  help='The run ID', default=0)
parser.add_argument('-n_mono', type=int,  help='The number of monomers', default=5)
parser.add_argument('-final_time', type=float,  help='Run up to this time (in MD units)', default=1000)


args = parser.parse_args()

run_ID=args.R
final_time=args.final_time
n_mono = args.n_mono

box_l=100

#===========================
#  File i/o
#===========================
sim_ID="poly-{}N-{}R".format(n_mono, run_ID)
h5md_file_name = "./{}.h5".format(sim_ID)
save_h5md = True

#===========================



#if no seed is provided espresso generates a seed

epsilon=1.0
sigma=1.0

## counters
type_counter=0
bond_counter=0
inter_counter=0


# set up the integrator and thermostat
system.time_step = 0.01
skin=10.0
system.cell_system.skin = skin
system.box_l = [box_l, box_l, box_l]
system.thermostat.set_langevin(kT=1.0, gamma=1.0)


### Particle types
type_counter = 0
type_mono = type_counter
type_counter+=1

system.non_bonded_inter[type_mono, type_mono].lennard_jones.set_params( epsilon=1., sigma=1., cutoff=2.**(1. / 6.), shift="auto")

# See PRL-A, Volume 33, number 5, May 1986, Gary S. Grest and Kurt Kremer
k_fene = (30.*epsilon)/(sigma*sigma)
r_fene = 1.5*sigma 
fene = interactions.FeneBond(k=k_fene, d_r_max=r_fene)
system.bonded_inter.add(fene)


poly_start=0
polymer.create_polymer(N_P = 1, bond_length = 0.97, MPC=n_mono, start_id=poly_start, bond=fene, type_poly_neutral=type_mono, type_poly_charged=type_mono, mode=1)
poly_end=poly_start+n_mono


# this is to init a polymer a s straight line
#for p in range(0,n_mono):
#    system.part.add(id=p, pos=[p+0.5, 1,1], type=type_mono )
#for p in range(1,n_mono):
#    system.part[p].add_bond((fene,p-1))
print("\t** Placed {} monomers\n".format(n_mono))




cap=1
system.non_bonded_inter.set_force_cap(cap)
print( "Warming up... ",)
dist=system.analysis.mindist()

for t in range(4000):
    print( "Warming step: {} min_dist={} cap={}\r".format(t, dist, cap))
    system.integrator.run(200)
    dist=system.analysis.mindist()
    cap = cap + 1.
    system.non_bonded_inter.set_force_cap(cap)
    if (dist >= 0.9):
        break
print( "Done.")
print( "Remove capping of LJ-interactions." )
system.non_bonded_inter.set_force_cap(0)
print( "Warmup finished.")



    

if save_h5md:
    h5_file = h5md.H5md(filename=h5md_file_name, write_pos=True, write_vel=False,
                write_force=False, write_species=True, write_mass=False, write_charge=True, write_ordered=True)
    h5_file.write()
    h5_file.flush()

system.time=0

print("integrating...")
while (system.time <final_time):
    system.integrator.run(5000)

    if save_h5md:
        h5_file.write()
        h5_file.flush()

if save_h5md:  h5_file.close()
print("**Simulation finished")



