from __future__ import print_function

import espressomd
from espressomd import thermostat
from espressomd import interactions
from espressomd.io.writer import h5md
from espressomd import electrostatics
from espressomd import lb

import numpy as np
import argparse
import os
#import h5py




# System parameters
#############################################################
print(espressomd.code_info.features())

system = espressomd.System()


parser = argparse.ArgumentParser(description='Read simulation parameters')
parser.add_argument('-R', type=int,  help='The run ID', default=0)
parser.add_argument('-n_mono', type=int,  help='The number of monomers', default=5)
parser.add_argument('-E_ext', type=float,  help='The E field in the x-direction', default=0.10)
parser.add_argument('-final_time', type=float,  help='Run up to this time (in MD units)', default=1000)
parser.add_argument('-l_B', type=float,  help='The Bjerrum length', default=1.0)


args = parser.parse_args()

E_ext=args.E_ext
run_ID=args.R

final_time=args.final_time
l_bjerrum = args.l_B
n_mono = args.n_mono


box_l=100




#===========================
#  File i/o
#===========================
sim_ID="electrophoresis-{}N-{}R".format(n_mono, run_ID)
h5md_file_name = "./{}.h5".format(sim_ID)
save_h5md = True

#===========================





#===========================
#  Thermostat
#===========================
## Magic Lattice-Boltzmann parameters
LB_agrid  = 1
LB_viscosity = 1.0
LB_kT = 1.0
LB_time_step = 0.01
LB_density = 1.0
LB_gamma = 20


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

#cell_system.set_domain_decomposition(use_verlet_lists=False)
#print(system.cell_system.tune_skin())



### Particle types
type_counter = 0
type_mono = type_counter
type_counter+=1
type_cation = type_counter
type_counter+=1



for i in range(type_counter):
  for j in range(i+1):
    system.non_bonded_inter[i, j].lennard_jones.set_params(
      epsilon=1., sigma=1.,
      cutoff=2.**(1. / 6.), shift="auto")

for p in np.arange(0,n_mono):
    system.part.add(id=p, pos=np.random.random(3)*system.box_l, type=type_mono, q=1)

# monomer charge
for p in np.arange(0,n_mono):
    system.part[p].q=-1
print("\t** Placed {} monomers\n".format(n_mono))


cation_start=n_mono
cation_end=cation_start+n_mono
# counter-ions (from placed monomers)
for p in np.arange(cation_start,cation_end):
    system.part.add(id=p, pos=np.random.random(3)*system.box_l, type=type_cation, q=1)
print("\t** Placed {} countertions\n".format(cation_end-cation_start))



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



# apply electric force, and perform a sanity check
tot_q=0
for p in system.part:
    tot_q+=p.q
    p.ext_force = [p.q*E_ext, 0, 0] 
print ("total charge: ",tot_q)
if tot_q!=0: exit()




if True:
  #############################################################
  print("** Activating P3M...",)
  #############################################################
  p3m = electrostatics.P3M(bjerrum_length=l_bjerrum, accuracy=1e-4)
  system.actors.add(p3m)
  print("** Done P3M")

    

if True:
  #############################################################
  print("** Activating LB...",)
  ############################################################# 
  system.thermostat.turn_off()
  system.thermostat.set_lb(kT=1.0)
  # remove system momentum (usually from LD)
  system.galilei.galilei_transform()

  lb_fluid = lb.LBFluid_GPU(agrid=LB_agrid, visc=LB_viscosity, tau=LB_time_step, dens=LB_density, fric=LB_gamma)

  system.actors.add(lb_fluid)

  print("** Done LB")

print("Equilibrate...",)

system.integrator.run(5000)

print("** Done equilibration")


system.time=0
if save_h5md:
    h5_file = h5md.H5md(filename=h5md_file_name, write_pos=True, write_vel=False,
                write_force=False, write_species=True, write_mass=False, write_charge=True, write_ordered=True)
    h5_file.write()
    h5_file.flush()

print("integrating...")
while (system.time <final_time):
    system.integrator.run(500)

    if save_h5md:
        h5_file.write()
        h5_file.flush()
if save_h5md:  h5_file.close()

print("**Simulation finished")


