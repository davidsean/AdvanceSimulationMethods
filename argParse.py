from __future__ import print_function

import espressomd

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
parser.add_argument('-E_ext', type=float,  help='The E field in the x-direction', default=0.10)
parser.add_argument('-final_time', type=float,  help='Run up to this time (in MD units)', default=1000)


args = parser.parse_args()

E_ext=args.E_ext
run_ID=args.R
final_time=args.final_time
n_mono = args.n_mono

print ("Current simulation parameters:")
print ("E_ext: {}".format(E_ext))
print ("run_ID: {}".format(run_ID))
print ("final_time: {}".format(final_time))
print ("n_mono: {}".format(n_mono))

sim_ID="{}N-{}R".format(n_mono, run_ID)
print ("siulation output file: {}".format(sim_ID))


