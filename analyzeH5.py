from __future__ import print_function

import numpy as np
import sys
import h5py
import matplotlib.pyplot as plt


c=['#000000','#E69F00','#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7']
c=c+c+c+c+c
m=[ 'o', 's', 'v', '<' ,'>', 'D', 'p', (8, 1, 0)]
m=m+m+m+m+m+m

ls=['-', ':', '-.']
ls=ls+ls+ls+ls+ls+ls+ls+ls



def analyze(sim_ID):

  print("opening ./{}.h5 ...".format(sim_ID))
  h5_file = h5py.File("./{}.h5".format(sim_ID), "r")

  time=h5_file['particles/atoms/charge/time']
  pos=h5_file['particles/atoms/position/value']
  # unwrap particle positions (they are wrapped in the box by default)
  pos=np.array(pos) + np.array(h5_file['particles/atoms/image/value'])*np.array(h5_file['particles/atoms/box/edges'])

  cations=np.arange(0,n_mono)
  cm_x=np.mean(pos[:,cations,0],axis=1)
  t=time
  
  return np.array(cm_x), np.array(t)



fig_1 = plt.figure()
ax_1 = fig_1.add_subplot(111)

n_mono=5
if True:
    for run_ID in [0]:

      sim_ID="electrophoresis-{}N-{}R".format(n_mono, run_ID)
      x,t=analyze(sim_ID)
      print (x,t)
      ax_1.plot(t, x, ls='-')
      

ax_1.legend(loc='best')
ax_1.set_ylabel(r'$x$ position: $x/\sigma$')
ax_1.set_xlabel(r'time: $ t $')


fig_1.savefig('migrate.pdf')
plt.close(fig_1)







