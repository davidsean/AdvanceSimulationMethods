import numpy as np
import matplotlib.pyplot as plt

eta=1.
f=0.001

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

data=np.loadtxt("center_velocity.dat")

ax1.plot(data, marker='o', label=r'Sims')

d_y=28.
u_max=(1./(8.*eta))*f*(d_y)**2
ax1.plot((data*0+1)*u_max, marker='.', label=r'Theory: ${u_{max}}=%.3f$' %(u_max))

ax1.legend(loc='best')
ax1.set_xlabel(r'time')
ax1.set_ylabel(r'fluid velocity $v_x$')
#ax1.set_xscale('log')
#ax1.set_xlim(1, 30)
#ax1.set_ylim(0.0,1.0)

fig1.savefig('center.png')


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

data=np.loadtxt("profile_velocity.dat")
ax1.plot(data[:,0], data[:,1], marker='o', label=r'Sims')
y=np.arange(0,32,1)
v=-(1./(2.*eta)*f*(y-1.5)*(y-29.5))
ax1.plot(y, v, marker='o', label=r'Theory')

ax1.legend(loc='best')
ax1.set_xlabel(r'$y$-position: ($y$)')
ax1.set_ylabel(r'$x$-fluid velocity: $v_x$')
#ax1.set_xscale('log')
#ax1.set_xlim(1, 30)
#ax1.set_ylim(0.0,1.0)

fig1.savefig('profile.png')



exit()


