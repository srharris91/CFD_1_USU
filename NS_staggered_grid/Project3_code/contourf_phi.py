import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab,cm
from scipy.integrate import quad

def get_data(filename):
#    header=np.genfromtxt(filename,dtype=str)[0,:]
    data=np.genfromtxt(filename,skip_header=0)
    return data

def plot_data(xlabel,x,ylabel,y,filename):
    fig=plt.figure(figsize=(3,3))
    ax1=fig.add_subplot(111)
    ax1.plot(x,y,'o')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    #ax1.axis('equal')
    fig.savefig(filename,bbox_inches='tight')


x=get_data('output/x.txt')
y=get_data('output/y.txt')
u=get_data('output/u.txt')
#error = get_data('error.txt')
#conv = get_data('output/convergence.txt')
# line = get_data('output/line.txt')
#iterations = get_data('Iteration.txt')
# print error,error.max(),error.min()

# X,Y = np.meshgrid(x,y)



fig=plt.figure(figsize=(3,3))
ax1=fig.add_subplot(111)
ax1.contourf(x,y,u,20,cmap=plt.cm.gray,origin='lower')
ax1.set_xlabel(r'$x$')
ax1.set_ylabel(r'$y$')
im = plt.imshow(u, interpolation='none', extent=(0,1,0,1), cmap=plt.cm.gray, norm=cm.colors.Normalize(vmax=abs(u).max(), vmin=abs(u).min()))
plt.colorbar(im)
# plt.show()
# plt.legend(loc='best',frameon=False,numpoints=1)
fig.savefig('output/u.pdf',bbox_inches='tight')

#fig=plt.figure()
#ax1=fig.add_subplot(111)
#ax1.contourf(x,y,error,20,cmap=plt.cm.gray,origin='lower')
#ax1.set_xlabel(r'$x$',fontsize=24)
#ax1.set_ylabel(r'$y$',fontsize=24)
#im = plt.imshow(error, interpolation=None, extent=(0,1,0,1), cmap=plt.cm.gray, norm=cm.colors.Normalize(vmax=abs(error).max(), vmin=abs(error).min()))
#plt.colorbar(im)
## plt.show()
## plt.legend(loc='best',frameon=False,numpoints=1)
#fig.savefig('error.pdf',bbox_inches='tight')

#fig=plt.figure(figsize=(3,3))
#ax1=fig.add_subplot(111)
#ax1.loglog(conv[:,0],conv[:,1],linewidth=4)
#ax1.set_xlabel('Iteration')
#ax1.set_ylabel('RSS error')
#fig.savefig('output/convergence.pdf',bbox_inches='tight')

# fig =plt.figure(figsize=(3,3))
# ax1=fig.add_subplot(111)
# ax1.plot(line[:,0],line[:,1],'ko-')
# ax1.set_xlabel('Distance along diagonal')
# ax1.set_ylabel(r'$\u$')
# fig.savefig('output/line.pdf',bbox_inches='tight')
# fig=plt.figure()
# ax1=fig.add_subplot(111)
# ax1.plot(iterations[:,0],iterations[:,1],linewidth=4)
# ax1.set_xlabel(r'$\Omega$')
# ax1.set_ylabel(r'Iterations to convergence')
# fig.savefig('Iterations.pdf',bbox_inches='tight')

# fig=plt.figure()
# ax1=fig.add_subplot(111)
# ax1.plot(CL2,CDi2,'o-',linewidth=4,ms=8,label='AR = 2')
# ax1.plot(CL6,CDi6,'--o',linewidth=4,ms=8,label='AR = 6')
# ax1.plot(CL12,CDi12,'-*',linewidth=4,ms=8,label='AR = 12')
# ax1.set_xlabel(r'$C_L$',fontsize=24)
# ax1.set_ylabel(r'$C_{Di}$',fontsize=24)
# plt.legend(loc='best',frameon=False,numpoints=1)
# fig.savefig('Drag_Polar.pdf',bbox_inches='tight')

