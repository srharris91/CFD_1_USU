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
xu=get_data('output/xu.txt')
yv=get_data('output/yv.txt')
u=get_data('output/u.txt')
v=get_data('output/v.txt')
p=get_data('output/P.txt')
#error = get_data('error.txt')
#conv = get_data('output/convergence.txt')
# line = get_data('output/line.txt')
#iterations = get_data('Iteration.txt')
# print error,error.max(),error.min()

# X,Y = np.meshgrid(x,y)



fig=plt.figure(figsize=(4,15))
ax1=fig.add_subplot(411)
ax1.contourf(xu,y,u,20,cmap=plt.cm.jet,origin='lower')
ax1.set_xlabel(r'$x$')
ax1.set_ylabel(r'$y$')
ax1.set_title(r'$u$')
im = plt.imshow(u, interpolation='none', extent=(0,1,0,1), cmap=plt.cm.jet, norm=cm.colors.Normalize(vmax=(u).max(), vmin=(u).min()))
plt.colorbar(im)
# plt.show()
# plt.legend(loc='best',frameon=False,numpoints=1)
#fig.savefig('output/u.pdf',bbox_inches='tight')

#fig=plt.figure(figsize=(3,3))
ax1=fig.add_subplot(412)
ax1.contourf(x,yv,v,20,cmap=plt.cm.jet,origin='lower')
ax1.set_xlabel(r'$x$')
ax1.set_ylabel(r'$y$')
ax1.set_title(r'$v$')
im = plt.imshow(v, interpolation='none', extent=(0,1,0,1), cmap=plt.cm.jet, norm=cm.colors.Normalize(vmax=(v).max(), vmin=(v).min()))
plt.colorbar(im)
# plt.show()
# plt.legend(loc='best',frameon=False,numpoints=1)
#fig.savefig('output/uv.pdf',bbox_inches='tight')



ax1=fig.add_subplot(413)
#ax1.contourf(x,y,p,20,cmap=plt.cm.gray,origin='lower')
ax1.contourf(x,y,p,20,cmap=plt.cm.jet,origin='lower')
ax1.set_xlabel(r'$x$')
ax1.set_ylabel(r'$y$')
ax1.set_title(r'$P$')
im = plt.imshow(p, interpolation='none', extent=(0,1,0,1), cmap=plt.cm.jet, norm=cm.colors.Normalize(vmax=(p).max(), vmin=(p).min()))
plt.colorbar(im)

ax1=fig.add_subplot(414)
#ax1.contourf(x,y,p,20,cmap=plt.cm.gray,origin='lower')
ax1.quiver(x,y,u,v)
ax1.set_xlabel(r'$x$')
ax1.set_ylabel(r'$y$')
ax1.set_title(r'$quiver$')
#im = plt.imshow(p, interpolation='none', extent=(0,1,0,1), cmap=plt.cm.jet, norm=cm.colors.Normalize(vmax=(p).max(), vmin=(p).min()))
#plt.colorbar(im)
fig.savefig('output/uvp.pdf',bbox_inches='tight')
