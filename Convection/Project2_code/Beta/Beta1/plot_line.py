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


# get numerical data
line_10 = get_data('line_0_M10.txt')
line_20 = get_data('line_0_M20.txt')
line_40 = get_data('line_0_M40.txt')
# exact value
exact       = np.zeros((1000,2))
exact[:,0]  = np.linspace(0,np.sqrt(2),1000)
for i in range(0,len(exact[:,0])):
    if (exact[i,0]<=np.sqrt(2)/2.): # to get in between the element values
        exact[i,1]  = 100
    else:
        exact[i,1]  = 0

fig =plt.figure(figsize=(6,6))
ax1=fig.add_subplot(111)
ax1.plot(line_10[:,0],line_10[:,1],'bo-',label='Mesh = 10X10')
ax1.plot(line_20[:,0],line_20[:,1],'rs-',label='Mesh = 20X20')
ax1.plot(line_40[:,0],line_40[:,1],'g^-',label='Mesh = 40X40')
ax1.plot(exact   [:,0],exact   [:,1],'k-',label='exact')
ax1.set_xlabel('Distance along diagonal')
ax1.set_ylabel(r'$\phi$')
ax1.legend(loc='best',numpoints=1)
fig.savefig('lines.pdf',bbox_inches='tight')

