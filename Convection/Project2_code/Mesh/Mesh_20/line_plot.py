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
line_0  = get_data('line_0.txt')
line_0_9= get_data('line_0_9.txt')
line_1  = get_data('line_1.txt')
# exact value
exact       = np.zeros((1000,2))
exact[:,0]  = np.linspace(0,np.sqrt(2),1000)
for i in range(0,len(exact[:,0])):
    if (exact[i,0]<=np.sqrt(2)/2.): # to get in between the element values
        exact[i,1]  = 100
    else:
        exact[i,1]  = 0

fig =plt.figure(figsize=(4,4))
ax1=fig.add_subplot(111)
ax1.plot(line_0  [:,0],line_0  [:,1],'bo-',label=r'$\beta=0.$')
ax1.plot(line_0_9[:,0],line_0_9[:,1],'ro-',label=r'$\beta=0.9$')
ax1.plot(line_1  [:,0],line_1  [:,1],'go-',label=r'$\beta=1.$')
ax1.plot(exact   [:,0],exact   [:,1],'k-',label='exact')
ax1.set_xlabel('Distance along diagonal')
ax1.set_ylabel(r'$\phi$')
ax1.legend(loc='best',numpoints=1)
fig.savefig('lines.pdf',bbox_inches='tight')

