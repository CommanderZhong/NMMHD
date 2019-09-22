#!/usr/bin/env python
# -*- coding: UTF-8 -*-


#run with python2
import  numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.pylab as pylab

theta=np.arange(0,2*np.pi,0.02)
bnob=np.abs(np.cos(theta))

plt.figure(figsize=[10,3])
fig=1

titles=["(a)","(b)","(c)"]

for s in [0.5,1.0,2.0]:

    cfob=np.sqrt(1.0/2*(1+s+np.sqrt((1+s)**2-4*s*(np.cos(theta))**2)))
    csob=np.sqrt(1.0/2*(1+s-np.sqrt((1+s)**2-4*s*(np.cos(theta))**2)))

    plt.subplot(1,3,fig,projection='polar')
    fig+=1

    plt.ylim((0,2))

    
    plt.plot(theta,cfob)    #label='cf'
    plt.plot(theta,csob)    #label='cs'
    plt.plot(theta,bnob)    #label='bn'
    plt.plot(np.arange(0,2,0.02)*0+np.pi/12,np.arange(0,2,0.02),'--',color='black')
    
    
    plt.axis('off')
    nmark=theta.size/8
    plt.text(theta[nmark],cfob[nmark]+0.02,'f',fontsize=15)
    plt.text(theta[nmark],csob[nmark]-0.3,'s',fontsize=15)
    plt.text(theta[nmark],bnob[nmark]+0.02,'t',fontsize=15)
    plt.text(-np.pi/15,1,'1',fontsize=15)
    plt.text(-np.pi/30,2,'H',fontsize=15)
    plt.text(np.pi/36,1.7,r'$\theta$',fontsize=15)
    plt.text(-np.pi/1.90,2.05,titles[fig-2],fontsize=15)

    plt.annotate("",xy=(0,2),xytext=(np.pi,2),arrowprops=dict(arrowstyle="-|>"))
    plt.annotate("",xy=(np.pi/11.5,1.5),xytext=(np.pi/24,1.6),arrowprops=dict(arrowstyle="->"))
    plt.annotate("",xy=(-0.01,1.5),xytext=(np.pi/24,1.6),arrowprops=dict(arrowstyle="->"))

plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.subplots_adjust(top = 1, bottom = 0.05, right = 1, left = 0, hspace = 0, wspace = 0)

plt.savefig("figure1.eps")
#plt.show()
