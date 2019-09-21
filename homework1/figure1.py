#!/usr/bin/env python
# -*- coding: UTF-8 -*-


#run with python2
import  numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.pylab as pylab

theta=np.arange(0,2*np.pi,0.02)
bnob=np.abs(np.cos(theta))

plt.figure(figsize=[15,5])
fig=1
for s in [0.5,1.0,2.0]:

    cfob=np.sqrt(1.0/2*(1+s+np.sqrt((1+s)**2-4*s*(np.cos(theta))**2)))
    csob=np.sqrt(1.0/2*(1+s-np.sqrt((1+s)**2-4*s*(np.cos(theta))**2)))

    plt.subplot(1,3,fig,projection='polar')
    fig+=1

    plt.ylim((0,2))

    
    plt.plot(theta,cfob,label='cf')
    plt.plot(theta,csob,label='cs')
    plt.plot(theta,bnob,label='bn')
    plt.plot(np.arange(0,2,0.02)*0+np.pi/12,np.arange(0,2,0.02),'--',color='black')
    
    #plt.xticks([])
    #plt.yticks([])
    plt.axis('off')
    nmark=theta.size/8
    plt.text(theta[nmark],cfob[nmark]+0.02,'f')
    plt.text(theta[nmark],csob[nmark]-0.02,'s')
    plt.text(theta[nmark],bnob[nmark]+0.02,'t')
    plt.text(0,1,'1')
    plt.text(0,2,'H')
    plt.text(np.pi/36,1.6,r'$\theta$')
    plt.title('s='+str(s))

    plt.annotate("",xy=(0,2),xytext=(np.pi,2),arrowprops=dict(arrowstyle="-|>"))
    plt.annotate("",xy=(np.pi/12,1.5),xytext=(np.pi/18,1.6),arrowprops=dict(arrowstyle="->"))
    plt.annotate("",xy=(0,1.5),xytext=(np.pi/36,1.6),arrowprops=dict(arrowstyle="->"))

plt.savefig("figure1.eps")
plt.show()
