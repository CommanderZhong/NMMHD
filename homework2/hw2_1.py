##homework2_problem1

import numpy as np 
import matplotlib.pyplot as plt 

def upwind(nx,c,t):

    dx=3.0/nx 
    dt=c*dx
    nt=int(round(t/dt))
    # t=nt*dt-dt
    
    x=np.arange( -1.0 , 2.0 , dx ) 

    u=np.zeros( nx )
    for i in np.arange(nx):
        if x[i] < -0.4 :
            u[i]= 0.0
        elif x[i] <-0.2 :
            u[i]= 1.0 - np.abs( x[i] + 0.3 )/0.1
        elif x[i] < -0.1 :
            u[i]= 0.0
        elif x[i] < -0.0 :
            u[i]= 1.0
        else :
            u[i]= 0.0
    
    #analytic solution
    nx_ture=5000
    dx_ture=3.0/nx_ture 
    
    x_ture=np.arange( -1.0 , 2.0 , dx_ture ) 

    u_ture=np.zeros( nx_ture )
    for i in np.arange(nx_ture):
        if x_ture[i] < -0.4 :
            u_ture[i]= 0.0
        elif x_ture[i] <-0.2 :
            u_ture[i]= 1.0 - np.abs( x_ture[i] + 0.3 )/0.1
        elif x_ture[i] < -0.1 :
            u_ture[i]= 0.0
        elif x_ture[i] < -0.0 :
            u_ture[i]= 1.0
        else :
            u_ture[i]= 0.0

    plt.plot(x_ture,u_ture,'b:')

    plt.plot(x_ture+t,u_ture,'b-')

    #numerical solution
    ut=np.zeros((nt,nx))

    for i in np.arange(nt):
        if i==0 :
            ut[i,:]=u
        else:
            
            for j in np.arange(nx):
                if j==0 :
                    ut[i,j]= ut[i-1,j]
                else:
                    ut[i,j]= ut[i-1,j] - dt/dx*(ut[i-1,j] - ut[i-1,j-1])

    plt.plot(x,ut[nt-1,:],'k--')
    plt.scatter(x,ut[nt-1,:],facecolors='none',edgecolors='k',s=50)

    plt.xlim(-1.0,2.0)
    plt.ylim(-0.2,1.2)
    plt.title('time='+str(t)+'  nx='+str(nx)+'  c='+str(c)+'    /Upwind')
    #plt.show()
    
    pass

def minmod(a,b):
    if np.abs(a)<=np.abs(b) and a*b>0:
        x=a
    elif np.abs(a)>np.abs(b) and a*b>0:
        x=b
    elif a*b<=0:
        x=0

    return x

def minmodd(nx,c,t):


    #analytic solution
    nx_ture=5000
    dx_ture=3.0/nx_ture 
    
    x_ture=np.arange( -1.0 , 2.0 , dx_ture ) 

    u_ture=np.zeros( nx_ture )
    for i in np.arange(nx_ture):
        if x_ture[i] < -0.4 :
            u_ture[i]= 0.0
        elif x_ture[i] <-0.2 :
            u_ture[i]= 1.0 - np.abs( x_ture[i] + 0.3 )/0.1
        elif x_ture[i] < -0.1 :
            u_ture[i]= 0.0
        elif x_ture[i] < -0.0 :
            u_ture[i]= 1.0
        else :
            u_ture[i]= 0.0

    plt.plot(x_ture,u_ture,'b:')

    plt.plot(x_ture+t,u_ture,'b-')


    c=0.95/max(u_ture)
    dx=3.0/nx 
    dt=c*dx
    nt=int(round(t/dt))
    # t=nt*dt-dt
    
    x=np.arange( -1.0 , 2.0 , dx ) 

    u=np.zeros( nx )
    for i in np.arange(nx):
        if x[i] < -0.4 :
            u[i]= 0.0
        elif x[i] <-0.2 :
            u[i]= 1.0 - np.abs( x[i] + 0.3 )/0.1
        elif x[i] < -0.1 :
            u[i]= 0.0
        elif x[i] < -0.0 :
            u[i]= 1.0
        else :
            u[i]= 0.0
    
    

    #numerical solution
    ut=np.zeros((nt,nx))

    for i in np.arange(nt):
        if i==0 :
            ut[i,:]=u
        else:
            
            for j in np.arange(nx):
                if j==0 :
                    ut[i,j]= ut[i-1,j]
                elif j==1 or j==nx-1:
                    ut[i,j]= ut[i-1,j] - dt/dx*(ut[i-1,j] - ut[i-1,j-1]) 
                else:
                    sigma1=minmod((ut[i-1,j] - ut[i-1,j-1])/dx,(ut[i-1,j+1] - ut[i-1,j])/dx)
                    sigma2=minmod((ut[i-1,j] - ut[i-1,j-1])/dx,(ut[i-1,j-1] - ut[i-1,j-2])/dx)
                    ut[i,j]= ut[i-1,j] - dt/dx*(ut[i-1,j] - ut[i-1,j-1]) + 0.5*dt/dx*(dx-dt)*(sigma1-sigma2)

    plt.plot(x,ut[nt-1,:],'k--')
    plt.scatter(x,ut[nt-1,:],facecolors='none',edgecolors='k',s=50)

    plt.xlim(-1.0,2.0)
    plt.ylim(-0.2,1.2)
    plt.title('time='+str(t)+'  nx='+str(nx)+'  c='+str(c)+'    /Minmod')
    #plt.show()
    
    pass

def laxW(nx,c,t):

    dx=3.0/nx 
    dt=c*dx
    nt=int(round(t/dt))
    # t=nt*dt-dt
    
    x=np.arange( -1.0 , 2.0 , dx ) 

    u=np.zeros( nx )
    for i in np.arange(nx):
        if x[i] < -0.4 :
            u[i]= 0.0
        elif x[i] <-0.2 :
            u[i]= 1.0 - np.abs( x[i] + 0.3 )/0.1
        elif x[i] < -0.1 :
            u[i]= 0.0
        elif x[i] < -0.0 :
            u[i]= 1.0
        else :
            u[i]= 0.0
    
    #analytic solution
    nx_ture=5000
    dx_ture=3.0/nx_ture 
    
    x_ture=np.arange( -1.0 , 2.0 , dx_ture ) 

    u_ture=np.zeros( nx_ture )
    for i in np.arange(nx_ture):
        if x_ture[i] < -0.4 :
            u_ture[i]= 0.0
        elif x_ture[i] <-0.2 :
            u_ture[i]= 1.0 - np.abs( x_ture[i] + 0.3 )/0.1
        elif x_ture[i] < -0.1 :
            u_ture[i]= 0.0
        elif x_ture[i] < -0.0 :
            u_ture[i]= 1.0
        else :
            u_ture[i]= 0.0

    plt.plot(x_ture,u_ture,'b:')

    plt.plot(x_ture+t,u_ture,'b-')

    #numerical solution
    ut=np.zeros((nt,nx))

    for i in np.arange(nt):
        if i==0 :
            ut[i,:]=u
        else:
            
            for j in np.arange(nx):
                if j==0 :
                    ut[i,j]= ut[i-1,j]
                elif j==nx-1:
                    ut[i,j]= ut[i-1,j] 
                else:
                    ut[i,j]= ut[i-1,j] - dt/dx*0.5*(ut[i-1,j+1] - ut[i-1,j-1]) + (dt/dx)**2*0.5*(ut[i-1,j+1] - 2*ut[i-1,j]+ ut[i-1,j-1])

    plt.plot(x,ut[nt-1,:],'k--')
    plt.scatter(x,ut[nt-1,:],facecolors='none',edgecolors='k',s=50)

    plt.xlim(-1.0,2.0)
    plt.ylim(-0.2,1.2)
    plt.title('time='+str(t)+'  nx='+str(nx)+'  c='+str(c)+'    /Lax-Wendroff')
    #plt.show()
    
    pass


def euler(nx,c,t,step):

    dx=3.0/nx 
    dt=c*dx
    nt=int(round(t/dt))
    # t=nt*dt-dt
    
    x=np.arange( -1.0 , 2.0 , dx ) 

    #origin data
    u=np.zeros( nx )
    for i in np.arange(nx):
        if x[i] < -0.4 :
            u[i]= 0.0
        elif x[i] <-0.2 :
            u[i]= 1.0 - np.abs( x[i] + 0.3 )/0.1
        elif x[i] < -0.1 :
            u[i]= 0.0
        elif x[i] < -0.0 :
            u[i]= 1.0
        else :
            u[i]= 0.0
    
    #analytic solution
    nx_ture=5000
    dx_ture=3.0/nx_ture 
    
    x_ture=np.arange( -1.0 , 2.0 , dx_ture ) 

    u_ture=np.zeros( nx_ture )
    for i in np.arange(nx_ture):
        if x_ture[i] < -0.4 :
            u_ture[i]= 0.0
        elif x_ture[i] <-0.2 :
            u_ture[i]= 1.0 - np.abs( x_ture[i] + 0.3 )/0.1
        elif x_ture[i] < -0.1 :
            u_ture[i]= 0.0
        elif x_ture[i] < -0.0 :
            u_ture[i]= 1.0
        else :
            u_ture[i]= 0.0

    plt.plot(x_ture,u_ture,'b:')

    plt.plot(x_ture+step*dt,u_ture,'b-')

    #numerical solution
    ut=np.zeros((nt,nx))
    ut[0,:]=u
    for i in np.arange(nt):
        if i==nt-1 :
            ut[i,:]=u

        else:
            
            for j in np.arange(nx):
                if j==nx-1 or j==0:
                    ut[i+1,j]= ut[i,j]
                else:
                    ut[i+1,j]= ut[i,j] - dt/dx*0.5*(ut[i,j+1] - ut[i,j-1])

    plt.plot(x,ut[step,:],'k--')
    plt.scatter(x,ut[step,:],facecolors='none',edgecolors='k',s=50)

    plt.xlim(-1.0,2.0)
    plt.ylim(-0.2,1.2)
    plt.title('step='+str(step)+'   time='+str(step*dt)+'  nx='+str(nx)+'  c='+str(c)+'    /Euler')
    
    pass


if __name__ == "__main__":

    

    nx=500
    c=0.95
    t=0.5

    size=16

    titles=['(a)','(b)','(c)','(d)','(e)','(f)']

    #different nx[200,500,1000,2500] at c=0.95,
    nxs=np.array([200,500,1000,2500])
    
    plt.figure(figsize=[size,size])
    for i in np.arange(4):
        
        plt.subplot(2,2,i+1)
        upwind(nxs[i],c,t)
        plt.text(0.45,-0.35,titles[i],fontsize=15)

    plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.95, left = 0.05)
    plt.savefig("hw2_1_nx.eps")
    

    # diffrent c[0.05,0.5,0.95,1.0] at nx=500,
    cs=np.array([0.05,0.5,0.95,1.0])
    
    plt.figure(figsize=[size,size])
    for i in np.arange(6):
        plt.subplot(3,2,i+1)
        if i==4:
            laxW(nx,c,t)
            plt.text(0.45,-0.35,titles[i],fontsize=15)
        elif i==5:
            minmodd(nx,c,t)
            plt.text(0.45,-0.35,titles[i],fontsize=15)
        else:   
            
            upwind(nx,cs[i],t)
            plt.text(0.45,-0.35,titles[i],fontsize=15)

    plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.95, left = 0.05)
    plt.savefig("hw2_1_c.eps")

    # diffrent time=[0.25,0.5,1.0,1.5]
    time=np.array([0.25,0.5,1.0,1.5])
    
    plt.figure(figsize=[size,size])
    for i in np.arange(4):
        
        plt.subplot(2,2,i+1)
        upwind(nx,c,time[i])
        plt.text(0.45,-0.35,titles[i],fontsize=15)

    plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.95, left = 0.05)
    plt.savefig("hw2_1_time.eps")
    
    #different scheme
    plt.figure(figsize=[size,size/2])
    for i in np.arange(2):
        
        plt.subplot(1,2,i+1)
        if i==0:
            laxW(nx,c,t)
            plt.text(0.45,-0.35,titles[i],fontsize=15)
        else:
            minmodd(nx,c,t)
            plt.text(0.45,-0.35,titles[i],fontsize=15)

    plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.05)
    plt.savefig("hw2_1_scheme.eps")

    # Euler
    steps=np.array([1,5,10,50])
    
    plt.figure(figsize=[size,size])
    for i in np.arange(4):
        
        plt.subplot(2,2,i+1)
        euler(nx,1,t,steps[i])
        plt.text(0.45,-0.35,titles[i],fontsize=15)

    plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.95, left = 0.05)
    plt.savefig("hw2_1_euler.eps")
    
    


    pass