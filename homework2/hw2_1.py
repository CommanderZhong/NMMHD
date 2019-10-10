##homework2_problem1

import numpy as np 
import matplotlib.pyplot as plt 

def upwind(u0,nx,nt):

    # dx=dx
    # dt=c*dx
    # t=time
    # nt=int(round(t/dt))+1
    # print nt,round(t/dt),int(t/dt),t/dt
    
    u=np.zeros((nt,nx))

    for i in np.arange(nt):
        if i==0 :
            u[i,:]=u0
        else:
            
            for j in np.arange(nx):
                if j==0 :
                    pass
                else:
                    u[i,j]= u[i-1,j] - dt/dx*(u[i-1,j] - u[i-1,j-1])

    return u[nt-1,:]

if __name__ == "__main__":

    t=0.5
    c=1
    nt=100
    dt=t/nt

    dx=dt/c
    nx=int(round(3.0/dx))
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

    plt.plot(x,u,'b:')
    
    
    ut=upwind(u,nx,nt)

    plt.plot(x+t,u,'b-')
    plt.plot(x,ut,'k--')
    plt.scatter(x,ut,facecolors='none',edgecolors='k',s=50)
    plt.show()
    
    pass