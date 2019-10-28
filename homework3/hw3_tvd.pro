PRO hw3_tvd

;;-----------Initial value------------
gam=1.4   ;gamma
WL=[[0.445],[0.311],[8.928]]
WR=[[0.5],[0.0],[1.4275]]
num=2001L
t0=1.4
c=0.5  ;Courant coefficient
x=FINDGEN(num)/(num-1)*2-1    ;x from -1 to 1
dx=2./(num-1)
dt=dx*c
nt=t0/dt+1L

;;-----------Density-----------
;initial rho
rho0=FLTARR(num)
rho0[WHERE(x LT 0)]=WL[0,0]
rho0[WHERE(x GE 0)]=WR[0,0]

;Analytical solution rho
rho1=FLTARR(num)
rho1[WHERE(x LT -0.369)]=WL[0,0]
rho1[WHERE((x GE -0.369) AND (x LT -0.229))]=WL[0,0]-(WL[0,0]-0.345)/(0.369-0.229)*(x[WHERE((x GE -0.369) AND (x LT -0.229))]+0.369)  ;linear interpolation
rho1[WHERE((x GE -0.229) AND (x LT 0.214))]=0.345 ;from the Excel
rho1[WHERE((x GE 0.214) AND (x LT 0.347))]=1.304  ;from the Excel
rho1[WHERE(x GE 0.347)]=0.500                     ;from the Excel

;;-----------Mass Flux-----------
;initial m
m0=FLTARR(num)
m0[WHERE(x LT 0)]=WL[0,1]
m0[WHERE(x GE 0)]=WR[0,1]

;Analytical solution m
m1=FLTARR(num)
m1[WHERE(x LT -0.369)]=WL[0,1]
m1[WHERE((x GE -0.369) AND (x LT -0.229))]=WL[0,1]-(WL[0,1]-0.527)/(0.369-0.229)*(x[WHERE((x GE -0.369) AND (x LT -0.229))]+0.369)  ;linear interpolation
m1[WHERE((x GE -0.229) AND (x LT 0.214))]=0.527 ;from the Excel
m1[WHERE((x GE 0.214) AND (x LT 0.347))]=1.994  ;from the Excel
m1[WHERE(x GE 0.347)]=0.                        ;from the Excel

;;-----------Energy-----------
;initial E
E0=FLTARR(num)
E0[WHERE(x LT 0)]=WL[0,2]
E0[WHERE(x GE 0)]=WR[0,2]

;Analytical solution E
E1=FLTARR(num)
E1[WHERE(x LT -0.369)]=WL[0,2]
E1[WHERE((x GE -0.369) AND (x LT -0.229))]=WL[0,2]-(WL[0,2]-6.570)/(0.369-0.229)*(x[WHERE((x GE -0.369) AND (x LT -0.229))]+0.369)  ;linear interpolation
E1[WHERE((x GE -0.229) AND (x LT 0.214))]=6.570 ;from the Excel
E1[WHERE((x GE 0.214) AND (x LT 0.347))]=7.691  ;from the Excel
E1[WHERE(x GE 0.347)]=1.428                     ;from the Excel

;;Numerical solution with TVD scheme
;initial velocity
u=FLTARR(num)
u[WHERE(x LT 0)]=WL[0,1]/WL[0,0]
u[WHERE(x GE 0)]=WR[0,1]/WL[0,0]
;initial pressure
p=FLTARR(num)
p[WHERE(x LT 0)]=(gam-1)*(WL[0,2]-0.5*WL[0,0]*u[WHERE(x LT 0)]^2)
p[WHERE(x GE 0)]=(gam-1)*(WR[0,2]-0.5*WL[0,0]*u[WHERE(x GE 0)]^2)
;initial acoustic
a=FLTARR(num)
a[WHERE(x LT 0)]=SQRT(gam*p[WHERE(x LT 0)]/WL[0,0])
a[WHERE(x GE 0)]=SQRT(gam*p[WHERE(x GE 0)]/WR[0,0])
;density
rho=FLTARR(num,nt)
rho[*,0]=rho0
;mass
m=FLTARR(num,nt)
m[*,0]=m0
;Energy
E=FLTARR(num,nt)
E[*,0]=E0

FOR n=0,nt-2 DO BEGIN
  FOR j=1,num-2 DO BEGIN
    rho[j,n+1]=rho[j,n]-0.5*(tvd(rho,m,E,u,p,a,n,j,gam,c,num,0)-tvd(rho,m,E,u,p,a,n,j-1,gam,c,num,0))
    m[j,n+1]=m[j,n]-0.5*(tvd(rho,m,E,u,p,a,n,j,gam,c,num,1)-tvd(rho,m,E,u,p,a,n,j-1,gam,c,num,1))
    E[j,n+1]=m[j,n]-0.5*(tvd(rho,m,E,u,p,a,n,j,gam,c,num,2)-tvd(rho,m,E,u,p,a,n,j-1,gam,c,num,2))
  ENDFOR
  rho[num-1,n+1]=rho[num-1,n]
  m[num-1,n+1]=m[num-1,n]
  E[num-1,n+1]=E[num-1,n]
  u=m[*,n+1]/rho[*,n+1]
  p=(gam-1)*(E[*,n+1]-0.5*rho[*,n+1]*u^2)
  a=SQRT(gam*p/rho[*,n+1])
ENDFOR

;;--------Plot and Save Image-----------------
;rho
fig_tvd=PLOT(x,rho0,YTITLE='$\rho$',TITLE='Time=0.1400  /TVD',YRANGE=[0.2,1.4],':',XTICKFORMAT='(A6)')
fig_tvd.POSITION=[0.1,0.68,0.95,0.95]
fig_tvd=PLOT(x,rho1,/OVERPLOT)
fig_tvd=PLOT(x,rho[*,nt-1],/OVERPLOT,'--')
fig_tvd.SYMBOL='o'
fig_tvd.SYM_SIZE=1.2
;mass
fig_tvd=PLOT(x,m0,XTITLE='x',YTITLE='m',YRANGE=[-0.5,2.5],':',/CURR)
fig_tvd.POSITION=[0.1,0.08,0.95,0.34]
fig_tvd=PLOT(x,m1,/OVERPLOT)
fig_tvd=PLOT(x,m[*,nt-1],/OVERPLOT,'--')
fig_tvd.SYMBOL='o'
fig_tvd.SYM_SIZE=1.2
;energy
fig_tvd=PLOT(x,E0,YTITLE='E',YRANGE=[0,10],':',/CURR,XTICKFORMAT='(A6)')
fig_tvd.POSITION=[0.1,0.38,0.95,0.64]
fig_tvd=PLOT(x,E1,/OVERPLOT)
fig_tvd=PLOT(x,E[*,nt-1],/OVERPLOT,'--')
fig_tvd.SYMBOL='o'
fig_tvd.SYM_SIZE=1.2

fig_tvd.SAVE,'fig_tvd.pdf',RESOLUTION=512,/TRANSPARENT
fig_tvd.CLOSE
END