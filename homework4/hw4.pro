PRO hw4
;----------intial value----------
nx=401
gm=5./3 ;gamma
mu=1 ;mu
Hx=5./SQRT(4*!PI) ;Hx or Bx
x0=0.
bt=2 ;beta
uv=22
CFL=0.95
c=CFL/uv ;courant coeffitient

dx=4./(nx-1) ;-3---1
dt=dx*c
nt=FLOOR(0.1/dt)+1
x=INDGEN(nx)/(nx-1.)*4-3
loc=WHERE(x LE 0,count)

;------------------fast shock-------------
Wf=MAKE_ARRAY(nx,nt,7,/DOUBLE) 
Wf[0:count-1,0,*]=TRANSPOSE(REBIN([2.121,4.981,-13.27,-0.163,-0.6521, 2.572/SQRT(4*!PI), 10.29/SQRT(4*!PI)],7,count))
Wf[count:nx-1,0,*]=TRANSPOSE(REBIN([1, 1,-15.3, 0, 0, 1./SQRT(4*!PI), 4./SQRT(4*!PI)],7,nx-count))
U=MAKE_ARRAY(nx,nt,7,/DOUBLE)
U[*,0,0]=Wf[*,0,0] &  U[*,0,5:6]=Wf[*,0,5:6]
FOR j=0,nx-1 DO BEGIN
  U[j,0,2:4]=Wf[j,0,2:4]*Wf[j,0,0]
ENDFOR
U[*,0,1]=Wf[*,0,0]*TOTAL(REFORM(Wf[*,0,2:4]^2),2)+TOTAL(REFORM(Wf[*,0,5:6]^2),2)+bt*Wf[*,0,1]/(gm-1)
FOR n=0,nt-2 DO BEGIN
  U[0,n+1,*]=U[0,n,*]
  U[nx-1,n+1,*]=U[nx-1,n,*]
  Wf[0,n+1,*]=U2W(REFORM(U[0,n+1,*]),gm) & Wf[nx-1,n+1,*]=U2W(REFORM(U[nx-1,n+1,*]),gm)
  FOR j=1,nx-2 DO BEGIN
;    U[j,n+1,*]=(U[j+1,n,*]+U[j-1,n,*])*0.5-c*0.5*(U2F(REFORM(U[j+1,n,*]),gm)-U2F(REFORM(U[j-1,n,*]),gm)) ;Lax
    Um=0.5*(U[j,n,*]+U[j+1,n,*])-0.5*c*(U2F(REFORM(U[j+1,n,*]),gm)-U2F(REFORM(U[j,n,*]),gm))
    Un=0.5*(U[j,n,*]+U[j-1,n,*])-0.5*c*(U2F(REFORM(U[j,n,*]),gm)-U2F(REFORM(U[j-1,n,*]),gm))
    U[j,n+1,*]=U[j,n,*]-c*(U2F(REFORM(Um),gm)-U2F(REFORM(Un),gm)) ;Lax-Wendroff
    Wf[j,n+1,*]=U2W(REFORM(U[j,n+1,*]),gm)
  ENDFOR
ENDFOR

;------------------plot image------------
fig_hw4=PLOT(x,Wf[*,0,0],XTICKFORMAT='(A6)',':',YTITLE='$\rho$');
fig_hw4.POSITION=[0.1,0.765,0.47,0.95]
fig_hw4.YRANGE=[0.4,2.4]
fig_hw4=PLOT(x,Wf[*,nt-1,0],/OVERPLOT,'--')
fig_hw4.SYMBOL='o'
fig_hw4=PLOT(x,Wf[*,0,1],/CURR,XTICKFORMAT='(A6)',':',YTITLE='$p$');
fig_hw4.POSITION=[0.1,0.545,0.47,0.735]
fig_hw4.YRANGE=[-1,6]
fig_hw4=PLOT(x,Wf[*,nt-1,1],/OVERPLOT,'--')
fig_hw4.SYMBOL='o'
fig_hw4=PLOT(x,Wf[*,0,2],/CURR,XTICKFORMAT='(A6)',':',YTITLE='$v_x$');
fig_hw4.POSITION=[0.1,0.325,0.47,0.515]
fig_hw4.YRANGE=[-19,-12]
fig_hw4=PLOT(x,Wf[*,nt-1,2],/OVERPLOT,'--')
fig_hw4.SYMBOL='o'
fig_hw4=PLOT(x,Wf[*,0,3],/CURR,':',YTITLE='$v_y$',XTITLE='$x$');
fig_hw4.POSITION=[0.1,0.105,0.47,0.295]
fig_hw4.YRANGE=[-0.25,0.3]
fig_hw4=PLOT(x,Wf[*,nt-1,3],/OVERPLOT,'--')
fig_hw4.SYMBOL='o'
fig_hw4=PLOT(x,Wf[*,0,4],/CURR,XTICKFORMAT='(A6)',':',YTITLE='$v_z$');
fig_hw4.POSITION=[0.58,0.765,0.99,0.95]
fig_hw4.YRANGE=[-1,1]
fig_hw4=PLOT(x,Wf[*,nt-1,4],/OVERPLOT,'--')
fig_hw4.SYMBOL='o'
fig_hw4=PLOT(x,Wf[*,0,5],/CURR,XTICKFORMAT='(A6)',':',YTITLE='$H_y$');
fig_hw4.POSITION=[0.58,0.545,0.99,0.735]
fig_hw4.YRANGE=[0,0.8]
fig_hw4=PLOT(x,Wf[*,nt-1,5],/OVERPLOT,'--')
fig_hw4.SYMBOL='o'
fig_hw4=PLOT(x,Wf[*,0,6],/CURR,XTICKFORMAT='(A6)',':',YTITLE='$H_z$');
fig_hw4.POSITION=[0.58,0.325,0.99,0.515]
fig_hw4.YRANGE=[0.,3.5]
fig_hw4=PLOT(x,Wf[*,nt-1,6],/OVERPLOT,'--')
fig_hw4.SYMBOL='o'
fig_hw4=PLOT(x,U[*,0,1],/CURR,':',YTITLE='$E$',XTITLE='$x$');
fig_hw4.POSITION=[0.58,0.105,0.99,0.295]
fig_hw4.YRANGE=[100,450]
fig_hw4=PLOT(x,U[*,nt-1,1],/OVERPLOT,'--')
fig_hw4.SYMBOL='o'
t1=TEXT(315,497,'Time=0.01  /Lax-Wendroff  Nx='+STRMID(STRING(nx),5,3)+'  CFL='+STRMID(STRING(CFL),5,4),/DEVICE,ALIGNMENT=0.5)
fig_hw4.SAVE,'hw4_lw_1f.pdf',RESOLUTION=512,/TRANSPARENT
fig_hw4.CLOSE


;------------------slow shock-------------
Ws=MAKE_ARRAY(nx,nt,7) ;slow shock
Ws[0:count-1,0,*]=TRANSPOSE(REBIN([2.219, 0.4442, 0.5048, 0.0961, 0.0961, 1/SQRT(4*!PI), 1/SQRT(4*!PI)],7,count))
Ws[count:nx-1,0,*]=TRANSPOSE(REBIN([1, 0.1,-0.9225, 0, 0, 1/SQRT(4*!PI), 1/SQRT(4*!PI)],7,nx-count))



END