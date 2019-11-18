PRO hw4
;----------intial value----------
nx=401
gm=5./3 ;gamma
mu=1 ;mu
Hx=5 ;Hx or Bx
x0=0.
c=0.5 ;courant coeffitient

dx=4./(nx-1) ;-3---1
dt=dx*c
nt=0.1/dt+1
x=INDGEN(nx)/(nx-1.)*4-3
loc=WHERE(x LE 0,count)

;------------------fast shock-------------
Wf=MAKE_ARRAY(nx,nt,7,/DOUBLE) 
Wf[0:count-1,0,*]=TRANSPOSE(REBIN([2.121,4.981,-13.27,-0.163,-0.6521, 2.572, 10.29],7,count))
Wf[count:nx-1,0,*]=TRANSPOSE(REBIN([1, 1,-15.3, 0, 0, 1, 4],7,nx-count))
U=MAKE_ARRAY(nx,nt,7,/DOUBLE)
U[*,0,0]=Wf[*,0,0] &  U[*,0,5:6]=Wf[*,0,5:6]
FOR j=0,nx-1 DO BEGIN
  U[j,0,2:4]=Wf[j,0,2:4]*Wf[j,0,0]
ENDFOR
bt=4*!PI*Wf[*,0,1]/(TOTAL(REFORM(Wf[*,0,5:6]^2),2)+25)
U[*,0,1]=0.5*Wf[*,0,0]*TOTAL(Wf[*,0,2:4]^2,3)+0.5*TOTAL(Wf[*,0,5:6]^2,3)+bt*Wf[*,0,1]/(gm-1)
FOR n=0,nt-2 DO BEGIN
  U[0,n+1,*]=U[0,n,*]
  U[nx-1,n+1,*]=U[nx-1,n,*]
  Wf[0,n+1,*]=U2W(REFORM(U[0,n+1,*]),gm) & Wf[nx-1,n+1,*]=U2W(REFORM(U[nx-1,n+1,*]),gm)
  FOR j=1,nx-2 DO BEGIN
    ;U[j,n+1,*]=U[j,n,*]-c*(Ftvd(REFORM(U[*,n,*]),gm,c,nx,j)-Ftvd(REFORM(U[*,n,*]),gm,c,nx,j-1))
    U[j,n+1,*]=(U[j+1,n,*]+U[j-1,n,*])*0.5-c*0.5*(U2F(REFORM(U[j+1,n,*]),gm)-U2F(REFORM(U[j-1,n,*]),gm))
    Wf[j,n+1,*]=U2W(REFORM(U[j,n+1,*]),gm)
  ENDFOR
ENDFOR

;------------------slow shock-------------
Ws=MAKE_ARRAY(nx,nt,7) ;slow shock
Ws[0:count-1,0,*]=TRANSPOSE(REBIN([2.219, 0.4442, 0.5048, 0.0961, 0.0961, 1, 1],7,count))
Ws[count:nx-1,0,*]=TRANSPOSE(REBIN([1, 0.1,-0.9225, 0, 0, 1, 1],7,nx-count))


;------------------plot image------------
fig_hw4=PLOT(x,Wf[*,0,0],XTICKFORMAT='(A6)',':')
fig_hw4.YRANGE=[0.5,2.5]
fig_hw4=PLOT(x,Wf[*,nt-1,0],/OVERPLOT,'--')
fig_hw4.SYMBOL='o'

END