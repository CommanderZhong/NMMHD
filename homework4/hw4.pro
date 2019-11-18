PRO hw4
;----------intial value----------
nx=201
gm=5./3 ;gamma
mu=1 ;mu
Hx=5 ;Hx or Bx
x0=0.
c=0.1 ;courant coeffitient

dx=4./(nx-1) ;-3---1
dt=dx*c
nt=0.1/dt+1
x=INDGEN(nx)/(nx-1.)*4-3
loc=WHERE(x LE 0,count)

;------------------slow shock-------------
Ws=MAKE_ARRAY(nx,nt,7,/DOUBLE) 
Ws[0:count-1,0,*]=TRANSPOSE(REBIN([2.121,4.981,-13.27,-0.163,-0.6521, 2.572, 10.29],7,count))
Ws[count:nx-1,0,*]=TRANSPOSE(REBIN([1, 1,-15.3, 0, 0, 1, 4],7,nx-count))
U=MAKE_ARRAY(nx,nt,7,/DOUBLE)
U[*,*,0]=Ws[*,*,0] & U[*,*,2:4]=Ws[*,*,2:4]*Ws[*,*,0] & U[*,*,5:6]=Ws[*,*,5:6]
bt=4*!PI*Ws[*,*,1]/TOTAL(Ws[*,0,5:6]^2,2)
U[*,*,1]=0.5*Ws[*,*,0]*TOTAL(Ws[*,*,2:4]^2,3)+0.5*TOTAL(Ws[*,*,5:6]^2,3)+bt*Ws[*,*,1]/(gm-1)
FOR n=0,nt-2 DO BEGIN
  U[0,n+1,*]=U[0,n,*]
  U[nx-1,n+1,*]=U[nx-1,n+1,*]
  FOR j=1,nx-2 DO BEGIN
    U[j,n+1,*]=U[j,n,*]-c*(Ftvd(U[*,n,*],gm,c,nx,j)-Ftvd(U[*,n,*],gm,c,nx,j-1))
  ENDFOR
ENDFOR

;------------------fast shock-------------
Wf=MAKE_ARRAY(nx,nt,7) ;fast shock
Wf[0:count-1,0,*]=TRANSPOSE(REBIN([2.219, 0.4442, 0.5048, 0.0961, 0.0961, 1, 1],7,count))
Wf[count:nx-1,0,*]=TRANSPOSE(REBIN([1, 0.1,-0.9225, 0, 0, 1, 1],7,nx-count))

END