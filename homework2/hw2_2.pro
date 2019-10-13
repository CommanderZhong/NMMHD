Function minmod,a,b
  IF (ABS(a) LE ABS(b)) AND (a*b GT 0.) THEN RETURN,a
  IF (ABS(a) GT ABS(b)) AND (a*b GT 0.) THEN RETURN,b
  IF a*b LE 0. THEN RETURN,0.
END

;;for t=0
;; u=  1.8, x<-0.8
;;     1.4+0.4*cos(2pi(x+0.8)), -0.8<=x<-0.3
;;     1.0, -0.3<=x<0.0
;;     1.8, x>=0.0
PRO hw2_2


c=1. ;c=dt/dx
dx=0.05
num=5./dx+1
dt=c*dx
x=INDGEN(num)*dx-3
t=INDGEN(num)*dt
u=FLTARR(num,num)
u[WHERE((x LT -0.8) OR (x GE 0.0)),0]=1.8
u[WHERE((x GE -0.8) AND (x LT -0.3)),0]=1.4+0.4*COS(2*!PI*(x[WHERE((x GE -0.8) AND (x LT -0.3))]+0.8))
u[WHERE((x GE -0.3) AND (x LT 0.0)),0]=1.0
;plot,x,u[*,0]
FOR n=0,2./dt-2 DO BEGIN
  FOR j=0,num-1 DO BEGIN
    IF j EQ 0 THEN BEGIN
      u[j,n+1]=u[j,n]
    ENDIF ELSE BEGIN
      u[j,n+1]=u[j,n]-c*(u[j,n]-u[j-1,n])
    ENDELSE
  ENDFOR
ENDFOR
plot,x,u[*,0]
oplot,x,u[*,29],linestyle=2
END