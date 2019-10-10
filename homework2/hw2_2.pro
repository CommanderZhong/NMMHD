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

num=61
np=5 ;number of pictures
dx=3./(num-1)
x=INDGEN(num)*dx-1.
u=FLTARR(num,np)
u[WHERE((x LT -0.8) OR (x GE 0.0)),0]=1.8
u[WHERE((x GE -0.8) AND (x LT -0.3)),0]=1.4+0.4*COS(2*!PI*(x[WHERE((x GE -0.8) AND (x LT -0.3))]+0.8))
u[WHERE((x GE -0.3) AND (x LT 0.0)),0]=1.0
;plot,x,u[*,0]

END