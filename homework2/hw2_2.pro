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
PRO hw2_2,t0,dx
;;t0=0.25,0.5,0.75,1.0

IF NOT KEYWORD_SET(t0) THEN t0=0.25
IF NOT KEYWORD_SET(dx) THEN dx=0.02

num=5./dx+1
x=INDGEN(num)*dx-3
u=FLTARR(num,num)
u[WHERE((x LT -0.8) OR (x GE 0.0)),0]=1.8
u[WHERE((x GE -0.8) AND (x LT -0.3)),0]=1.4+0.4*COS(2*!PI*(x[WHERE((x GE -0.8) AND (x LT -0.3))]+0.8))
u[WHERE((x GE -0.3) AND (x LT 0.0)),0]=1.0
c=0.95/MAX(u) ;dt<=dx/max(u)
dt=c*dx
nt=ROUND(1./dt)+1
t=INDGEN(nt)*dt
u1=u
c1=1./MAX(u)
dt1=c1*dx
nt1=ROUND(1./dt1)+1
t1=INDGEN(nt1)*dt1

FOR n=0,nt-2 DO BEGIN
  FOR j=0,num-1 DO BEGIN
    IF j EQ 0 THEN BEGIN
      u[j,n+1]=u[j,n]
    ENDIF ELSE BEGIN
      IF (j EQ 1) or (j EQ num-1) THEN BEGIN
        u[j,n+1]=u[j,n]-0.5*c*(u[j,n]^2-u[j-1,n]^2)
      ENDIF ELSE BEGIN
        u[j,n+1]=u[j,n]-0.5*c*(u[j,n]^2-u[j-1,n]^2)-0.25*c*(dx-MAX(u)*dt)*(minmod((u[j,n]-u[j-1,n])/dx,(u[j+1,n]-u[j,n])/dx)-minmod((u[j-1,n]-u[j-2,n])/dx,(u[j,n]-u[j-1,n])/dx))
      ENDELSE
    ENDELSE
  ENDFOR
ENDFOR
FOR n=0,nt1-2 DO BEGIN
  FOR j=0,num-1 DO BEGIN
    IF j EQ 0 THEN BEGIN
      u1[j,n+1]=u1[j,n]
    ENDIF ELSE BEGIN
      u1[j,n+1]=u1[j,n]-0.5*c1*(u1[j,n]^2-u1[j-1,n]^2)
    ENDELSE
  ENDFOR
ENDFOR

k=(t0/dt)
k1=(t0/dt1)

;A = FINDGEN(17) * (!PI*2/16.)
;USERSYM, 1.5*COS(A), 1.5*SIN(A)
;PLOT,x,u[*,0],YRANGE=[0.8,2.0],XRANGE=[-1.0,2.0],LINESTYLE=1,TITLE='Time='+STRMID(STRING(k*dt),5,6)
;FOR n=2./dx,N_ELEMENTS(u[*,k])-1 DO BEGIN
;  PLOTS,x[n],u[n,k],/DATA,PSYM=8
;ENDFOR
fig=PLOT(x,u[*,0],YRANGE=[0.8,2.0],XRANGE=[-1.0,2.0],':',TITLE='Time='+STRMID(STRING(k*dt),5,6)+'/  Minmod')
fig=PLOT(x,u[*,ROUND(k)],' ',/OVERPLOT,/CURR)
fig.SYMBOL='o'
fig.SYM_SIZE=1.5
fig=PLOT(x,u1[*,ROUND(k1)],/OVERPLOT,/CURR)
fig.SAVE,'fig2_'+STRMID(STRING(t0/0.25),6,1)+'.pdf',resolution=512,/TRANSPARENT
fig.CLOSE
END