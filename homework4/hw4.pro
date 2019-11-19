PRO hw4,fs,ws,nx,cfl
;; keyword fs to judge the fast or slow shock
; fs=1 fast, fs=2 slow
;; keyword ws to judge the weak or strong shock
; ws=1 weak, ws=2 strong

IF NOT KEYWORD_SET(fs) THEN fs=1
IF NOT KEYWORD_SET(ws) THEN ws=1
IF NOT KEYWORD_SET(nx) THEN nx=201
IF NOT KEYWORD_SET(cfl) THEN cfl=0.8

;----------intial value----------
gm=5./3 ;gamma
mu=1 ;mu
Hx=5./SQRT(4*!PI) ;Hx or Bx
bt=2 ;beta
 
IF ws EQ 1 THEN BEGIN
  dx=4./(nx-1) ;-3---1
  x=INDGEN(nx)/(nx-1.)*4-3
  loc=WHERE(x LE 0,count)
  IF fs EQ 1 THEN BEGIN
    uv=11.4639
    c=CFL/uv ;courant coeffitient
    dt=dx*c
    nt=FLOOR(0.1/dt)+1
    W=MAKE_ARRAY(nx,nt,7,/DOUBLE)
    W[0:count-1,0,*]=TRANSPOSE(REBIN([2.121,4.981,-13.27,-0.163,-0.6521, 2.572/SQRT(4*!PI), 10.29/SQRT(4*!PI)],7,count))
    W[count:nx-1,0,*]=TRANSPOSE(REBIN([1, 1,-15.3, 0, 0, 1./SQRT(4*!PI), 4./SQRT(4*!PI)],7,nx-count))
  ENDIF ELSE BEGIN
    uv=0.1621
    c=CFL/uv ;courant coeffitient
    dt=dx*c
    nt=FLOOR(1.0/dt)+1
    W=MAKE_ARRAY(nx,nt,7,/DOUBLE)
    W[0:count-1,0,*]=TRANSPOSE(REBIN([2.219, 0.4442, 0.5048, 0.0961, 0.0961, 1/SQRT(4*!PI), 1/SQRT(4*!PI)],7,count))
    W[count:nx-1,0,*]=TRANSPOSE(REBIN([1, 0.1,-0.9225, 0, 0, 1/SQRT(4*!PI), 1/SQRT(4*!PI)],7,nx-count))
  ENDELSE
ENDIF
IF ws EQ 2 THEN BEGIN
  dx=1./(nx-1) ;-3---1
  x=INDGEN(nx)/(nx-1.)
  loc=WHERE(x LE 0.2,count)
  IF fs EQ 1 THEN BEGIN
    uv=-5.2049
    c=CFL/ABS(uv) ;courant coeffitient
    dt=dx*c
    nt=FLOOR(0.1/dt)+1
    nt1=FLOOR(0.05/dt)+1
    W=MAKE_ARRAY(nx,nt,7,/DOUBLE)
    W[0:count-1,0,*]=TRANSPOSE(REBIN([3.896,305.9,0,-0.058,-0.226,3.951/SQRT(4*!PI),15.8/SQRT(4*!PI)],7,count))
    W[count:nx-1,0,*]=TRANSPOSE(REBIN([1, 1,-15.3, 0, 0, 1./SQRT(4*!PI), 4./SQRT(4*!PI)],7,nx-count))
  ENDIF ELSE BEGIN
    uv=-0.4377
    c=CFL/ABS(uv) ;courant coeffitient
    dt=dx*c
    nt=FLOOR(1.6/dt)+1
    nt1=FLOOR(0.8/dt)+1
    W=MAKE_ARRAY(nx,nt,7,/DOUBLE)
    W[0:count-1,0,*]=TRANSPOSE(REBIN([3.108, 1.4336,0,0.2633,0.2633, 0.1/SQRT(4*!PI), 0.1/SQRT(4*!PI)],7,count))
    W[count:nx-1,0,*]=TRANSPOSE(REBIN([1, 0.1,-0.9225, 0, 0, 1/SQRT(4*!PI), 1/SQRT(4*!PI)],7,nx-count))
  ENDELSE
ENDIF

U=MAKE_ARRAY(nx,nt,7,/DOUBLE)
U[*,0,0]=W[*,0,0] &  U[*,0,5:6]=W[*,0,5:6]
FOR j=0,nx-1 DO BEGIN
  U[j,0,2:4]=W[j,0,2:4]*W[j,0,0]
ENDFOR
U[*,0,1]=W[*,0,0]*TOTAL(REFORM(W[*,0,2:4]^2),2)+TOTAL(REFORM(W[*,0,5:6]^2),2)+bt*W[*,0,1]/(gm-1)
FOR n=0,nt-2 DO BEGIN
  U[0,n+1,*]=U[0,n,*]
  U[nx-1,n+1,*]=U[nx-1,n,*]
  W[0,n+1,*]=U2W(REFORM(U[0,n+1,*]),gm) & W[nx-1,n+1,*]=U2W(REFORM(U[nx-1,n+1,*]),gm)
  FOR j=1,nx-2 DO BEGIN
;    U[j,n+1,*]=(U[j+1,n,*]+U[j-1,n,*])*0.5-c*0.5*(U2F(REFORM(U[j+1,n,*]),gm)-U2F(REFORM(U[j-1,n,*]),gm)) ;Lax
    Um=0.5*(U[j,n,*]+U[j+1,n,*])-0.5*c*(U2F(REFORM(U[j+1,n,*]),gm)-U2F(REFORM(U[j,n,*]),gm))
    Un=0.5*(U[j,n,*]+U[j-1,n,*])-0.5*c*(U2F(REFORM(U[j,n,*]),gm)-U2F(REFORM(U[j-1,n,*]),gm))
    U[j,n+1,*]=U[j,n,*]-c*(U2F(REFORM(Um),gm)-U2F(REFORM(Un),gm)) ;Lax-Wendroff
    W[j,n+1,*]=U2W(REFORM(U[j,n+1,*]),gm)
  ENDFOR
ENDFOR

;------------------plot image------------
fig_hw4=PLOT(x,W[*,0,0],XTICKFORMAT='(A6)',':',YTITLE='$\rho$');
fig_hw4.POSITION=[0.1,0.765,0.47,0.95]
fig_hw4.YRANGE=[MIN(W[*,0,0])-0.3*MAX(ABS(W[*,0,0])),MAX(W[*,0,0])+0.3*MAX(ABS(W[*,0,0]))]
fig_hw4=PLOT(x,W[*,nt-1,0],/OVERPLOT,'--')
fig_hw4.SYMBOL='o'
fig_hw4=PLOT(x-(nt-1)*dt*uv,W[*,0,0],/OVERPLOT,XRANGE=[MIN(x),MAX(x)])
IF ws EQ 2 THEN BEGIN
  fig_hw4=PLOT(x,W[*,nt1-1,0],/OVERPLOT,'--')
  fig_hw4.SYMBOL='x'
  fig_hw4=PLOT(x-(nt1-1)*dt*uv,W[*,0,0],/OVERPLOT,XRANGE=[MIN(x),MAX(x)])
ENDIF
fig_hw4=PLOT(x,W[*,0,1],/CURR,XTICKFORMAT='(A6)',':',YTITLE='$p$');
fig_hw4.POSITION=[0.1,0.545,0.47,0.735]
fig_hw4.YRANGE=[MIN(W[*,0,1])-0.3*MAX(ABS(W[*,0,1])),MAX(W[*,0,1])+0.3*MAX(ABS(W[*,0,1]))]
fig_hw4=PLOT(x,W[*,nt-1,1],/OVERPLOT,'--')
fig_hw4.SYMBOL='o'
fig_hw4=PLOT(x-(nt-1)*dt*uv,W[*,0,1],/OVERPLOT,XRANGE=[MIN(x),MAX(x)])
IF ws EQ 2 THEN BEGIN
  fig_hw4=PLOT(x,W[*,nt1-1,1],/OVERPLOT,'--')
  fig_hw4.SYMBOL='x'
  fig_hw4=PLOT(x-(nt1-1)*dt*uv,W[*,0,1],/OVERPLOT,XRANGE=[MIN(x),MAX(x)])
ENDIF
fig_hw4=PLOT(x,W[*,0,2],/CURR,XTICKFORMAT='(A6)',':',YTITLE='$v_x$');
fig_hw4.POSITION=[0.1,0.325,0.47,0.515]
fig_hw4.YRANGE=[MIN(W[*,0,2])-0.3*MAX(ABS(W[*,0,2])),MAX(W[*,0,2])+0.3*MAX(ABS(W[*,0,2]))]
fig_hw4=PLOT(x,W[*,nt-1,2],/OVERPLOT,'--')
fig_hw4.SYMBOL='o'
fig_hw4=PLOT(x-(nt-1)*dt*uv,W[*,0,2],/OVERPLOT,XRANGE=[MIN(x),MAX(x)])
IF ws EQ 2 THEN BEGIN
  fig_hw4=PLOT(x,W[*,nt1-1,2],/OVERPLOT,'--')
  fig_hw4.SYMBOL='x'
  fig_hw4=PLOT(x-(nt1-1)*dt*uv,W[*,0,2],/OVERPLOT,XRANGE=[MIN(x),MAX(x)])
ENDIF
fig_hw4=PLOT(x,W[*,0,3],/CURR,':',YTITLE='$v_y$',XTITLE='$x$');
fig_hw4.POSITION=[0.1,0.105,0.47,0.295]
fig_hw4.YRANGE=[MIN(W[*,0,3])-0.3*MAX(ABS(W[*,0,3])),MAX(W[*,0,3])+0.3*MAX(ABS(W[*,0,3]))]
fig_hw4=PLOT(x,W[*,nt-1,3],/OVERPLOT,'--')
fig_hw4.SYMBOL='o'
fig_hw4=PLOT(x-(nt-1)*dt*uv,W[*,0,3],/OVERPLOT,XRANGE=[MIN(x),MAX(x)])
IF ws EQ 2 THEN BEGIN
  fig_hw4=PLOT(x,W[*,nt1-1,3],/OVERPLOT,'--')
  fig_hw4.SYMBOL='x'
  fig_hw4=PLOT(x-(nt1-1)*dt*uv,W[*,0,3],/OVERPLOT,XRANGE=[MIN(x),MAX(x)])
ENDIF
fig_hw4=PLOT(x,W[*,0,4],/CURR,XTICKFORMAT='(A6)',':',YTITLE='$v_z$');
fig_hw4.POSITION=[0.58,0.765,0.99,0.95]
fig_hw4.YRANGE=[MIN(W[*,0,4])-0.3*MAX(ABS(W[*,0,4])),MAX(W[*,0,4])+0.3*MAX(ABS(W[*,0,4]))]
fig_hw4=PLOT(x,W[*,nt-1,4],/OVERPLOT,'--')
fig_hw4.SYMBOL='o'
fig_hw4=PLOT(x-(nt-1)*dt*uv,W[*,0,4],/OVERPLOT,XRANGE=[MIN(x),MAX(x)])
IF ws EQ 2 THEN BEGIN
  fig_hw4=PLOT(x,W[*,nt1-1,4],/OVERPLOT,'--')
  fig_hw4.SYMBOL='x'
  fig_hw4=PLOT(x-(nt1-1)*dt*uv,W[*,0,4],/OVERPLOT,XRANGE=[MIN(x),MAX(x)])
ENDIF
fig_hw4=PLOT(x,W[*,0,5],/CURR,XTICKFORMAT='(A6)',':',YTITLE='$H_y$');
fig_hw4.POSITION=[0.58,0.545,0.99,0.735]
fig_hw4.YRANGE=[MIN(W[*,0,5])-0.3*MAX(ABS(W[*,0,5])),MAX(W[*,0,5])+0.3*MAX(ABS(W[*,0,5]))]
fig_hw4=PLOT(x,W[*,nt-1,5],/OVERPLOT,'--')
fig_hw4.SYMBOL='o'
fig_hw4=PLOT(x-(nt-1)*dt*uv,W[*,0,5],/OVERPLOT,XRANGE=[MIN(x),MAX(x)])
IF ws EQ 2 THEN BEGIN
  fig_hw4=PLOT(x,W[*,nt1-1,5],/OVERPLOT,'--')
  fig_hw4.SYMBOL='x'
  fig_hw4=PLOT(x-(nt1-1)*dt*uv,W[*,0,5],/OVERPLOT,XRANGE=[MIN(x),MAX(x)])
ENDIF
fig_hw4=PLOT(x,W[*,0,6],/CURR,XTICKFORMAT='(A6)',':',YTITLE='$H_z$');
fig_hw4.POSITION=[0.58,0.325,0.99,0.515]
fig_hw4.YRANGE=[MIN(W[*,0,6])-0.3*MAX(ABS(W[*,0,6])),MAX(W[*,0,6])+0.3*MAX(ABS(W[*,0,6]))]
fig_hw4=PLOT(x,W[*,nt-1,6],/OVERPLOT,'--')
fig_hw4.SYMBOL='o'
fig_hw4=PLOT(x-(nt-1)*dt*uv,W[*,0,6],/OVERPLOT,XRANGE=[MIN(x),MAX(x)])
IF ws EQ 2 THEN BEGIN
  fig_hw4=PLOT(x,W[*,nt1-1,6],/OVERPLOT,'--')
  fig_hw4.SYMBOL='x'
  fig_hw4=PLOT(x-(nt1-1)*dt*uv,W[*,0,6],/OVERPLOT,XRANGE=[MIN(x),MAX(x)])
ENDIF
fig_hw4=PLOT(x,U[*,0,1],/CURR,':',YTITLE='$E$',XTITLE='$x$');
fig_hw4.POSITION=[0.58,0.105,0.99,0.295]
fig_hw4.YRANGE=[MIN(U[*,0,1])-0.3*MAX(ABS(U[*,0,1])),MAX(U[*,0,1])+0.3*MAX(ABS(U[*,0,1]))]
fig_hw4=PLOT(x,U[*,nt-1,1],/OVERPLOT,'--')
fig_hw4.SYMBOL='o'
fig_hw4=PLOT(x-(nt-1)*dt*uv,U[*,0,1],/OVERPLOT,XRANGE=[MIN(x),MAX(x)])
IF ws EQ 2 THEN BEGIN
  fig_hw4=PLOT(x,U[*,nt1-1,1],/OVERPLOT,'--')
  fig_hw4.SYMBOL='x'
  fig_hw4=PLOT(x-(nt1-1)*dt*uv,U[*,0,1],/OVERPLOT,XRANGE=[MIN(x),MAX(x)])
ENDIF
t1=TEXT(315,497,'Nx='+STRMID(STRING(nx),5,3)+'  CFL='+STRMID(STRING(CFL),5,4)+'  /Lax-Wendroff',/DEVICE,ALIGNMENT=0.5)
IF (fs EQ 1) AND (ws EQ 1) AND (nx=131) THEN fig_hw4.SAVE,'hw4_lw_1f1.pdf',RESOLUTION=512,/TRANSPARENT
IF (fs EQ 1) AND (ws EQ 1) AND (nx=261) THEN fig_hw4.SAVE,'hw4_lw_1f2.pdf',RESOLUTION=512,/TRANSPARENT
IF (fs EQ 2) AND (ws EQ 1) THEN fig_hw4.SAVE,'hw4_lw_1s.pdf',RESOLUTION=512,/TRANSPARENT
IF (fs EQ 1) AND (ws EQ 2) THEN fig_hw4.SAVE,'hw4_lw_2f.pdf',RESOLUTION=512,/TRANSPARENT
IF (fs EQ 2) AND (ws EQ 2) THEN fig_hw4.SAVE,'hw4_lw_2s.pdf',RESOLUTION=512,/TRANSPARENT
fig_hw4.CLOSE
END

PRO use_hw4
;;To use procedure hw4

  hw4,1,1,133,0.65
  hw4,1,1,261,0.65
; hw4,2,1,201,0.8  ;bad result
  hw4,1,2,133,0.3
  hw4,2,2,133,0.2
END