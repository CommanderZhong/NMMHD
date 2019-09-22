PRO dispersion,ita,xmax,ymax		;ita=omep/omec
;parameters setting 

ome=(INDGEN(501)+1)*ymax/500		;ome/omec			
r=1-ita^2/(ome^2+ome)
l=1-ita^2/(ome^2-ome)
p=1-ita^2/ome^2
s=0.5*(r+l)
nc=INDGEN(501)*xmax/500
SET_PLOT,'ps'
DEVICE,FILENAME='figure2_'+NUM2STR(ita)+'.eps',/COLOR,/ENCAPSULATED

tansq=FLTARR(501,501)
theta=tansq
;caculating theta

FOR i=0,500 DO BEGIN
	tansq[i,*]=-p*(nc[i]^2/ome^2-r)*(nc[i]^2/ome^2-l)/((s*nc[i]^2/ome^2-r*l)*(nc[i]^2/ome^2-p))
ENDFOR
post=WHERE(tansq LT 0)
theta=ATAN(SQRT(tansq))
theta[post]=-1

;plotting parameters setting and etc.
LOADCT,13
!P.THICK=0.2
CONTOUR,theta,nc,ome,NLEVELS=128,XRANGE=[0.0,xmax],YRANGE=[0.0,ymax],C_CHARTHICK=2,ZRANGE=[0,3.14159/2],XTITLE='!3kc/!4x!3!dc!n!x',YTITLE='!4x!3/!4x!3!dc!n!x',/FILL
CONTOUR,theta,nc,ome,NLEVELS=30,/OVERPLOT,C_CHARTHICK=2
OPLOT,nc,nc,LINESTYLE=1
XYOUTS,0.9*xmax,0.8*ymax,'!4h!3=90!uo!n!x',CHARTHICK=2
XYOUTS,0.9*xmax,0.6*ymax,'!4h!3=0!uo!n!x',CHARTHICK=2
XYOUTS,0.5*xmax,1.1*ymax,'!4x!3/!4x!3!dc!n!x='+NUM2STR(ita),CHARTHICK=2
DEVICE,/CLOSE

END
