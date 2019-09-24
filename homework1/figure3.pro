PRO hw1_3
  th0=15*!DTOR  ;theta_0 with rad
  g=5./3. ;gama
  s1=1.4 ;1-g*SIN(th0)^2/(g-1) ~ 0.833
  s2=0.3
  n=201 ;number of data
  d=4e-4
  hf1=(2./(g-1))*SIN(th0)
  
  ;case 1
  s0=s1
  
  hf=INDGEN(n)/(n-1.)*hf1-d
  B=g/2*hf*SIN(th0)-(1-s0)
  C=2*SIN(th0)-(g-1)*hf
  Rx=B^2+C*(hf+2*s0*SIN(th0))
  
  y1=(B+SQRT(Rx))/C ;Xf+/hf
  out=PLOT(hf,y1,XRANGE=[0,1.1],YRANGE=[-MAX(y1/200.),MAX(y1)/60.])
  out.XTICKFORMAT='(A6)'
  out.YTICKFORMAT='(A6)'
  out.AXIS_STYLE=1
  out.XTICKLEN=0
  out.YTICKLEN=0
  out.THICK=2
  
  ;case 2
  s0=s2
  hf2=(SIN(th0)*(2-g)*(1+s0)+2*COS(th0)*SQRT((g-1)*(1-s0)^2+s0*g^2*SIN(th0)^2))/(2*(g-1)-g^2*SIN(th0)^2/2)
  y0=(1+s0*(g-1))*SIN(th0)/((1-s0)*(g-1)-g*SIN(th0)^2)
  
  hf=hf1+INDGEN(n)/(n-1.)*(hf2-hf1)+d
  B=g/2*hf*SIN(th0)-(1-s0)
  C=2*SIN(th0)-(g-1)*hf
  Rx=B^2+C*(hf+2*s0*SIN(th0))
  y2=(B-SQRT(Rx))/C ;Xf-/hf
  out=PLOT(hf,y2,/OVERPLOT,/CURR,'--')
  out.THICK=2
  
  hf=INDGEN(n)/(n-1.)*hf2
  B=g/2*hf*SIN(th0)-(1-s0)
  C=2*SIN(th0)-(g-1)*hf
  Rx=B^2+C*(hf+2*s0*SIN(th0))
  y3=(B+SQRT(Rx))/C ;Xf+/hf
  out=PLOT(hf,y3,/OVERPLOT,/CURR,'--')
  out.THICK=2
  
  line1y=INDGEN(n)/(n-1.)*(MAX(y1)/10.+50)-50
  line1x=REPLICATE(hf1,n)
  line2x=REPLICATE(hf2,n)
  line3x=INDGEN(n)/(n-1.)*1.2
  line3y=REPLICATE(y0,n)
  out=PLOT(line1x,line1y,/OVERPLOT,/CURR,'.')
  out=PLOT(line2x,line1y,/OVERPLOT,/CURR,'.')
  out=PLOT(line3x,line3y,/OVERPLOT,/CURR,'.')
;  t1=TEXT(384,355,'$\frac{X^+_f}{h_f}$',/DEVICE)
;  t1=TEXT(40,310,'$\frac{X^{\pm}_f}{h_f}$',/DEVICE)
;  t1=TEXT(440,155,'$\frac{X^+_f}{h_f}$',/DEVICE)
;  t1=TEXT(470,200,'$\frac{X^+_f}{h_f}$',/DEVICE)
;  t2=TEXT(300,40,'$h_f$',/DEVICE)
;  t2=TEXT(310,40,'$\hat{h_f}$',/DEVICE)
;  t2=TEXT(320,40,'$\hat{\hat{h_f}}$',/DEVICE)
  
  out.SAVE,'figure3.pdf',resolution=512,/TRANSPARENT
  out.CLOSE
END