pro hw1_3
  th0=15*!DTOR  ;theta_0 with rad
  g=5./3. ;gama
  s1=1.4 ;1-g*SIN(th0)^2/(g-1) ~ 0.833
  s2=0.3
  n=401 ;number of data
  d=2e-4
  hf1=(2./(g-1))*SIN(th0)
  
  ;case 1
  s0=s1
  
  hf=INDGEN(n)/(n-1.)*hf1-d
  B=g/2*hf*SIN(th0)-(1-s0)
  C=2*SIN(th0)-(g-1)*hf
  Rx=B^2+C*(hf+2*s0*SIN(th0))
  
  y1=(B+SQRT(Rx))/C ;Xf+/hf
  out=PLOT(hf,y1,XRANGE=[0,MAX(hf)*1.6],YRANGE=[-MAX(y1/200.),MAX(y1)/50.])
  out.XTICKFORMAT='(A6)'
  out.YTICKFORMAT='(A6)'
  out.AXIS_STYLE=1
  out.XTICKLEN=0
  out.YTICKLEN=0
  
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
  
  hf=INDGEN(n)/(n-1.)*hf2-0.001
  B=g/2*hf*SIN(th0)-(1-s0)
  C=2*SIN(th0)-(g-1)*hf
  Rx=B^2+C*(hf+2*s0*SIN(th0))
  y3=(B+SQRT(Rx))/C ;Xf+/hf
  out=PLOT(hf,y3,/OVERPLOT,/CURR,'--')
end