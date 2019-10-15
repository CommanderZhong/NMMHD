PRO use_hw2_2
  t0=[0.25,0.5,0.75,1.0]
  FOR i=0,N_ELEMENTS(t0)-1 DO BEGIN
    hw2_2,t0[i],0.01
  ENDFOR
END