FUNCTION U2W,U,gm
  Hx=5./SQRT(4*!PI)
  bt=2.
  W=MAKE_ARRAY(7,/DOUBLE)
  W[0]=U[0]
  W[5:6]=U[5:6]
  W[2:4]=U[2:4]/U[0]
  W[1]=(U[1]-W[0]*TOTAL(W[2:4]^2)-TOTAL(U[5:6]^2))*(gm-1)/bt
  RETURN,W
END

FUNCTION U2F,U,gm
  Hx=5./SQRT(4*!PI)
  bt=2
  F=MAKE_ARRAY(7,/DOUBLE)
  W=MAKE_ARRAY(7,/DOUBLE)
  W=U2W(U,gm)
  F[0]=W[0]*W[2]
  F[1]=W[0]*W[2]*(TOTAL(W[2:4]^2)+gm/(gm-1)*bt*W[1]/W[0])+$
  2*(W[5]^2*W[2]+W[6]^2*W[2]-Hx*W[5]*W[3]-Hx*W[6]*W[4])
  F[2]=W[0]*W[2]^2+bt*0.5*W[1]+0.5*(W[5]^2+W[6]^2)
  F[3]=W[0]*W[2]*W[3]-W[5]*Hx
  F[4]=W[0]*W[2]*W[4]-W[6]*Hx
  F[5]=W[2]*W[5]-W[3]*Hx
  F[6]=W[2]*W[6]-W[4]*Hx
  loo=where(~finite(F) eq 1,count)
  if count gt 1 then begin
    print,f
    stop
  endif
  RETURN,F
END