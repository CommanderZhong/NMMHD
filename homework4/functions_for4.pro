FUNCTION Q,x,ep
;;Equation (9.62) of the lecture note

  Q0=x
  coor=WHERE(ABS(x) LT 2*ep,count)
  IF count GT 0 THEN BEGIN
    Q0[coor]=x[coor]^2/(4*ep)+ep
  ENDIF
  Q0[WHERE(ABS(x) GE 2*ep)]=ABS(x[WHERE(ABS(x) GE 2*ep)])
  RETURN,Q0
END

FUNCTION U2W,U,gm
  Bx=5
  
  W=MAKE_ARRAY(7,/DOUBLE)
  W[0]=U[0]
  W[5:6]=U[5:6]
  W[2:4]=U[2:4]/U[0]
  W[1]=SQRT((U[1]-0.5*TOTAL(U[2:4]^2)/U[0]-0.5*TOTAL(U[5:6]^2))*$
    (gm-1)*(TOTAL(U[5:6]^2)+Bx^2)/(4*!PI))
  RETURN,W
END

FUNCTION U2F,U,gm
  Bx=5
  
  F=MAKE_ARRAY(7,/DOUBLE)
  W=MAKE_ARRAY(7,/DOUBLE)
  W=U2W(U,gm)
  bt=4*!PI*W[1]/(TOTAL(W[5:6]^2)+Bx^2)
  F[0]=W[0]*W[2]
  F[1]=W[0]*W[2]*(0.5*TOTAL(W[2:4]^2)+gm/(gm-1)*bt*W[1]/W[0])+$
  (W[5]^2*W[2]+W[6]^2*W[2]-Bx*W[5]*W[3]-Bx*W[5]*W[4])
  F[2]=W[0]*W[2]^2+bt*W[1]+0.5*(W[5]^2+W[6]^2)
  F[3]=W[0]*W[2]*W[3]-W[5]*Bx
  F[4]=W[0]*W[2]*W[4]-W[6]*Bx
  F[5]=W[2]*W[5]-W[3]*Bx
  F[6]=W[2]*W[6]-W[4]*Bx
  RETURN,F
END

FUNCTION epyz,By,Bz
;first for epsilon_y, second for epsilon_z
  IF By^2+Bz^2 NE 0 THEN RETURN,[By/SQRT(By^2+Bz^2),Bz/SQRT(By^2+Bz^2)]
  IF By^2+Bz^2 EQ 0 THEN RETURN,[1/SQRT(2),1/SQRT(2)]
END

FUNCTION alsf,a,cf,cs
;first for alpha_s, second for alpha_f
  RETURN,[SQRT((cf^2-a^2)/(cf^2-cs^2)),SQRT((a^2-cs^2)/(cf^2-cs^2))]
END

FUNCTION Evelocity,W,gm
  Bx=5
  v=MAKE_ARRAY(5) 
  bt=4*!PI*W[1]/(TOTAL(W[5:6]^2)+Bx^2)
  a=SQRT(bt*gm*W[1]/W[0])
  bx=Bx/SQRT(W[0])
  bv=SQRT(TOTAL(W[5:6]^2)/W[0])
  b2=bx^2+bv^2
  cs=SQRT(0.5*(a^2+b2-SQRT((a^2+b2)^2-4*a^2*bx^2)))
  cf=SQRT(0.5*(a^2+b2+SQRT((a^2+b2)^2-4*a^2*bx^2)))
  vx=W[2]
  v=[a,bx,cs,cf,vx]
  RETURN,v
END

FUNCTION LeftEV,W,gm
  
END
