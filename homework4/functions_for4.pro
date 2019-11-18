FUNCTION Q,x
;;Equation (9.62) of the lecture note
  ep=0.1
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
  W[0]=U(0)
  W[5:6]=U(5:6)
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
  v=[a,bx,cs,cf,vx,bv]
  RETURN,v
END

FUNCTION LeftEV,W,gm
  Bx=5
  v=MAKE_ARRAY(6,/DOUBLE) ;a,bx,cs,cf,vx,bv
  v=Evelocity(W,gm)
  al=alsf(v[0],v[3],v[2]) ;alpha s $ f
  ep=epyz(W[5],W[6])  ;epsilon y $ z
  L=MAKE_ARRAY(7,7,/DOUBLE)
  
  ;for lambda=vx-+cf
  L[*,0]=[(gm-1)*0.5*al[1]*TOTAL(W[2:4]^2)+al[1]*v[3]*W[2]-al[0]*v[2]*(ep[0]*W[3]+ep[1]*W[4]),$
  (gm-1)*0.5*al[1],al[1]*((1-gm)*W[2]-v[3]),(1-gm)*al[1]*W[3]+al[0]*v[2]*ep[0],$
  (1-gm)*al[1]*W[4]+al[0]*v[2]*ep[1],SQRT(W[0])*ep[0]*((1-gm)*al[1]*v[5]+al[0]*v[0]),$
  SQRT(W[0])*ep[0]*((1-gm)*al[1]*v[5]+al[0]*v[0])]/(W[0]*v[0]^2)
  
  L[*,6]=[(gm-1)*0.5*al[1]*TOTAL(W[2:4]^2)-al[1]*v[3]*W[2]+al[0]*v[2]*(ep[0]*W[3]+ep[1]*W[4]),$
  (gm-1)*0.5*al[1],al[1]*((1-gm)*W[2]+v[3]),(1-gm)*al[1]*W[3]-al[0]*v[2]*ep[0],$
  (1-gm)*al[1]*W[4]-al[0]*v[2]*ep[1],SQRT(W[0])*ep[0]*((1-gm)*al[1]*v[5]+al[0]*v[0]),$
  SQRT(W[0])*ep[0]*((1-gm)*al[1]*v[5]+al[0]*v[0])]/(W[0]*v[0]^2)
  
  ;for lambda=vx-+cs
  L[*,2]=[(gm-1)*0.5*al[0]*TOTAL(W[2:4]^2)+al[0]*v[2]*W[2]+al[1]*v[3]*(ep[0]*W[3]+ep[1]*W[4]),$
  (gm-1)*0.5*al[0],al[0]*((1-gm)*W[2]-v[2]),(1-gm)*al[0]*W[3]-al[1]*v[3]*ep[0],$
  (1-gm)*al[0]*W[4]-al[1]*v[3]*ep[1],SQRT(W[0])*ep[0]*((1-gm)*al[0]*v[5]-al[1]*v[0]),$
  SQRT(W[0])*ep[1]*((1-gm)*al[0]*v[5]-al[1]*v[0])]/(W[0]*v[0]^2)
  
  L[*,4]=[(gm-1)*0.5*al[0]*TOTAL(W[2:4]^2)-al[0]*v[2]*W[2]-al[1]*v[3]*(ep[0]*W[3]+ep[1]*W[4]),$
  (gm-1)*0.5*al[0],al[0]*((1-gm)*W[2]+v[2]),(1-gm)*al[0]*W[3]+al[1]*v[3]*ep[0],$
  (1-gm)*al[0]*W[4]+al[1]*v[3]*ep[1],SQRT(W[0])*ep[0]*((1-gm)*al[0]*v[5]-al[1]*v[0]),$
  SQRT(W[0])*ep[1]*((1-gm)*al[0]*v[5]-al[1]*v[0])]/(W[0]*v[0]^2)
  
  ;for lambda=vx-+bx
  L[*,1]=[0,0,0,ep[1]/SQRT(2),-ep[0]/SQRT(2),ep(1)/SQRT(2)*SQRT(W[0]),-ep[0]/SQRT(2)*SQRT(W[0])]/W[0]
  L[*,5]=[0,0,0,ep[1]/SQRT(2),-ep[0]/SQRT(2),-ep(1)/SQRT(2)*SQRT(W[0]),+ep[0]/SQRT(2)*SQRT(W[0])]/W[0]
  
  ;for lambda=vx
  L[*,3]=(gm-1)/v[0]^2*[v[0]^2/(gm-1)-0.5*TOTAL(W[2:4]^2),-0.5,W[2],W[3],W[4],SQRT(W[0])*v[5]*ep[0],SQRT(W[0])*v[5]*ep[1]]
  
  RETURN,TRANSPOSE(L)
END

FUNCTION RightEV,W,gm
  Bx=5
  v=MAKE_ARRAY(6,/DOUBLE) ;a,bx,cs,cf,vx,bv
  v=Evelocity(W,gm)
  al=alsf(v[0],v[3],v[2]) ;alpha s $ f
  ep=epyz(W[5],W[6])  ;epsilon y $ z
  R=MAKE_ARRAY(7,7,/DOUBLE)
  
  ;for lambda=-+cf
  R[0,*]=W[0]*[al[1]*0.5,0.5*al[1]*TOTAL(W[2:4])+al[1]*v[0]^2/(gm-1)+al[0]*v[0]*v[5]-al[1]*W[2]*v[3]+al[0]*v[2]*(ep[0]*W[3]+ep[1]*W[4]),$
    al[1]*(W[2]-v[3])*0.5,0.5*(al[1]*W[3]+al[0]*v[2]*ep[0]),$
    0.5*(al[1]*W[4]+al[0]*v[2]*ep[1]),0.5*v[0]*al[0]*ep[0]/SQRT(W[0]),0.5*v[0]*al[0]*ep[1]/SQRT(W[0])]
    
  R[6,*]=W[0]*[al[1]*0.5,0.5*al[1]*TOTAL(W[2:4])+al[1]*v[0]^2/(gm-1)+al[0]*v[0]*v[5]+al[1]*W[2]*v[3]-al[0]*v[2]*(ep[0]*W[3]+ep[1]*W[4]),$
    al[1]*(W[2]+v[3])*0.5,0.5*(al[1]*W[3]-al[0]*v[2]*ep[0]),$
    0.5*(al[1]*W[4]-al[0]*v[2]*ep[1]),0.5*v[0]*al[0]*ep[0]/SQRT(W[0]),0.5*v[0]*al[0]*ep[1]/SQRT(W[0])]
    
  ;for lambda=vx-+cs
  R[2,*]=W[0]*[al[0]*0.5,0.5*al[0]*TOTAL(W[2:4])+al[0]*v[0]^2/(gm-1)-al[1]*v[0]*v[5]-al[0]*W[2]*v[2]-al[1]*v[3]*(ep[0]*W[3]+ep[1]*W[4]),$
    al[0]*(W[2]-v[2])*0.5,0.5*(al[0]*W[3]-al[1]*v[3]*ep[0]),$
    0.5*(al[0]*W[4]-al[1]*v[3]*ep[1]),-0.5*v[0]*al[1]*ep[0]/SQRT(W[0]),-0.5*v[0]*al[1]*ep[1]/SQRT(W[0])]

  R[4,*]=W[0]*[al[0]*0.5,0.5*al[0]*TOTAL(W[2:4])+al[0]*v[0]^2/(gm-1)-al[1]*v[0]*v[5]+al[0]*W[2]*v[2]+al[1]*v[3]*(ep[0]*W[3]+ep[1]*W[4]),$
    al[0]*(W[2]+v[2])*0.5,0.5*(al[0]*W[3]+al[1]*v[3]*ep[0]),$
    0.5*(al[0]*W[4]+al[1]*v[3]*ep[1]),-0.5*v[0]*al[1]*ep[0]/SQRT(W[0]),-0.5*v[0]*al[1]*ep[1]/SQRT(W[0])]
    
  ;for lambda=vx-+bx
  R[1,*]=W[0]*[0,0,0,ep[1]/SQRT(2),-ep[0]/SQRT(2),ep[1]/SQRT(2*W[0]),-ep[0]/SQRT(2*W[0])]
  R[5,*]=W[0]*[0,0,0,ep[1]/SQRT(2),-ep[0]/SQRT(2),-ep[1]/SQRT(2*W[0]),ep[0]/SQRT(2*W[0])]
  
  ;for lambda=vx
  R[3,*]=[1,TOTAL(W[2:4]^2),W[2],W[3],W[4],0,0]
  
  RETURN,TRANSPOSE(R)
END

FUNCTION gm1_2,alpha,g2,g1
  gmk=alpha
  coor=WHERE(ABS(alpha) NE 0,count)
  IF count GT 0 THEN BEGIN
    gmk[coor]=(g2[coor]-g1[coor])/alpha[coor]
  ENDIF
  gmk[WHERE(alpha EQ 0)]=0
  RETURN,gmk
END

FUNCTION Ftvd,U,gm,c,nx,j
  ;U [nx,7]
  L1=MAKE_ARRAY(7,7,/DOUBLE) & L2=L1 & L3=L1 ;1 for j-1/2, 2 for j+1/2, 3 for j+3/2
  R1=MAKE_ARRAY(7,7,/DOUBLE) & R2=R1 & R3=R1 ;1 for j-1/2, 2 for j+1/2, 3 for j+3/2
  alpha1=MAKE_ARRAY(7,/DOUBLE) & alpha2=alpha1 & alpha3=alpha1 ;1 for j-1/2, 2 for j+1/2, 3 for j+3/2
  s1=alpha1 & s2=s1 & s3=s1 ;1 for j-1/2, 2 for j+1/2, 3 for j+3/2
  gp1=s1 & gp2=s1 & gp3=s1 ;1 for j-1/2, 2 for j+1/2, 3 for j+3/2
  g1=s1 & g2=s1 ;1 for j, 2 for j+1
  g1_2=MAKE_ARRAY(7,/DOUBLE)  ;gamma 1/2
  EV1=MAKE_ARRAY(7,/DOUBLE) ;EigenValue
  EV2=EV1 & EV3=EV1 & v1=EV1 & v2=EV1 & v3=EV1  ;1 for j-1/2, 2 for j+1/2, 3 for j+3/2
  
  FU1_2=0.5*(U2F(REFORM(U[j,*]),gm)+U2F(REFORM(U[j+1,*]),gm))  ;1/2*(F(U_j)+F(U_j+1)
  ;for j+1/2
  v2=Evelocity(U2W(0.5*(REFORM(U[j,*])+REFORM(U[j+1,*])),gm),gm)
  EV2[0]=v2[4]-v2[3] & EV2[1]=v2[4]-v2[1] & EV2[2]=v2[4]-v2[2] & EV2[3]=v2[4] & EV2[6]=v2[4]+v2[3] & EV2[5]=v2[4]+v2[1] & EV2[4]=v2[4]+v2[2]
  nu2=c*EV2
  L2=LeftEV(U2W(0.5*(REFORM(U[j,*]+U[j+1,*])),gm),gm)
  R2=RightEV(U2W(0.5*(REFORM(U[j,*]+U[j+1,*])),gm),gm)
  alpha2=L2#REFORM(U[j+1,*]-U[j,*])
  gp2=0.5*(Q(nu2)-nu2^2)*alpha2
  s2=SIGNUM(gp2)
  
  IF j NE nx-2 THEN BEGIN
    ;for j+3/2
    v3=Evelocity(U2W(0.5*REFORM(U[j+2,*]+U[j+1,*]),gm),gm)
    EV3[0]=v3[4]-v3[3] & EV3[1]=v3[4]-v3[1] & EV3[2]=v3[4]-v3[2] & EV3[3]=v3[4] & EV3[6]=v3[4]+v3[3] & EV3[5]=v3[4]+v3[1] & EV3[4]=v3[4]+v3[2]
    nu3=3*EV2
    L3=LeftEV(U2W(0.5*REFORM(U[j+2,*]+U[j+1,*]),gm),gm)
    ;R3=RightEV(U2W(0.5*(U[j+2,*]+U[j+1,*]),gm),gm)
    alpha3=L3#REFORM(U[j+2,*]-U[j+1,*])
    gp3=0.5*(Q(nu3)-nu3^2)*alpha3
    s3=SIGNUM(gp3)
  ENDIF
  IF j NE 0 THEN BEGIN
    ;for j-1/2
    v1=Evelocity(U2W(0.5*REFORM(U[j,*]+U[j-1,*]),gm),gm)
    EV1[0]=v1[4]-v1[3] & EV1[1]=v1[4]-v1[1] & EV1[2]=v1[4]-v1[2] & EV1[3]=v1[4] & EV1[6]=v1[4]+v1[3] & EV1[5]=v1[4]+v1[1] & EV1[4]=v1[4]+v1[2]
    nu1=c*EV1
    L1=LeftEV(U2W(0.5*REFORM(U[j,*]+U[j-1,*]),gm),gm)
    ;R1=RightEV(U2W(0.5*(U[j,*]+U[j-1,*]),gm),gm)
    alpha1=L1#REFORM(U[j,*]-U[j-1,*])
    gp1=0.5*(Q(nu1)-nu1^2)*alpha1
    s1=SIGNUM(gp1)
  ENDIF
  
  IF j EQ 0 THEN BEGIN
    g2=s3*MAX([0,MIN([ABS(gp3),gp2*s3])])
    g1=g2
  ENDIF ELSE BEGIN
    IF j EQ nx-2 THEN BEGIN
      g1=s2*MAX([0,MIN([ABS(gp2),gp1*s2])])
      g2=g1
    ENDIF ELSE BEGIN
      g1=s2*MAX([0,MIN([ABS(gp2),gp1*s2])])
      g2=s3*MAX([0,MIN([ABS(gp3),gp2*s3])])
    ENDELSE
  ENDELSE
  g1_2=gm1_2(alpha2,g2,g1)
  RETURN,FU1_2+0.5/c*(R2#[g1+g2-Q(nu2+g1_2)*alpha2])
END