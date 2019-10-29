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

FUNCTION V_aver,u,v
;; a smooth function V(u,v) for average of u and v

  RETURN,0.5*(u+v)
END

FUNCTION L_fun,u,a,gam
;; Left eigenvectors computation

L=FLTARR(3,3)
L[*,*]=0
H=a^2/(gam-1)+0.5*u^2
L[0,0]=0.5*u*(u+2*a/(gam-1))
L[1,0]=-(u+a/(gam-1))
L[2,0]=1
L[0,1]=2*(H-u^2)
L[1,1]=2*u
L[2,1]=-2
L[0,2]=0.5*u*(u-2*a/(gam-1))
L[1,2]=-(u-a/(gam-1))
L[2,2]=1
L=L*(gam-1)/(2*a^2)
RETURN,L
END

FUNCTION R_fun,u,a,gam
;; Right eigenvectors computation

R=FLTARR(3,3)
H=a^2/(gam-1)+0.5*u^2
R[0,0]=1
R[0,1]=u-a
R[0,2]=H-u*a
R[1,0]=1
R[1,1]=u
R[1,2]=0.5*u^2
R[2,0]=1
R[2,1]=u+a
R[2,2]=H+u*a
RETURN,R
END

FUNCTION TVD,rho,m,E,u,p,a,n,j,gam,c,num,no
;;For computation of TVD scheme

;d=0 for j-1/2 if no other notes 
;d=1 for j+1/2 if no other notes 
;d=2 for j+3/2 if no other notes
a0=FLTARR(3,3)
v=FLTARR(3,3)
L=FLTARR(3,3,3)
R=FLTARR(3,3,3)
alph=FLTARR(3,3)
g_p=FLTARR(3,3)
s=FLTARR(3,3)
g=FLTARR(3,3)
y=FLTARR(3)
ep=0.1

;------------------density------------------
IF no EQ 0 THEN BEGIN
  f0=m[j,n] ;f_j
  f1=m[j+1,n]  ;f_j+1
  IF (j+2 LT num) AND (j GT 0) THEN BEGIN
    FOR d=0,2 DO BEGIN
       a0[0:2,d]=V_aver(u[j-1+d],u[j+d])+([0:2]-1)*V_aver(a[j-1+d],a[j+d])
       v[*,d]=c*a0[*,d]
       L[*,*,d]=V_aver(L_fun(u[j-1+d],a[j-1+d],gam),L_fun(u[j+d],a[j+d],gam))
       alph[*,d]=L[0,*,d]*(rho[j+d,n]-rho[j-1+d,n])+$
         L[1,*,d]*(m[j+d,n]-m[j-1+d,n])+$
         L[2,*,d]*(E[j+d,n]-E[j-1+d,n])
       R[*,*,d]=V_aver(R_fun(u[j-1+d],a[j-1+d],gam),R_fun(u[j+d],a[j+d],gam))
       g_p[*,d]=0.5*(Q(v[*,d],ep)-v[*,d]^2)*alph[*,d]
    ENDFOR
    s=SIGNUM(g_p)
    g0=s[*,1]*MAX([0,MIN([ABS(g_p[*,1]),g_p[*,0]*s[*,1]])]) ;g_j
    g1=s[*,2]*MAX([0,MIN([ABS(g_p[*,2]),g_p[*,1]*s[*,2]])]) ;g_j+1
  ENDIF
  IF j+2 EQ num THEN BEGIN
    FOR d=0,1 DO BEGIN
      a0[0:2,d]=V_aver(u[j-1+d],u[j+d])+([0:2]-1)*V_aver(a[j-1+d],a[j+d])
       v[*,d]=c*a0[*,d]
       L[*,*,d]=V_aver(L_fun(u[j-1+d],a[j-1+d],gam),L_fun(u[j+d],a[j+d],gam))
       alph[*,d]=L[0,*,d]*(rho[j+d,n]-rho[j-1+d,n])+$
         L[1,*,d]*(m[j+d,n]-m[j-1+d,n])+$
         L[2,*,d]*(E[j+d,n]-E[j-1+d,n])
       R[*,*,d]=V_aver(R_fun(u[j-1+d],a[j-1+d],gam),R_fun(u[j+d],a[j+d],gam))
       g_p[*,d]=0.5*(Q(v[*,d],ep)-v[*,d]^2)*alph[*,d]
    ENDFOR
    s=SIGNUM(g_p)
    g0=s[*,1]*MAX([0,MIN([ABS(g_p[*,1]),g_p[*,0]*s[*,1]])]) ;g_j
    g1=g0 ;g_j+1
  ENDIF
  IF j EQ 0 THEN BEGIN
    FOR d=1,2 DO BEGIN
      a0[0:2,d]=V_aver(u[j-1+d],u[j+d])+([0:2]-1)*V_aver(a[j-1+d],a[j+d])
       v[*,d]=c*a0[*,d]
       L[*,*,d]=V_aver(L_fun(u[j-1+d],a[j-1+d],gam),L_fun(u[j+d],a[j+d],gam))
       alph[*,d]=L[0,*,d]*(rho[j+d,n]-rho[j-1+d,n])+$
         L[1,*,d]*(m[j+d,n]-m[j-1+d,n])+$
         L[2,*,d]*(E[j+d,n]-E[j-1+d,n])
       R[*,*,d]=V_aver(R_fun(u[j-1+d],a[j-1+d],gam),R_fun(u[j+d],a[j+d],gam))
       g_p[*,d]=0.5*(Q(v[*,d],ep)-v[*,d]^2)*alph[*,d]
    ENDFOR
    s=SIGNUM(g_p)
    g1=s[*,2]*MAX([0,MIN([ABS(g_p[*,2]),g_p[*,1]*s[*,2]])]) ;g_j+1
    g0=g1 ;g_j
  ENDIF

  FOR k=0,2 DO BEGIN
    IF alph[k,1] NE 0 THEN BEGIN
      y[k]=(g1[k]-g0[k])/alph[k,1]
    ENDIF ELSE BEGIN
      y[k]=0
    ENDELSE
  ENDFOR
  Q1=Q(v[*,1]+y,ep)
  RETURN,0.5*(f0+f1)+1/(2*c)*TOTAL((g1+g0-Q1*alph[*,1])*R[*,no,1])
ENDIF

;--------------mass------------------------
IF no EQ 1 THEN BEGIN
  f0=(gam-1)*E[j,n]+(3-gam)*m[j,n]^2/2./rho[j,n] ;f_j
  f1=(gam-1)*E[j+1,n]+(3-gam)*m[j+1,n]^2/2./rho[j+1,n]  ;f_j+1
  IF (j+2 LT num) AND (j GT 0) THEN BEGIN
    FOR d=0,2 DO BEGIN
      a0[0:2,d]=V_aver(u[j-1+d],u[j+d])+([0:2]-1)*V_aver(a[j-1+d],a[j+d])
       v[*,d]=c*a0[*,d]
       L[*,*,d]=V_aver(L_fun(u[j-1+d],a[j-1+d],gam),L_fun(u[j+d],a[j+d],gam))
       alph[*,d]=L[0,*,d]*(rho[j+d,n]-rho[j-1+d,n])+$
         L[1,*,d]*(m[j+d,n]-m[j-1+d,n])+$
         L[2,*,d]*(E[j+d,n]-E[j-1+d,n])
       R[*,*,d]=V_aver(R_fun(u[j-1+d],a[j-1+d],gam),R_fun(u[j+d],a[j+d],gam))
       g_p[*,d]=0.5*(Q(v[*,d],ep)-v[*,d]^2)*alph[*,d]
    ENDFOR
    s=SIGNUM(g_p)
    g0=s[*,1]*MAX([0,MIN([ABS(g_p[*,1]),g_p[*,0]*s[*,1]])]) ;g_j
    g1=s[*,2]*MAX([0,MIN([ABS(g_p[*,2]),g_p[*,1]*s[*,2]])]) ;g_j+1
  ENDIF
  IF j+2 EQ num THEN BEGIN
    FOR d=0,1 DO BEGIN
      a0[0:2,d]=V_aver(u[j-1+d],u[j+d])+([0:2]-1)*V_aver(a[j-1+d],a[j+d])
       v[*,d]=c*a0[*,d]
       L[*,*,d]=V_aver(L_fun(u[j-1+d],a[j-1+d],gam),L_fun(u[j+d],a[j+d],gam))
       alph[*,d]=L[0,*,d]*(rho[j+d,n]-rho[j-1+d,n])+$
         L[1,*,d]*(m[j+d,n]-m[j-1+d,n])+$
         L[2,*,d]*(E[j+d,n]-E[j-1+d,n])
       R[*,*,d]=V_aver(R_fun(u[j-1+d],a[j-1+d],gam),R_fun(u[j+d],a[j+d],gam))
       g_p[*,d]=0.5*(Q(v[*,d],ep)-v[*,d]^2)*alph[*,d]
    ENDFOR
    s=SIGNUM(g_p)
    g0=s[*,1]*MAX([0,MIN([ABS(g_p[*,1]),g_p[*,0]*s[*,1]])]) ;g_j
    g1=g0 ;g_j+1
  ENDIF
  IF j EQ 0 THEN BEGIN
    FOR d=1,2 DO BEGIN
      a0[0:2,d]=V_aver(u[j-1+d],u[j+d])+([0:2]-1)*V_aver(a[j-1+d],a[j+d])
       v[*,d]=c*a0[*,d]
       L[*,*,d]=V_aver(L_fun(u[j-1+d],a[j-1+d],gam),L_fun(u[j+d],a[j+d],gam))
       alph[*,d]=L[0,*,d]*(rho[j+d,n]-rho[j-1+d,n])+$
         L[1,*,d]*(m[j+d,n]-m[j-1+d,n])+$
         L[2,*,d]*(E[j+d,n]-E[j-1+d,n])
       R[*,*,d]=V_aver(R_fun(u[j-1+d],a[j-1+d],gam),R_fun(u[j+d],a[j+d],gam))
       g_p[*,d]=0.5*(Q(v[*,d],ep)-v[*,d]^2)*alph[*,d]
    ENDFOR
    s=SIGNUM(g_p)
    g1=s[*,2]*MAX([0,MIN([ABS(g_p[*,2]),g_p[*,1]*s[*,2]])]) ;g_j+1
    g0=g1 ;g_j
  ENDIF
  FOR k=0,2 DO BEGIN
    IF alph[k,1] NE 0 THEN BEGIN
      y[k]=(g1[k]-g0[k])/alph[k,1]
    ENDIF ELSE BEGIN
      y[k]=0
    ENDELSE
  ENDFOR
  Q1=Q(v[*,1]+y,ep)
  RETURN,0.5*(f0+f1)+1/(2*c)*TOTAL((g1+g0-Q1*alph[*,1])*R[*,no,1])
ENDIF

;--------------energy-----------------------
IF no EQ 2 THEN BEGIN
  f0=(gam*E[j,n]-(gam-1)/2*m[j,n]^2/rho[j,n])*m[j,n]/rho[j,n] ;f_j
  f1=(gam*E[j+1,n]-(gam-1)/2*m[j+1,n]^2/rho[j+1,n])*m[j+1,n]/rho[j+1,n]  ;f_j+1
  IF (j+2 LT num) AND (j GT 0) THEN BEGIN
    FOR d=0,2 DO BEGIN
      a0[0:2,d]=V_aver(u[j-1+d],u[j+d])+([0:2]-1)*V_aver(a[j-1+d],a[j+d])
       v[*,d]=c*a0[*,d]
       L[*,*,d]=V_aver(L_fun(u[j-1+d],a[j-1+d],gam),L_fun(u[j+d],a[j+d],gam))
       alph[*,d]=L[0,*,d]*(rho[j+d,n]-rho[j-1+d,n])+$
         L[1,*,d]*(m[j+d,n]-m[j-1+d,n])+$
         L[2,*,d]*(E[j+d,n]-E[j-1+d,n])
       R[*,*,d]=V_aver(R_fun(u[j-1+d],a[j-1+d],gam),R_fun(u[j+d],a[j+d],gam))
       g_p[*,d]=0.5*(Q(v[*,d],ep)-v[*,d]^2)*alph[*,d]
    ENDFOR
    s=SIGNUM(g_p)
    g0=s[*,1]*MAX([0,MIN([ABS(g_p[*,1]),g_p[*,0]*s[*,1]])]) ;g_j
    g1=s[*,2]*MAX([0,MIN([ABS(g_p[*,2]),g_p[*,1]*s[*,2]])]) ;g_j+1
  ENDIF
  IF j+2 EQ num THEN BEGIN
    FOR d=0,1 DO BEGIN
      a0[0:2,d]=V_aver(u[j-1+d],u[j+d])+([0:2]-1)*V_aver(a[j-1+d],a[j+d])
       v[*,d]=c*a0[*,d]
       L[*,*,d]=V_aver(L_fun(u[j-1+d],a[j-1+d],gam),L_fun(u[j+d],a[j+d],gam))
       alph[*,d]=L[0,*,d]*(rho[j+d,n]-rho[j-1+d,n])+$
         L[1,*,d]*(m[j+d,n]-m[j-1+d,n])+$
         L[2,*,d]*(E[j+d,n]-E[j-1+d,n])
       R[*,*,d]=V_aver(R_fun(u[j-1+d],a[j-1+d],gam),R_fun(u[j+d],a[j+d],gam))
       g_p[*,d]=0.5*(Q(v[*,d],ep)-v[*,d]^2)*alph[*,d]
    ENDFOR
    s=SIGNUM(g_p)
    g0=s[*,1]*MAX([0,MIN([ABS(g_p[*,1]),g_p[*,0]*s[*,1]])]) ;g_j
    g1=g0 ;g_j+1
  ENDIF
  IF j EQ 0 THEN BEGIN
    FOR d=1,2 DO BEGIN
      a0[0:2,d]=V_aver(u[j-1+d],u[j+d])+([0:2]-1)*V_aver(a[j-1+d],a[j+d])
       v[*,d]=c*a0[*,d]
       L[*,*,d]=V_aver(L_fun(u[j-1+d],a[j-1+d],gam),L_fun(u[j+d],a[j+d],gam))
       alph[*,d]=L[0,*,d]*(rho[j+d,n]-rho[j-1+d,n])+$
         L[1,*,d]*(m[j+d,n]-m[j-1+d,n])+$
         L[2,*,d]*(E[j+d,n]-E[j-1+d,n])
       R[*,*,d]=V_aver(R_fun(u[j-1+d],a[j-1+d],gam),R_fun(u[j+d],a[j+d],gam))
       g_p[*,d]=0.5*(Q(v[*,d],ep)-v[*,d]^2)*alph[*,d]
    ENDFOR
    s=SIGNUM(g_p)
    g1=s[*,2]*MAX([0,MIN([ABS(g_p[*,2]),g_p[*,1]*s[*,2]])]) ;g_j+1
    g0=g1 ;g_j
  ENDIF
  FOR k=0,2 DO BEGIN
    IF alph[k,1] NE 0 THEN BEGIN
      y[k]=(g1[k]-g0[k])/alph[k,1]
    ENDIF ELSE BEGIN
      y[k]=0
    ENDELSE
  ENDFOR
  Q1=Q(v[*,1]+y,ep)
  RETURN,0.5*(f0+f1)+1./(2*c)*TOTAL((g1+g0-Q1*alph[*,1])*R[*,no,1])
ENDIF
END