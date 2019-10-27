FUNCTION Q,x,ep
;;Equation (9.62) of the lecture note

  IF ABS(x) LT 2*ep THEN BEGIN
    Q=x^2/4/ep+ep
  ENDIF ELSE BEGIN
    Q=ABS(x)
  ENDELSE
  RETURN,Q
END

FUNCTION V_aver,u,v
;; a smooth function V(u,v) for average of u and v
  RETURN,2./(1./u+1./v)
END

PRO TVD,v,t0,gam
;;For computation of TVD scheme



END