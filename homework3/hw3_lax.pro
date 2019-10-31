pro hw3_lax

;setting parameters:dt/dx=0.18 and 300,400 or 600 points
cor=0.15
nx=300
nt=long(0.14*nx/cor)
w=dblarr(3,2*nx,nt+1)
gam=1.4
x=(indgen(2*nx)-(nx-1))*1.0/nx
dx=x[1]-x[0]

;setting initial values
for i=0,nx-2 do begin
    w[*,i,0]=[0.445,0.311,8.928]
endfor
for i=nx,nx*2-1 do begin
    w[*,i,0]=[0.5,0,1.4275]
endfor        

;computing lax form
w[*,nx-1,0]=(w[*,nx-2,0]+w[*,nx,0])*0.5    
for i=0,nt-1 do begin
    for j=1,2*nx-2 do begin
        aplus=mata(0.5*(w[*,j,i]+w[*,j+1,i]),gam) 
        aminus=mata(0.5*(w[*,j,i]+w[*,j-1,i]),gam)
        f=fw(w[*,j,i],gam)
        fplus=fw(w[*,j+1,i],gam)
        fminus=fw(w[*,j-1,i],gam)
        w[*,j,i+1]=w[*,j,i]-0.5*cor*(fplus-fminus)+0.5*cor^2*(aplus#(fplus-f)-aminus#(f-fminus))
    endfor
    w[*,2*nx-1,i+1]=w[*,2*nx-2,i+1]
    w[*,0,i+1]=w[*,1,i+1]
endfor
print,aplus
;print,w[0,*,50]

;setting parameters of analytic resolution

;plot rho,m,E
xref=[-0.9,-0.369,-0.229,0.214,0.214,0.347,0.347,0.6]
wref=[[0.445,0.445,0.345,0.345,1.304,1.304,0.500,0.500],$
      [0.311,0.311,0.527,0.527,1.994,1.994,0,0],$
      [8.928,8.928,6.570,6.570,7.691,7.691,1.428,1.428]]
strtt=['density','mass flux','energy']
set_plot,'ps'
loadct,39
device,filename='hw3_lax_600'+'.eps',/color,ENCAPSULATED=1
!p.multi=[0,0,3]
for k=0,2 do begin
    plot,x,w[k,*,nt],ytitle=strtt[k],xtitle='x',yrange=[1.1*min(w[k,*,nt])-0.1*max(w[k,*,nt]),1.1*max(w[k,*,nt])-0.1*min(w[k,*,nt])],charsize=1.6,charthick=2,psym=4,symsize=0.85
    oplot,x,w[k,*,nt],linestyle=4
    oplot,xref,wref[*,k]
    oplot,x,w[k,*,0],linestyle=1
endfor
xyouts,-0.5,55,'Lax-Windroff scheme with 600 points',charthick=2
device,/close
;print,w[0,*,nt]
end


function mata,wp,gamp

;computing matrix A
a=dblarr(3,3)
u=wp[1]/wp[0]
a[0,1]=1
a[1,0]=0.5*(gamp-3)*u^2
a[1,1]=-(gamp-3)*u
a[1,2]=gamp-1
a[2,0]=(gamp-1)*u^3-gamp*u/wp[0]*wp[2]
a[2,1]=gamp*wp[2]/wp[0]-1.5*(gamp-1)*u^2
a[2,2]=gamp*u

return,a

end

function fw,wp1,gamp1

;computing vector f(w)
a=dblarr(3)
a[0]=wp1[1]
a[1]=(gamp1-1)*wp1[2]+0.5*(3-gamp1)*wp1[1]^2/wp1[0]
a[2]=(gamp1*wp1[2]-0.5*(gamp1-1)*wp1[1]^2/wp1[0])*wp1[1]/wp1[0]
return,a
end