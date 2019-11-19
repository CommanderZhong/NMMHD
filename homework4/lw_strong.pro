;;assign4  strong shock fast & slow
;;
;;2019-11-17, Hanya Pan

common share1,bt,gm,Hx,nx,dx,dt

function inicon,xs,fs=fs,n=n,x0=x0
  if keyword_set(fs) then fs=fs else fs = 'fast'
  if keyword_set(n) then n=n else n = 1
  if keyword_set(x0) then x0=x0 else x0 = 0.
  if n eq 1 then begin
    if fs eq 'fast' then begin
      WL=[2.121,4.981,-13.27,-0.163,-0.6521,2.572,10.29]  ;fast shock
      WR=[1.,1.,-15.3,0.,0.,1.,4.]
    endif else begin
      WL=[2.219,0.4442,0.5048,0.0961,0.0961,1.,1.]  ;slow shock
      WR=[1.,0.1,-0.9225,0.,0.,1.,1.]
    endelse
  endif else begin
    if fs eq 'fast' then begin
      WL=[3.896,305.9,0.,-0.058,-0.226,3.951,15.8]  ;fast shock
      WR=[1.,1.,-15.3,0.,0.,1.,4.]
    endif else begin
      WL=[3.108,1.4336,0.,0.2633,0.2633,0.1,0.1]  ;slow shock
      WR=[1.,0.1,-0.9225,0.,0.,1.,1.]
    endelse
  endelse
  deunit=[1.,1.,1.,1.,1.,sqrt(4.*!pi),sqrt(4.*!pi)]
  wl=wl/deunit
  wr=wr/deunit

  indl=where(xs lt x0)
  indr=where(xs gt x0)
  w0=make_array(n_elements(xs),7,/float)
  w0[indl,0]=wl[0]
  w0[indl,1]=wl[1]
  w0[indl,2]=wl[2]
  w0[indl,3]=wl[3]
  w0[indl,4]=wl[4]
  w0[indl,5]=wl[5]
  w0[indl,6]=wl[6]

  w0[indr,0]=wr[0]
  w0[indr,1]=wr[1]
  w0[indr,2]=wr[2]
  w0[indr,3]=wr[3]
  w0[indr,4]=wr[4]
  w0[indr,5]=wr[5]
  w0[indr,6]=wr[6]
return,w0
end


function w2u,w ;[4,7],[lf rf ls rs]*7
common share1
  u=w
  u[*,0]=w[*,0]
  u[*,1]=w[*,0]*(w[*,2]^2+w[*,3]^2+w[*,4]^2)+w[*,5]^2+w[*,6]^2+bt*w[*,1]/(gm-1)
  u[*,2]=w[*,0]*w[*,2]
  u[*,3]=w[*,0]*w[*,3]
  u[*,4]=w[*,0]*w[*,4]
  u[*,5]=w[*,5]
  u[*,6]=w[*,6]
return,u
end


function u2f, u
common share1

  vx=u[*,2]/u[*,0] & vy=u[*,3]/u[*,0] & vz=u[*,4]/u[*,0]
  v2=(u[*,2]^2+u[*,3]^2+u[*,4]^2)/u[*,0]^2  ;v^2
  rho=u[*,0]
  Hy=u[*,5]
  Hz=u[*,6]
  p=(u[*,1]-rho*v2-Hy^2-Hz^2)*(gm-1)/bt

  f=u
  f[*,0]=rho*vx
  f[*,1]=rho*vx*(v2+gm/(gm-1)*bt*p/rho)+2*(Hy^2*vx+Hz^2*vx-Hx*Hy*vy-Hx*Hz*vz)
  f[*,2]=rho*vx^2+0.5*bt*p+0.5*(Hy^2+Hz^2)
  f[*,3]=rho*vx*vy-Hx*Hy
  f[*,4]=rho*vx*vz-Hx*Hz
  f[*,5]=vx*Hy-vy*Hx
  f[*,6]=vx*Hz-vz*Hx

;CFL----------
vel=sqrt(gm*p/u[*,0])+abs(vx)
CFL=real_part(dx/max(vel))
dt=0.5*CFL

  return, f
end


function lw,u0
common share1
f=u2f(u0)  ;;get initial value of F(U) [nx,7]

u1=u0
for j=1,nx-2 do begin
  u1[j,*]=-0.5*dt/dx*(f[j+1,*]-f[j-1,*])+0.5*(u0[j+1,*]+u0[j-1,*])
endfor
return,u1
end

c=0.005  ;dt/dx
nx=1000
tc=[0.05,0.1] & t=0.

bt=2.
gm=5./3.
Hx=5./sqrt(4.*!pi)
xst=0. & xrange=1.
xs=findgen(nx)/(nx-1)*xrange+xst
dx=xrange/(nx-1)
w0=inicon(xs,fs='fast',n=2,x0=0.2)  ;;get initial value of W of slow or fast shock
u0=w2u(w0)  ;;get initial value of U [nx,7]
u00=u0

while t lt tc[1] do begin
  u1=lw(u0)
  u0=u1
  t=t+dt
print,t
  if t lt tc[0] then u2=u1
endwhile

H2=Hx^2+u1[*,5]^2+u1[*,6]^2
vx=u1[*,2]/u1[*,0]
vy=u1[*,3]/u1[*,0]
vz=u1[*,4]/u1[*,0]
v2=vx^2+vy^2+vz^2
p=(u1[*,1]-u1[*,0]*v2-u1[*,5]^2-u1[*,6]^2)*(gm-1)/bt

w1=u1
w1[*,0]=u1[*,0]
w1[*,1]=p
w1[*,2]=vx
w1[*,3]=vy
w1[*,4]=vz
w1[*,5]=u1[*,5]
w1[*,6]=u1[*,6]

H2=Hx^2+u2[*,5]^2+u2[*,6]^2
vx=u2[*,2]/u2[*,0]
vy=u2[*,3]/u2[*,0]
vz=u2[*,4]/u2[*,0]
v2=vx^2+vy^2+vz^2
p=(u2[*,1]-u2[*,0]*v2-u2[*,5]^2-u2[*,6]^2)*(gm-1)/bt

w2=u2
w2[*,0]=u2[*,0]
w2[*,1]=p
w2[*,2]=vx
w2[*,3]=vy
w2[*,4]=vz
w2[*,5]=u2[*,5]
w2[*,6]=u2[*,6]

u0=u00

H2=Hx^2+u0[*,5]^2+u0[*,6]^2
vx=u0[*,2]/u0[*,0]
vy=u0[*,3]/u0[*,0]
vz=u0[*,4]/u0[*,0]
v2=vx^2+vy^2+vz^2
p=(u0[*,1]-u0[*,0]*v2-u0[*,5]^2-u0[*,6]^2)*(gm-1)/bt

w0=u0
w0[*,0]=u0[*,0]
w0[*,1]=p
w0[*,2]=vx
w0[*,3]=vy
w0[*,4]=vz
w0[*,5]=u0[*,5]
w0[*,6]=u0[*,6]

E0=0.5*w0[*,0]*(w0[*,2]^2+w0[*,3]^2+w0[*,4]^2)+w0[*,1]/(gm-1)
E1=0.5*w1[*,0]*(w1[*,2]^2+w1[*,3]^2+w1[*,4]^2)+w1[*,1]/(gm-1)
E2=0.5*w2[*,0]*(w2[*,2]^2+w2[*,3]^2+w2[*,4]^2)+w2[*,1]/(gm-1)



xsize=10. & ysize=6.				;;size of the figure;
;; --- position info
;; !--------------------------------------------
ncol=2						;;assume there's ncol*nrow panels
nrow=4
ncplus=0					;;consider different layout types.
nrplus=0					;;  ex:goes curve and several aia obs
width=0.82					;;width & length of the whole panels
length=0.88
left_margin=0.07					;;move the panels' position
down_margin=0.07
x_margin=0.08
y_margin=0.01
;; ---------------------------------------------

npanel=ncol*nrow				;;number of the panels
pos=make_array(4,npanel,/float)			;;create and calculate the positions
for i=0,npanel-1 do begin
	x0=(i mod ncol)*width/ncol+left_margin+(i mod ncol)*x_margin
	y0=(nrow-i/ncol-1)*length/nrow+down_margin+(nrow-i/ncol-1)*y_margin
	x1=(i mod ncol +1)*width/ncol+left_margin+(i mod ncol)*x_margin
	y1=(nrow-i/ncol)*length/nrow+down_margin+(nrow-i/ncol-1)*y_margin
	pos[*,i]=[x0,y0,x1,y1]			;;usage: !p.position=pos[*,i]
endfor

set_plot,'ps'
!p.font=-1
!p.MULTI=[0,ncol+ncplus,nrow+nrplus]		;;layout
dir='/home/hannah/Documents/2019fallsimulation/as4_m/'	;;directory
fn=str_replace(dir+'lw_strong_fast'+'.eps',' ','')	;;filename
device,filename=fn,xs=xsize,ys=ysize,/color,BITS_PER_PIXEL=24,set_font='Helvetica'
!p.charsize=0.7					;;size of characters

loadct,0
tvlct, 0,0,0,0
tvlct, 180,0,180,2
tvlct, 180,0,0,1

!p.position=pos[*,0]
plot,xs,w0[*,0],linest=1,ytitle='!4q!x',xtickformat='(A1)',yrange=[0,5]
oplot,xs,w1[*,0],color=1
oplot,xs,w2[*,0],color=2

!p.position=pos[*,1]
plot,xs,w0[*,1],linest=1,ytitle='p',xtickformat='(A1)',yrange=[-100,400]
oplot,xs,w1[*,1],color=1
oplot,xs,w2[*,1],color=2

!p.position=pos[*,2]
plot,xs,w0[*,2],linest=1,ytitle='v!dx!n',xtickformat='(A1)',yrange=[-20,5]
oplot,xs,w1[*,2],color=1
oplot,xs,w2[*,2],color=2

!p.position=pos[*,3]
plot,xs,w0[*,3],linest=1,ytitle='v!dy!n',xtickformat='(A1)',yrange=[-0.1,0.02]
oplot,xs,w1[*,3],color=1
oplot,xs,w2[*,3],color=2

!p.position=pos[*,4]
plot,xs,w0[*,4],linest=1,ytitle='v!dz!n',xtickformat='(A1)',yrange=[-0.3,0.05]
oplot,xs,w1[*,4],color=1
oplot,xs,w2[*,4],color=2

!p.position=pos[*,5]
plot,xs,w0[*,5],linest=1,ytitle='H!dy!n',xtickformat='(A1)',yrange=[0,1.5]
oplot,xs,w1[*,5],color=1
oplot,xs,w2[*,5],color=2

!p.position=pos[*,6]
plot,xs,w0[*,6],linest=1,ytitle='H!dz!n',xtitle='x',yrange=[0,5]
oplot,xs,w1[*,6],color=1
oplot,xs,w2[*,6],color=2

!p.position=pos[*,7]
plot,xs,E0,linest=1,ytitle='E',xtitle='x',yrange=[0,500]
oplot,xs,E1,color=1
oplot,xs,E2,color=2


device,/close
set_plot,'x'


end
