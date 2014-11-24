  subroutine physics(x,y,z,ij)
      real par(13)
      common/model/profmod(2500,38),garea(1000)
      common/constant/pi,sigmav,rmin,rmax,xdelta
      data par/5.,14.,.87,1.,6.,1.35,1.5,3.,1.7,2.7,1.5,0.,-.2/
      theta=par(1)*pi/180.
      omega=par(2)*pi/180.
      cbicone=par(3)
      vdisk=par(4)
      vbicone=par(5)
      vrot=par(6)
      rrot=par(7)
      adisk=par(8)
      delta=par(9)
      power=par(10)
      sig0=par(11)
      sig1=par(12)
      rho=0.
      v=0.
      dv=0.
      r2=x**2+y**2+z**2
      r=sqrt(r2)
      if(r.lt.rmin.or.r.gt.rmax)return
      sigmav=sig0+sig1*(rrot-r)/rrot
      ct=cos(theta)
      st=sin(theta)
      co=cos(omega)
      so=sin(omega)
      cphi=(x/r)*st+(y/r)*ct*so+(z/r)*ct*co 
       if(abs(cphi).gt.cbicone) return
c calpha for projection of radial velocity on line of sight
      calpha=x/r
c cvel for projection of tangential velocity on line of sight
      velx=(y/r)*ct*co-(z/r)*ct*so
      vely=(z/r)*st-(x/r)*ct*co
      velz=(x/r)*ct*so-(y/r)*st
      vel=sqrt(velx**2+vely**2+velz**2)
      cvel=velx/vel
c      power1=alog(vbicone/vdisk)/alog(rrot-rmin)
      v0=(abs(cphi)/cbicone)**power*(vbicone-vdisk)+vdisk
      v=v0*calpha 
      if(r.ge.rrot)rho=1./4./pi*(rrot/r)**2
      if(r.lt.rrot)then
         shape=exp(-.5*(r*cphi)**2/delta**2)
         rho=1/8./pi*(3.-(r/rrot)**2)*(1.+adisk*shape)
         dv=vrot*shape*(r/rrot)**(-.5)
         v=v+dv*cvel
      endif           
      rho=rho*(1.+par(13)*r/rrot)
      fabs=0.
      call fillprofmod(v,rho,ij)  
      return
      end
