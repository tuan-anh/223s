     subroutine physics(x,y,z,ij)
      real par(14)
      common/model/profmod(2500,38),garea(1000)
      common/constant/pi,sigmav,rmin,rmax,xdelta
      data par/12.,.01,0.,0.1,0.,0.,-.5,.2,.01,0.,0.3,0.,0.,0./
c parameters which cannot change much
      theta=3.*pi/180.
      omega=13.*pi/180.
      cjet=0.93
      rdisk=1.5
c jet expansion velocity (~10.)
      vjet=par(1)
c jet thickness (~.05)
      dcjet=par(2)
c jet velocity gradient at small r (~.2)
      pjet=par(3)
c jet density enhancement (~.5)
      rhojet=par(4)
c slow wind expansion velocity (~1.)
      vwind=par(5)
c disk rotation velocity (~1.,-.5)
      vrot=par(6)
      prot=par(7)
c disk thickness (~.5)
      adisk=par(8)
c disk limb (~.1)
      dfr=par(9)
c disk density (~1.)
      rhodisk=par(10)
c smearing
      smr0=par(11)
      smrdisk=par(12)
      smrjet=par(13)
      rho=0.	
      v=0.
      r2=x**2+y**2+z**2
      r=sqrt(r2)
      if(r.lt.rmin.or.r.gt.rmax)return
      ct=cos(theta)
      st=sin(theta)
      co=cos(omega)
      so=sin(omega)
c one gets cphi from dot product of OM and PPN
      cphi=(x/r)*st+(y/r)*ct*so+(z/r)*ct*co 
c calpha for projection of radial velocity on line of sight
      calpha=x/r
c one gets tangential velocity from cross product of OM and PPN
c cvel for projection of tangential velocity on line of sight
      velx=(y/r)*ct*co-(z/r)*ct*so
      vely=(z/r)*st-(x/r)*ct*co
      velz=(x/r)*ct*so-(y/r)*st
      vel=sqrt(velx**2+vely**2+velz**2)
      cvel=velx/vel
c
      shape=exp(-.5*(r*cphi)**2/adisk**2)
      fc=(abs(cphi)-cjet)/cjet
      if(fc.gt.1.) fc=1.
      fr=r/rdisk
      fdisk=shape/(1.+exp(-(fr-1.)/dfr))
      gr=1.
      if(fr.lt.1.)gr=fr**pjet
      fjet=exp(-.5*fc**2/dcjet**2)*gr      
      fjet=0.
      if(abs(cphi).gt..85.and.abs(cphi).lt..86)fjet=1.
      vexp=vwind+vjet*fjet
      dv=vrot*fr**prot
      rho=(rhojet*fjet)+rhodisk*fdisk
      v=vexp*calpha+dv*cvel
      sigmav=smr0+smrdisk*fdisk+smrjet*fjet
      fabs=0.
      call fillprofmod(v,rho,ij)  
      return
      end
