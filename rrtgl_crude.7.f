cc
c This version for fiting outer part (cone only, no disk)
c 
c
      program rrtgl_crude  
c
c this version is for the Red rectangle using ALMA data
c
c spectra include 38 velocity bins, 0.4 km/s wide, from -7.2 to 8 km/s
c sky maps are centered on the star in 50*50 bins of 0.1*0.1 arcsec**2
c from -2.5 to 2.5 arcsec; first index i increases from east to west,
c but y increases from west to east, we need to change its sign. 
c second index j increases from south to north; ij=i+50*(j-1), 1 to 2500.
c ico is 1 for CO(3-2) and 2 for CO(6-5)
c
      common/model/profmod(2500,38),garea(1000)
      common/constant/pi,sigmav,rmin,rmax,xdelta
      common/meas/profmeas(2,2500,38),smeas(2,2500)  
      data xdelta,sigmav/0.02,1.5/
      data rmin,rmax/.1,5./
      call HLIMAP(20000,'MAP')
c      call hbook1(101,'cphi', 2500.,0.,2500.,0.)
      pi=2.*acos(0.)


c
cccc calculate the array garea(1000) for Gaussian smearing in fillprofmod
c 
      x=-2.997
      garea(1)=.006*exp(-0.5*x**2)
      do ig=2,1000
         x=x+0.006
         garea(ig)=garea(ig-1)+.006*exp(-0.5*x**2)
      enddo
      do ig=1,1000
         garea(ig)=garea(ig)/garea(1000)
      enddo   
c
c all coordinates in arcseconds centered on the centre of the star
c rmin and rmax define the range in r over which we work
c profmod is the spectrum predicted by the model 
c xdelta is the scan step in x
c
       do 1 ij=1,2500
         do k=1,38
            profmod(ij,k)=0.
         enddo
         j=ij/50+1
         i=ij-50*(j-1)             
         y=-2.5+0.1*float(i-1)+0.05
         y=-y 
         z=-2.5+0.1*float(j-1)+0.05
c
c y is positive to the east and z to the north
c
         xmax2=rmax**2-z**2-y**2
         xmax=sqrt(xmax2)             
         nmax=2*ifix(xmax/xdelta)+1
         xmax=float(nmax-1)*xdelta/2.
c
c loop along the line of sight
c
         do 4 ix=1,nmax
            x=xmax-float(ix-1)*xdelta
            rchk=sqrt(x**2+y**2+z**2)
            if(rchk.lt.rmin) goto 4          
            call physics(x,y,z,ij)
 4       continue
    1 continue     
      call readdata
      call plot

      stop
      end 
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      
      subroutine physics(x,y,z,ij)
      real par(14)
      common/model/profmod(2500,38),garea(1000)
      common/constant/pi,sigmav,rmin,rmax,xdelta
      data par/13.,2.65,1.,.7,-2.,0.,-.5,.2,.01,0.,0.,0.,0.,0./
c parameters which cannot change much
      theta=1.*pi/180.
      omega=14.*pi/180.
      cjet=0.84
      rdisk=1.5
c jet expansion velocity (~10.)
      vjet=par(1)
c jet velocity gradient at small r (~.2)
      pjet=par(2)
c jet density enhancement (~.5)
      rhojet=par(3)
c slow wind expansion velocity (~1.)
      vwind=par(4)
c r power of density
      pexp=par(5)
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
      cphi=(x/r)*st+(y/r)*ct*so+(z/r)*ct*co 
      if(abs(cphi).gt.cjet) return
c calpha for projection of radial velocity on line of sight
      calpha=x/r
c cvel for projection of tangential velocity on line of sight
      velx=(y/r)*ct*co-(z/r)*ct*so
      vely=(z/r)*st-(x/r)*ct*co
      velz=(x/r)*ct*so-(y/r)*st
      vel=sqrt(velx**2+vely**2+velz**2)
      cvel=velx/vel
c
      shape=exp(-.5*(r*cphi)**2/adisk**2)
      fc=abs(cphi)/cjet
      fr=r/rdisk
      fdisk=shape/(1.+exp(-(fr-1.)/dfr))
      gr=1.
      if(fr.lt.1.)gr=fr**pjet
      vexp=(vwind+vjet*fc**pjet)
      dv=vrot*fr**prot
      rho=(1.+rhojet*fc**pjet)*fr**pexp
      v=vexp*calpha+dv*cvel
      sigmav=smr0+smrdisk*fdisk+smrjet*fc
      fabs=0.
      call fillprofmod(v,rho,ij)  
      return
      end
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      subroutine fillprofmod(v,rho,ij)
      common/model/profmod(2500,38),garea(1000)
      common/constant/pi,sigmav,rmin,rmax,xdelta
c	fill array profmod in 38 bins of 0.4 km/s between -7.2 and 8 km/s
c	smearing is over plus or minus 3 sigmas
c	deals with absorption
c	profmod(y,z) is reset after each x scan
c	profmod(y,z) is integrated along the line of sight, starting from far away and
c	moving toward the observer
c	absorption coefficient is 
c	kappa=(c**2*nl*gmu)*(8.*pi*freq**2*gl)**(-1)*Aul*(1.-gl*nu/gu/nl)
c	u and l stand for upper and lower levels; n for density of atoms; 
c	g for statistical weights 
c	the evolution of the outcoming flux is
c	I(x+dx)=I(x)+Iemitted*dx-kappa*I(x)*dx
c	with Iemitted=h*freq*(4.*pi)**(-1)*nu*Aul
c	with LTE Boltzman gives nu/nl=(gu/gl)exp(-h*freq/(kboltzman*temperature))
c	Hence I(x+dx)=I(x)+ h*freq*(4.*pi)**(-1)*nu*Aul*(1.-fabs*I(x))*dx
c 	with fabs=c**2(2.*h*freq**3)**(-1)*(exp(h*freq/(kb*T))-1.)
c	namely I(x+dx)=I(x)+rho*xdelta*(1.-fabs*I(x))
      vmin=v-3.*sigmav
      if(vmin.lt.-7.2) vmin=-7.2
      if(vmin.gt.8.) go to 2
      vmax=v+3.*sigmav
      if(vmax.gt.8.) vmax=8.
      if(vmax.lt.-7.2) go to 2
      if(vmin.gt.vmax) go to 2
      ivmin=ifix((vmin+7.2)/0.4)+1
      ivmax=ifix((vmax+7.2)/0.4)+1
      do iv=ivmin,ivmax
         vlo=-7.2+0.4*float(iv-1)
         vhi=vlo+0.4
         ihi=1
         ilo=1000
         if(sigmav.le.0.)then
            weight=1.
            goto 1
         endif
         ihi=ifix(((vhi-v)/sigmav+3.)/.006)
         if(ihi.gt.1000)ihi=1000
         ilo=ifix(((vlo-v)/sigmav+3.)/.006)
         if(ilo.lt.1)ilo=1
         weight=0.
         if(ihi.lt.1.or.ihi.gt.1000.or.ilo.lt.1.or.ilo.gt.1000)
     &   goto 1
         weight=garea(ihi)-garea(ilo)
 1       continue
         profmod(ij,iv)=profmod(ij,iv)
     &        +rho*xdelta*weight
      call hf1(101,float(ij)-.5,abs(vlo))
      enddo
 2    continue
      return
      end
c
c
c
c
      subroutine plot
      common/model/profmod(2500,38),garea(1000)
      common/meas/profmeas(2,2500,38),smeas(2,2500)  
      common/scale/f1,f2
      logical center, between
      dimension t(38),t32(38)
      call hbook1(11,'spectrum CO(3-2)mod',38.,0.,38.,0.)
      call hbook2(21,'mod',50.,0.,50.,50.,0.,50.,0.) 
      call hbook2(22,'map CO(3-2)',50.,0.,50.,50.,0.,50.,0.) 
      call hbook2(23,'map CO(6-5)',50.,0.,50.,50.,0.,50.,0.) 
c$$$      call hbook1(31,'check',38.,0.,38.,0.)
c$$$      call hbook1(32,'check',38.,0.,38.,0.)
c$$$      call hbook1(33,'check',38.,0.,38.,0.)
c$$$      call hbook1(34,'check',38.,0.,38.,0.)
      call hbook1(41,'check',50.,0.,50.,0.)
      call hbook1(42,'check',50.,0.,50.,0.)
      call hbook1(43,'check',50.,0.,50.,0.)
      call hbook1(44,'check',50.,0.,50.,0.)
      call hbook1(45,'check',38.,0.,38.,0.)
      call hbook1(46,'check',38.,0.,38.,0.)
      do ii=401,450
         call hbook1(ii,'',38.,0.,38.,0.)
      enddo


      center=.false.
      between=.false.
      do ij=1,2500
         j=ij/50+1
         i=ij-50*(j-1)    
         do k=1,38
            call hf2(21,float(i)-.5,float(j)-.5,profmod(ij,k))
            call hf1(11,float(k)-.5,profmod(ij,k))
         enddo
      enddo
         do kch=1,25
            mm=(kch-1)/5
            ii1=10*mm+1
            if(center)ii1=2*mm+21
            if(between)ii1=6*mm+11
            ii2=ii1+9
            if(center)ii2=ii1+1
            if(between)ii2=ii1+5
            nn=kch-5*mm-1
            jj1=10*nn+1
            if(center)jj1=2*nn+21
            if(between)jj1=6*nn+11
            jj2=jj1+9               
            if(center)jj2=jj1+1
            if(between)jj2=jj1+5
c            print*,kch,jj1,jj2,ii1,ii2
            do kk=1,38
               t(kk)=0.
               t32(kk)=0.
                do ij=1,2500
                   j=ij/50+1
                   i=ij-50*(j-1)  
                if(j.ge.ii1.and.j.le.ii2.and.i.ge.jj1.and.i.le.jj2)then
                   t(kk)=t(kk)+f1*profmod(ij,kk)
                   t32(kk)=t32(kk)+profmeas(1,ij,kk)
                endif
                enddo
               call hf1(400+kch,float(kk)-.5,t(kk))
               call hf1(425+kch,float(kk)-.5,t32(kk))
            enddo
         enddo

      return
      end


      subroutine readdata
      common/model/profmod(2500,38),garea(1000)
      common/meas/profmeas(2,2500,38),smeas(2,2500)  
      common/scale/f1,f2
c     
c read CO(3-2) data and store in co32
c ii, jj and kk are simply the numbers in the datacubes, incremented by 1
c
      open(2,file='co32-rebin3.txt')
      open(3,file='co65-rebin3.txt')
      sum=0.
      sum1=0.
      sum2=0.
      chi1=0.
      chi2=0.
      do i=1,50
         do j=1,50
            ij=i+50*(j-1)
            smeas(1,ij)=0.
            smeas(2,ij)=0.
            do k=1,38              
               read(2,*) i1,j1,k1,aa
               read(3,*) i2,j2,k2,bb
               if(i.ne.i1.or.i.ne.i2.or.j.ne.j1.or.j.ne.j2
     &              .or.k.ne.k1.or.k.ne.k2)then
                  print*,'read file problem: ',i,j,k,i1,j1,k1,i2,j2,k2
                  stop
               endif
               profmeas(1,ij,k)=aa
               profmeas(2,ij,k)=bb
c$$$               if(ij.eq.175)call hf1(31,float(k)-.5,profmod(ij,k))
c$$$               if(ij.eq.175)call hf1(32,float(k)-.5,profmeas(1,ij,k))
c$$$               if(ij.eq.1285)call hf1(33,float(k)-.5,profmod(ij,k))
c$$$               if(ij.eq.1285)call hf1(34,float(k)-.5,profmeas(1,ij,k))
               if(i.gt.10..and.i.lt.40..and.j.gt.15..and.j.lt.35.)goto 1
               if (i.ge.25.and.i.le.26.and.j.ge.25.and.j.le.26)goto 1
               smeas(1,ij)=smeas(1,ij)+aa
               smeas(2,ij)=smeas(2,ij)+bb
               sum1=sum1+aa
               sum2=sum2+bb
               sum=sum+profmod(ij,k)
 1             continue
            enddo
            call hf2(22,float(i)-.5,float(j)-.5,smeas(1,ij))
            call hf2(23,float(i)-.5,float(j)-.5,smeas(2,ij))               
         enddo
      enddo
      f1=sum1/sum
      f2=sum2/sum
      do ij=1,2500
         j=ij/50+1
         i=ij-50*(j-1) 
         do k=1,38
            if(i.gt.10..and.i.lt.40..and.j.gt.15..and.j.lt.35.)then
               if (i.ge.25.and.i.le.26.and.j.ge.25.and.j.le.26)goto 2
            error2=f1**2*profmod(ij,k)**2+profmeas(1,ij,k)**2
            if(error2.gt.0.)then
            chi1=chi1+(profmeas(1,ij,k)-f1*profmod(ij,k))**2
     &           /error2
            endif
            chi2=chi2+(profmeas(2,ij,k)-f2*profmod(ij,k))**2
            call hf1(41,float(i)-.5,profmeas(1,ij,k))
            call hf1(42,float(i)-.5,f1*profmod(ij,k))
            call hf1(43,float(j)-.5,profmeas(1,ij,k))
            call hf1(44,float(j)-.5,f1*profmod(ij,k))
            call hf1(45,float(k)-.5,profmeas(1,ij,k))
            call hf1(46,float(k)-.5,f1*profmod(ij,k))
         endif
 2          continue
         enddo
      enddo
      print*,'chi1, chi2: ', chi1/1.e3,chi2/1.e3

      close(2)
      close(3)
      return
      end

      subroutine test
      real par(10)
      data par/1.,14.,.87,.9,13.,.5,1.5,1.,.05,3.3/
      data rmin /.05/
      pi=2.*acos(0.)
c     
      theta=par(1)
      omega=par(2)
      cbicone=par(3)
      vdisk=par(4)
      vbicone=par(5)
      vrot=par(6)
      rrot=par(7)
      adisk=par(8)
      delta=par(9)
      power=par(10)

c$$$      do ir=1,50
c$$$         call hbook1(100+ir,'rho vs phi',20.,0.,1.,0.)
c$$$         call hbook1(200+ir,'v vs phi',20.,0.,1.,0.)
c$$$      enddo
c$$$      do iv=1,20
c$$$         call hbook1(150+iv,'rho vs r',50.,0.,4.,0.)
c$$$         call hbook1(250+iv,'v vs r',50.,0.,4.,0.)
c$$$      enddo
      call hbook2(300,'',30.,0.,2.5,30.,0.,2.5,0.)
      call hbook2(301,'',30.,0.,2.5,30.,0.,2.5,0.)

      do ir=1,50
         r=float(ir)*.08-.04
         do icphi=1,20
            cphi=float(icphi-1)/20.+.05
            if(abs(cphi).gt.cbicone)goto 1
            v0=(abs(cphi)/cbicone)**power*(vbicone-vdisk)+vdisk
            if(r.lt.rmin)v=v*(r/rmin)**2
            v=v0
            dv=0.
            if(r.ge.rrot)rho=1./4./pi*(rrot/r)**2
            if(r.lt.rrot)then
               shape=exp(-.5*(r*cphi/delta)**2)
               rho=1/8./pi*(3.-(r/rrot)**2)*(1.+adisk*shape)
               dv=vrot*shape/sqrt(r/rrot)
               if(r.lt.rmin)dv=dv*r/sqrt(rmin*rrot)
            endif           
            call hf1(100+ir,cphi,rho)
            call hf1(200+ir,cphi,v)
            call hf1(150+icphi,r,rho)
            call hf1(250+icphi,r,v)
            call hf2(300,r*abs(cphi),r*sqrt(1.-cphi**2),dv)
            call hf2(301,r*abs(cphi),r*sqrt(1.-cphi**2),rho)
 1          continue
         enddo
      enddo

      return
      end
