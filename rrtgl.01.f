      program COmodel  
c
c this version is for the Red rectangle using ALMA data
c
      external fcn       
      real*8 chi2
      real*8 pstart(30),pstep(30),pfit(30),errp(30)
      real par(30)
      real xpar(30),ifit(30)
c
c For CO(3-2)the beam is 500*490 mas^2, for CO(6-5), 270*240 mas^2
c namely 3.78 times smaller. 
c but Bujarrabal claims 480*440 and 310*250, namely 2.73 smaller: who is right?
c When calculating chisquare model is multiplied by beam size to be in Jy/beam.
c After having read the data using subroutine readdata,
c spectra include 40 velocity bins, 0.4 km/s wide, from -8 to 8 km/s
c sky maps are centered on the star in 50*50 bins of 0.1*0.1 arcsec**2
c from -2.5 to 2.5 arcsec; first index i increases from east to west,
c but y increases from west to east, we need to change its sign. 
c second index j increases from south to north; ij=i+50*(j-1), 1 to 2500.
c ico is 1 for CO(3-2) and 2 for CO(6-5)
c
      character*20 name(30)
      integer*4 timeArray(3)
      common/para/npar,par
      common/meas/profmeas(2,2500,40),smeas(2,2500)
      common/model/profmod(2,2500,40),smod(2,2500),garea(1000)
      common/constant/pi,frac(2),fJy,jup(2),eup(2),freq(2),
     &     sigmav,error(2),xdelta,beam(2)
      dimension probA(2) 
      data xdelta,sigmav,error/0.05,0.42,0.005,0.030/
      data freq,jup,eup,probA/345.562,691.474,3,6,33.2,116.2,
     &     7.4e-8,7.1e-7/
      data dstar,fracCO,fracH/710.,4.e-4,0.7/
      call HLIMAP(20000,'MAP')
      call itime(timeArray)
      print*,'start',timeArray(1),'h',timeArray(2),'m',timeArray(3),'s' 
      pi=2.*acos(0.)
      beam(1)=1.133*0.50*0.49
      beam(2)=1.133*0.27*0.24
! convert the unit of the density to m^-3. 
!     10^-7*Ms/mp/(pi*10^7)/(dstar*asec2m)^2/km2m
!! 1 arcsec=dstar*0.15E12 m
      xm2dens=1.41E17/dstar**2
      do ico=1,2
!! convert the intensity
!! =xm2dens*h[f]/4/pi/[deltaf]*fracH*fracCO*probA*[deltax]
         frac(ico)=xm2dens*2.37E-18*fracH*fracCO*probA(ico)*dstar/0.4 
         print*,'frac: ',frac,xm2dens
!! convert the flux density to Jy =[beamsize**2]/dstar**2*E26
         fJy=0.24E16
      enddo
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
c read CO data
c
      call readdata
c 
c minuit commands
c
      npar=14
      do i=1,npar
         pstart(i)=0.
      enddo
      pstart(1)=   5.
      pstart(2)=   11.
      pstart(3)=   0.03
      pstart(4)=   .1
      pstart(5)=   0.15
      pstart(6)=   .1
      pstart(7)=   100.
c mass loss rates (disk,bicone) in 10^-7 solar masses per year
      pstart(8)=   50.
      pstart(9)=   .001
c radial dependence of temperatures (disk,bicone)
      pstart(10)=   5.
      pstart(11)=   300.
      pstart(12)=   300.
      pstart(13)=   -2.
      pstart(14)=   -2.
      do i=1,npar
         pstep(i)=0.05*pstart(i)
      enddo
      pstep(1)=0.1
      pstep(2)=0.1
      pstep(4)=0.1
      pstep(10)=0.1
      pstep(13)=0.1
      pstep(14)=0.1
c
      call mninit(5,6,7)
      call mncomd(fcn,"SET PRI 0",ierr,0)         
c
      call mnparm(1,'x1',pstart(1),pstep(1),0.d0,0.d0,ierr)
      call mnparm(2,'x2',pstart(2),pstep(2),0.d0,0.d0,ierr)
      call mnparm(3,'x3',pstart(3),pstep(3),0.d0,0.d0,ierr)
      call mnparm(4,'x4',pstart(4),pstep(4),0.d0,0.d0,ierr)
      call mnparm(5,'x5',pstart(5),pstep(5),0.d0,0.d0,ierr)
      call mnparm(6,'x6',pstart(6),pstep(6),0.d0,0.d0,ierr)
      call mnparm(7,'x7',pstart(7),pstep(7),0.d0,0.d0,ierr)
      call mnparm(6,'x6',pstart(6),pstep(6),0.d0,0.d0,ierr)
      call mnparm(7,'x7',pstart(7),pstep(7),0.d0,0.d0,ierr)
      call mnparm(8,'x8',pstart(8),pstep(8),0.d0,0.d0,ierr)
      call mnparm(9,'x9',pstart(9),pstep(9),0.d0,0.d0,ierr)
      call mnparm(10,'x10',pstart(10),pstep(10),0.d0,0.d0,ierr)
      call mnparm(11,'x11',pstart(11),pstep(11),0.d0,0.d0,ierr)
      call mnparm(12,'x12',pstart(12),pstep(12),0.d0,0.d0,ierr)
      call mnparm(13,'x13',pstart(13),pstep(13),0.d0,0.d0,ierr)
      call mnparm(14,'x14',pstart(14),pstep(14),0.d0,0.d0,ierr)
c
      call mncomd(fcn,'FIX 1',ierr,0)
      call mncomd(fcn,'FIX 2',ierr,0)
      call mncomd(fcn,'FIX 3',ierr,0)
      call mncomd(fcn,'FIX 4',ierr,0)
      call mncomd(fcn,'FIX 5',ierr,0)
      call mncomd(fcn,'FIX 6',ierr,0)
      call mncomd(fcn,'FIX 7',ierr,0)
      call mncomd(fcn,'FIX 8',ierr,0)
      call mncomd(fcn,'FIX 9',ierr,0)
      call mncomd(fcn,'FIX 10',ierr,0)
      call mncomd(fcn,'FIX 11',ierr,0)
      call mncomd(fcn,'FIX 12',ierr,0) 
      call mncomd(fcn,'FIX 13',ierr,0)
      call mncomd(fcn,'FIX 14',ierr,0)
c     
c      call mncomd(fcn,'MINIMIZE 1000',ierr,0)
c    
      do i=1,npar
         call mnpout(i,name(i),pfit(i),errp(i),bnd1,bnd2,ivarbl)
         xpar(i)=pfit(i)
      enddo
      call mnstat(chi2,fedm,errdef,npari,nparx,istat)          
      call mncomd(fcn,'RET',ierr,0)
c 
      call itime(timeArray)
      print*,"finished at",timeArray(1),"h",timeArray(2),"min"
     n     ,timeArray(3),"s"
      call plot
      stop 
      end
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      subroutine fcn(npara,grad,chi2,xval,iflag,futil)
      real*8 xval,grad,chi2
      dimension grad(*),xval(*)
      real par(30)
      common/meas/profmeas(2,2500,40),smeas(2,2500)
      common/model/profmod(2,2500,40),smod(2,2500),garea(1000)
      common/para/npar,par
      common/constant/pi,frac(2),fJy,jup(2),eup(2),freq(2),
     &     sigmav,error(2),xdelta,beam(2)
      data rmin,rmax /0.1,4./
c
c all coordinates in arcseconds centered on the centre of the star
c profmeas are the measured spectra, smeas are their sums over velocity bins 
c they go from -8 to 8 km/s in 40 steps of 0.4 km/s
c sky coordinates go from -2.5 to 2.5 arcsec in 50 steps of 0.1 arcsec
c error is the measurement uncertainty (a single value for each 
c of CO(3-2) and CO(6-5) 
c rmin and rmax define the range in r over which we work
c profmod (same arguments as smeas) is the spectrum predicted by the model 
c xdelta is the scan step in x
c beam is the ratio between the CO(6-5) and CO(3-2) beams, 
c used to scale both sets of data to the CO(3-2) beam size     
      do i=1,npar
         par(i)=xval(i)
      enddo
      chi21=0.
      chi22=0.
      do 1 ij=1,2500
         do k=1,40
            profmod(1,ij,k)=0.
            profmod(2,ij,k)=0.
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
         do k=1,40
            chi21=chi21+((profmod(1,ij,k)-profmeas(1,ij,k))
     &           /error(1))**2/2500./40.
            chi22=chi22+((profmod(2,ij,k)-profmeas(2,ij,k))
     &           /error(2))**2/2500./40.
         enddo
    1 continue     
      chi2=chi21+chi22
      print*,'chi squares per dof',chi21,chi22,chi2
      print 998,(par(iprint),iprint=1,13)
      if(iflag.eq.3) print*,'final chi2 = ', chi21,chi22,chi2 
 998  format(13(f7.3,1x))
      return
      end 
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      
      subroutine physics(x,y,z,ij)
      real par(30),tag(2)
      common/para/npar,par
      common/model/profmod(2,2500,40),smod(2,2500),garea(1000)
      common/constant/pi,frac(2),fJy,jup(2),eup(2),freq(2),
     &     sigmav,error(2),xdelta,beam(2)
      data tag/1.,0./      
      r2=x**2+y**2+z**2
      r=sqrt(r2)
c	fCO is the UV dissociation of CO molecules at large r values
c	about .95 at 4 arcsec
      fCO=0.5**((r/11.)**2.5) 
c	theta, 5 degrees, angle between PPN axis and plane of the sky, par(1)
c	omega, 11 degrees, angle between north and projection of PPN axis on 
c	plane of the sky, measured eastward, par(2)
      theta=par(1)*pi/180.
      omega=par(2)*pi/180.
      ct=cos(theta)
      st=sin(theta)
      co=cos(omega)
      so=sin(omega)
      cphi=(x/r)*st+(y/r)*ct*so+(z/r)*ct*co 
c iphi=1 for the disk (abs(cphi)<.7), =2 for the bicone (else)
      iphi=1
      if(abs(cphi).ge.0.7)iphi=2
c calpha for projection of radial velocity on line of sight
      calpha=x/r
c cvel for projection of tangential velocity on line of sight
      velx=(y/r)*ct*co-(z/r)*ct*so
      vely=(z/r)*st-(x/r)*ct*co
      velz=(x/r)*ct*so-(y/r)*st
      vel=sqrt(velx**2+vely**2+velz**2)
      cvel=velx/vel
c par(3) is the sigma of the bicone wall in units of cphi, 0.08
c par(4) is the half aperture of the bicone in cphi
c par(5) is the sigma of the disk measured in cphi, 0.15
      shapeb=exp(-0.5*(abs(cphi)-1+par(4))**2/par(3)**2)
      shaped=exp(-0.5*(cphi)**2/par(5)**2)
      shape=tag(iphi)*shaped+(1.-tag(iphi))*shapeb
c par(6) is the expansion velocity of the disk in km/s,0.5
c par(7) is the expansion velocity of the bicone in km/s,50.
      v0=shape*par(5+iphi)
c par(8) is the mass loss rate on the disk
c par(9) is the mass loss rate on the bicone
      xmloss=shape*par(7+iphi)
c par(10) is rotation velocity of disk at r=1”
      v=v0*calpha      
c      if(iphi.eq.1) v=v+par(10)*(1./r)*cvel
      if(iphi.eq.1) v=par(10)*(1./r)*cvel
c     par(11,13) define temperature in the disk
c     par(12,14) define temperature in the bicone
      t=par(10+iphi)*r**(par(12+iphi))
      if(t.lt.2.) t=2.
      do ico=1,2
         prob_emis=(2.*float(jup(ico))+1)*exp(-eup(ico)/t)*2.8/t
         rho=xmloss/4./pi/r**2/v0*fCO*prob_emis*frac(ico)
c         if(r.lt..5.or.r.gt.1.5)rho=0.
         rho=rho*fJy*beam(ico)
c         fabs=0.68E23/freq(ico)**3*(exp(eup(ico)/t)-1.)
         fabs=0.
         call fillprofmod(v,rho,fabs,ico,ij)                 
      enddo
      return
      end

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      subroutine fillprofmod(v,rho,fabs,ico,ij)
      common/model/profmod(2,2500,40),smod(2,2500),garea(1000)
      common/constant/pi,frac(2),fJy,jup(2),eup(2),freq(2),
     &     sigmav,error(2),xdelta,beam(2)
c
c	fill array profmod in 40 bins of 0.4 km/s between -8 and 8 km/s
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
      if(vmin.lt.-8.) vmin=-8.
      if(vmin.gt.8.) go to 2
      vmax=v+3.*sigmav
      if(vmax.gt.8.) vmax=8.
      if(vmax.lt.-8.) go to 2
      if(vmin.gt.vmax) go to 2
      ivmin=ifix((vmin+8.)/0.4)+1
      ivmax=ifix((vmax+8.)/0.4)+1
      do iv=ivmin,ivmax
         vlo=-8.+0.4*float(iv-1)
         vhi=vlo+0.4
         weight=garea(ifix((vhi-v)/sigmav+3.)/.006)
     &        -garea(ifix((vlo-v)/sigmav+3.)/.006)
         profmod(ico,ij,iv)=profmod(ico,ij,iv)
     &        +rho*xdelta*weight*(1.-fabs*profmod(ico,ij,iv))
      enddo
 2    continue
      return
      end

      subroutine readdata
      real co32(70,360,360),co65(43,432,432)
      common/meas/profmeas(2,2500,40),smeas(2,2500)  
c
c read CO(3-2) data and store in co32
c ii, jj and kk are simply the numbers in the datacubes, incremented by 1
c
      open(2,file='12CO3-2_channel_image.txt')
      open(3,file='CO6-5_channel_image.txt')
      
      do k=1,70
         do j=1,360
            do i=1,360
               read(2,*) kk, jj, ii, aa
               co32(kk,jj,ii)=aa
            enddo
         enddo
      enddo
c
c read co65 data and store in co65
c ii, jj and kk are simply the numbers in the datacubes, incremented by 1
c we store in 44-kk rather than kk to have the velocity increase in the
c same direction for CO(3-2) and CO(6-5) data
c
      do k=1,43
         do j=1,432
            do i=1,432
               read(3,*) kk, jj, ii, aa
               co65(44-kk,jj,ii)= aa
            enddo
         enddo
      enddo
c
c convert velocity spectra in 40 bins of 0.4 km/s from -8 km/s to 8 km/s
c map pixels in 50 bins of .1 starting at -2.5 arcsec
c
      do i=1,50
         do j=1,50
            sco32=0.
            sco65=0.
            do k=1,40

               call transform(-181.5*.12,.12,i,-2.5,.1,floi,ilo)
               call transform(-177.5*.12,.12,j,-2.5,.1,floj,jlo)
               call transform(-36.5*.2115,.2115,k,-8.,.4,flok,klo)               
               profmeas(1,i+50*(j-1),k)=
     & floi*floj*(flok*co32(klo,jlo,ilo)+(1.-flok)*co32(klo+1,jlo,ilo))+
     & floi*(1.-floj)*(flok*co32(klo,jlo+1,ilo)+(1.-flok)
     & *co32(klo+1,jlo+1,ilo))+
     & (1.-floi)*floj*(flok*co32(klo,jlo,ilo+1)+(1.-flok)
     &             *co32(klo+1,jlo,ilo+1))+
     & (1.-floi)*(1.-floj)*(flok*co32(klo,jlo+1,ilo+1)+(1.-flok)
     &               *co32(klo+1,jlo+1,ilo+1))

               call transform(-220.5*.05,.05,i,-2.5,.1,floi,ilo)
               call transform(-210.5*.05,.05,j,-2.5,.1,floj,jlo)
               call transform(-20.7*.4,.4,k,-8.,.4,flok,klo)
               profmeas(2,i+50*(j-1),k)=
     & floi*floj*(flok*co65(klo,jlo,ilo)+(1.-flok)*co65(klo+1,jlo,ilo))+
     & floi*(1.-floj)*(flok*co65(klo,jlo+1,ilo)+(1.-flok)
     &               *co65(klo+1,jlo+1,ilo))+
     & (1.-floi)*floj*(flok*co65(klo,jlo,ilo+1)+(1.-flok)
     &               *co65(klo+1,jlo,ilo+1))+
     & (1.-floi)*(1.-floj)*(flok*co65(klo,jlo+1,ilo+1)+(1.-flok)
     &               *co65(klo+1,jlo+1,ilo+1))
		   sco32=sco32+profmeas(1,i+50*(j-1),k)
		   sco65=sco65+profmeas(2,i+50*(j-1),k)
             enddo
             smeas(1,i+50*(j-1))=sco32
             smeas(2,i+50*(j-1))=sco65
          enddo
      enddo
      return
      end

      subroutine transform(x0,a0,i1,x1,a1,flo,ilo)
c transform an array starting at x0 in bins of a0 into an array starting at x1 in bins of a1
c bin i1 will receive a content that is flo*the content of bin ilo+(1-flo)*the content of bin ilo+1
      x=x1+(float(i1)-.5)*a1
      xnew=(x-x0)/a0
      inew=ifix(xnew)+1
      dnew=xnew-float(inew)+.5
      if(dnew.ge.0.)then
         ilo=inew
         flo=1.-dnew 
      endif
      if(dnew.lt.0.) then
         ilo=inew-1
         flo=-dnew 
      endif
      return
      end
c
c
c
      subroutine plot
      common/meas/profmeas(2,2500,40),smeas(2,2500)
      common/model/profmod(2,2500,40),smod(2,2500),garea(1000)
      call hbook1(11,'spectrum CO(3-2)meas',40.,0.,40.,0.)
      call hbook1(12,'spectrum CO(6-5)meas',40.,0.,40.,0.)
      call hbook1(13,'spectrum CO(3-2)mod',40.,0.,40.,0.)
      call hbook1(14,'spectrum CO(6-5)mod',40.,0.,40.,0.)
      call hbook2(21,'map CO(3-2)meas',50.,0.,50.,50.,0.,50.,0.) 
      call hbook2(22,'map CO(6-5)meas',50.,0.,50.,50.,0.,50.,0.)
      call hbook2(23,'map CO(3-2)mod',50.,0.,50.,50.,0.,50.,0.) 
      call hbook2(24,'map CO(6-5)mod',50.,0.,50.,50.,0.,50.,0.)
      do ico=1,2
         do ij=1,2500
            j=ij/50+1
            i=ij-50*(j-1)    
            call hf2(20+ico,float(i)-.5,float(j)-.5,smeas(ico,ij))
            call hf2(22+ico,float(i)-.5,float(j)-.5,smod(ico,ij))
            do k=1,40
               call hf1(10+ico,float(k)-.5,profmeas(ico,ij,k))
               call hf1(12+ico,float(k)-.5,profmod(ico,ij,k))
            enddo
         enddo
      enddo
      return
      end



