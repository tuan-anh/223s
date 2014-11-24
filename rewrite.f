      program COmodel  
      external fcn       
      real*8 chi2
      real*8 pstart(30),pstep(30),plow(30),pup(30),pfit(30),errp(30)
      real profmod1(1001),profmod2(1001),smod(169),garea(1000)
      real ymeas(169),zmeas(169),profmeas1(169,101),profmeas2(169,101)
      real par(30)
      real pline(1001)
      real xpar(30),ifit(30)
      character*20 name(30),xname(30)
      integer*4 timeArray(3)
      common/para/npar,par
      common/meas/nmeas,ymeas,zmeas,profmeas1,profmeas2,vmin,vmax,weight
      common/model/kv,dv,xdelta,profmod1,profmod2,smod,garea
      common/profline/pline
      common/constant/pi,frac1,frac2,fJy,jup1,jup2,eup1,eup2,freq1,freq2
      data nmeas,vmin,vmax,kv,dv/169,-10.,10.,1001,0.02/
      data xdelta/0.05/      
      call itime(timeArray)
      print*,'start',timeArray(1),'h',timeArray(2),'m',timeArray(3),'s' 
      pi=2.*acos(0.)
!!! input
      dstar=143.                !! pc
      fracH=0.7
      fracCO=4.E-4
!! CO10
      jup1=1                    !! upper level
      eup1=5.5                  !! 1/K      
      probA1=7.4E-8
      freq1=115.3               ! GHz
      weight1=1./0.014**2
      open(1,file="CO10-2014.dat") !! input 
!! CO21
      jup2=2
      eup2=16.6                     
      probA2=7.1E-7
      freq2=230.5     
      weight2=1./(0.014*7.9)**2
      open(2,file="CO21-2014.dat")
      open(3,file="CO1n2-model-sym.dat")
      open(4,file="CO1n2-model-star.dat")
c      open(3,file="CO1n2-model-asym-offset.dat")
            
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! converse the unit of the density to m^-3. 
!     10^-7*Ms/mp/(pi*10^7)/(dstar*asec2m)^2/km2m
!! 1 arcsec=dstar*0.15E12 m

      xm2dens=1.41E17/dstar**2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      npar=24
      do i=1,npar
         pstart(i)=0.
      enddo

!!!! CO1n2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!
c asym fit
c$$$      pstart( 1)=   0.91396
c$$$      pstart( 2)=   0.18068
c$$$      pstart( 3)=   0.331
c$$$
c$$$      pstart( 4)=   2.582
c$$$      pstart( 5)=   4.833
c$$$      pstart( 6)=   0.215
c$$$      pstart( 7)=   0.108
c$$$      pstart( 8)=   0.53
c$$$      pstart( 9)=   0.90
c$$$
c$$$      pstart(10)=   161
c$$$      pstart(11)=   1.50
c$$$      pstart(12)=   0.88
c$$$      pstart(24)=   1.18
c$$$
c$$$      pstart(13)=   .626
c$$$      pstart(14)=   .296
c$$$      pstart(15)=   2.5
c$$$      pstart(16)=   2.5
c$$$
c$$$      pstart(20)=   -0.96
c$$$      pstart(21)=   0.95
c$$$      pstart(22)=   -0.40
c$$$      pstart(23)=   0.69

c scan
c$$$      pstart( 1)=   0.920
c$$$      pstart( 2)=   0.15
c$$$      pstart( 3)=   0.33
c$$$
c$$$      pstart( 4)=   2.58
c$$$      pstart( 5)=   4.83
c$$$      pstart( 6)=   0.22
c$$$      pstart( 7)=   0.11
c$$$      pstart( 8)=   0.53
c$$$      pstart( 9)=   0.89
c$$$
c$$$      pstart(10)=   162
c$$$      pstart(11)=   1.6
c$$$      pstart(12)=   0.88
c$$$      pstart(24)=   1.18
c$$$
c$$$      pstart(13)=   .64
c$$$      pstart(14)=   .30
c$$$      pstart(15)=   2.5
c$$$      pstart(16)=   2.5
c$$$
c$$$      pstart(20)=   -0.96
c$$$      pstart(21)=   0.95
c$$$      pstart(22)=   -0.40
c$$$      pstart(23)=   0.69

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c sym fit
      pstart( 1)=   0.91396
      pstart( 2)=   0.18068
      pstart( 3)=   0.331

      pstart( 4)=   2.56
      pstart( 5)=   4.91
      pstart( 6)=   0.21
      pstart( 7)=   0.11

      pstart( 8)=   0.52
      pstart( 9)=   0.91

      pstart(10)=   98.
      pstart(11)=   1.12
      pstart(12)=   0.
      pstart(24)=   0.

c$$$      pstart(10)=   84.5
c$$$      pstart(11)=   0.942
c$$$      pstart(12)=  13.1

      pstart(13)=   .62
      pstart(14)=   .36
      pstart(15)=   2.5
      pstart(16)=   2.5

      pstart(20)=   0.
      pstart(21)=   0.
      pstart(22)=   0.
      pstart(23)=   0.

      do i=1,npar
         pstep(i)=0.05
      enddo

      pstep(17)=1.E-17      
      
cccc calculate the array garea(1000) 
      x=-2.997
      garea(1)=.006*exp(-0.5*x**2)
      do ig=2,1000
         x=x+0.006
         garea(ig)=garea(ig-1)+.006*exp(-0.5*x**2)
      enddo
      do ig=1,1000
         garea(ig)=garea(ig)/garea(1000)
      enddo
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i=1,nmeas
         do k=1,101
            read(1,*) ii,v,profmeas1(ii,k)
            read(2,*) ii,v,profmeas2(ii,k)
         enddo
      enddo
      frac1=xm2dens*2.37E-18*fracH*fracCO*probA1*dstar/dv  !! convert the intensity 
      frac2=xm2dens*2.37E-18*fracH*fracCO*probA2*dstar/dv !! convert the intensity 
      print*,frac,xm2dens
                                                       !! =xm2dens*h[f]/4/pi/[deltaf]*fracH*fracCO*probA*[deltax]
      fJy=0.24E16                                      !! convert the flux density to Jy
                                                       !! = [psize**2]/dstar**2*E26

      close(1)
      close(2)
      do ide=1,13
         do ira=1,13
            k=(ide-1)*13+ira
            ymeas(k)=-float(ira-7)*1.4
            zmeas(k)=-float(ide-7)*1.4          
         enddo
      enddo
      call mninit(5,6,7)
      call mncomd(fcn,"SET PRI 0",ierr,0)         

      call mnparm(1,'x1',pstart(1),pstep(1),0.d0,1.57d0,ierr)
      call mnparm(2,'x2',pstart(2),pstep(2),0.d0,1.57d0,ierr)
      call mnparm(3,'x3',pstart(3),pstep(3),0.d0,1.d0,ierr)
      call mnparm(4,'x4',pstart(4),pstep(4),0.1d0,10.d0,ierr)
      call mnparm(5,'x5',pstart(5),pstep(5),0.1d0,10.d0,ierr)
c      call mnparm(6,'x6',pstart(6),pstep(6),0.d0,1.d0,ierr)
c      call mnparm(7,'x7',pstart(7),pstep(7),0.d0,1.d0,ierr)
      call mnparm(6,'x6',pstart(6),pstep(6),0.d0,5.d0,ierr)
      call mnparm(7,'x7',pstart(7),pstep(7),0.d0,5.d0,ierr)

      call mnparm(8,'x8',pstart(8),pstep(8),0.d0,10.d0,ierr)
      call mnparm(9,'x9',pstart(9),pstep(9),0.d0,10.d0,ierr)
      call mnparm(10,'x10',pstart(10),pstep(10),0.d0,0.d0,ierr)
      call mnparm(11,'x11',pstart(11),pstep(11),-2.d0,0.d0,ierr)
      call mnparm(12,'x12',pstart(12),pstep(12),-1.d0,1.d0,ierr)
      call mnparm(13,'x13',pstart(13),pstep(13),0.d0,0.999d0,ierr)
      call mnparm(14,'x14',pstart(14),pstep(14),0.d0,0.999d0,ierr)
      call mnparm(15,'x15',pstart(15),pstep(15),0.1d0,3.d0,ierr)
      call mnparm(16,'x16',pstart(16),pstep(16),0.1d0,3.d0,ierr)
      call mnparm(17,'x17',pstart(17),pstep(17),0.d0,0.d0,ierr)
      call mnparm(18,'x18',pstart(18),pstep(18),0.d0,0.d0,ierr)
      call mnparm(19,'x19',pstart(19),pstep(19),0.d0,0.d0,ierr)
      call mnparm(20,'x20',pstart(20),pstep(20),0.d0,0.d0,ierr)
      call mnparm(21,'x21',pstart(21),pstep(21),0.d0,0.d0,ierr)
      call mnparm(22,'x22',pstart(22),pstep(22),-.999d0,.999d0,ierr)
      call mnparm(23,'x23',pstart(23),pstep(23),-.999d0,.999d0,ierr)
      call mnparm(24,'x24',pstart(24),pstep(24),0.d0,5.d0,ierr)

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
      call mncomd(fcn,'FIX 15',ierr,0)
      call mncomd(fcn,'FIX 16',ierr,0)
      call mncomd(fcn,'FIX 17',ierr,0)
      call mncomd(fcn,'FIX 18',ierr,0)
      call mncomd(fcn,'FIX 19',ierr,0)
      call mncomd(fcn,'FIX 20',ierr,0)
      call mncomd(fcn,'FIX 21',ierr,0)
      call mncomd(fcn,'FIX 22',ierr,0)
      call mncomd(fcn,'FIX 23',ierr,0)
      call mncomd(fcn,'FIX 24',ierr,0)
     
      call mncomd(fcn,'MINIMIZE 1000',ierr,0)    
      do i=1,npar
         call mnpout(i,name(i),pfit(i),errp(i),bnd1,bnd2,ivarbl)
         xpar(i)=pfit(i)
      enddo
      call mnstat(chi2,fedm,errdef,npari,nparx,istat)          
      call mncomd(fcn,'RET',ierr,0)

      do m=1,npar
         print 999,'pstart(',m,')=',xpar(m)
      enddo
 999  format(6X,A7,I2,A2,f10.5)

      call itime(timeArray)
      print*,"finished at",timeArray(1),"h",timeArray(2),"min"
     n     ,timeArray(3),"s"

      stop 
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fcn(npara,grad,chi2,xval,iflag,futil)
      real*8 xval,grad,chi2
      dimension grad(*),xval(*)
      integer kv
      real dv
      real profmod1(1001),profmod2(1001),smod(169),garea(1000)
      real ymeas(169),zmeas(169),profmeas1(169,101),profmeas2(169,101)
      real par(30)
      common/meas/nmeas,ymeas,zmeas,profmeas1,profmeas2,vmin,vmax,weight
      common/model/kv,dv,xdelta,profmod1,profmod2,smod,garea
      common/profline/pline(1001)
      common/para/npar,par
      common /wind/v,rho
      common/constant/pi,frac1,frac2,fJy,jup1,jup2,eup1,eup2,freq1,freq2
      data alp1s,alp2s/0.1,20./
      dimension prof1(101),prof2(101)  

c nmeas is the number of measured velocity spectra
c ymeas, zmeas are the coordinates of the source associated with each spectrum 
c all coordinates in arcsecons centered on the centre of the star
c profmeas are the measured spectra, smeas are their sums over bins 
c they go from v=vmin to vmax=vmin+kv*dv in kv steps of dv
c errmeas is the measurement uncertainty (a single value for all measurements ..? 
c up to us to define) 
c all these data have been read earlier and stored in common /meas/
c sigma(x,y)beam are the sigmas of the (Gaussian) beam
c alp1s and alp2s define the range in r over which we work
c profmod is the spectrum predicted by the model for the current position
c xdelta is the scan step in x
c the constants in common /constant/ have been stored earlier

      do i=1,npar
         par(i)=xval(i)
      enddo
      chi2=0.
      chi21=0.
      chi22=0.
      do 1 imeas=1,169
         do iprof=1,kv
            profmod1(iprof)=0.
            profmod2(iprof)=0.
         enddo
         ide=(imeas-1)/13+1
         ira=mod(imeas-1,13)+1
         if(ide.ge.4.and.ide.le.10.and.ira.ge.4.and.ira.le.10)then
c         if(ide.ge.1.and.ide.le.13.and.ira.ge.1.and.ira.le.13)then
            npix=10
            psize=1.4/10.
            if(abs(ymeas(imeas)).le.1..and.abs(zmeas(imeas)).le.1.)then
               npix=10
               psize=1.4/10.
            endif
            do i=1,npix
               do j=1,npix
                  do icase=1,2                     
                  y=ymeas(imeas)+(float(i)-float(npix+1)/2.)*psize
                  z=zmeas(imeas)+(float(j)-float(npix+1)/2.)*psize
                  xmax2=alp2s**2-z**2-y**2
                  xmax=sqrt(xmax2)             
                  nmax=2*ifix(xmax/xdelta)+1
                  xmax=float(nmax-1)*xdelta/2.
                  
                  do iprof=1,kv
                     pline(iprof)=0.
                  enddo
                  do 4 ix=1,nmax
                     x=xmax-float(ix-1)*xdelta
                     rchk=sqrt(x**2+y**2+z**2)
                     if(rchk.lt.alp1s) goto 4
                     call physics(x,y,z,imeas,icase)
 4                continue
                  do iprof=1,kv
                     if(icase.eq.1)profmod1(iprof)=profmod1(iprof)
     n                    +pline(iprof)
                     if(icase.eq.2)profmod2(iprof)=profmod2(iprof)
     n                    +pline(iprof)
                  enddo
               enddo
            enddo
            enddo

            do iv=1,kv
                profmod1(iv)=profmod1(iv)*fJy*psize**2
                profmod2(iv)=profmod2(iv)*fJy*psize**2
            enddo
            do iv=1,101
               prof1(iv)=0.
               prof2(iv)=0.
            enddo
             iw=ifix(0.2/dv/2.)
             do iv=1,101
                v=float(iv-1)*0.2-10.
                ih=(iv-1)*2*iw+iw+1
                il=(iv-1)*2*iw-iw+1
                n=0
                if(il.lt.1) il=1
                if(ih.gt.kv) ih=kv
                do iprof=il,ih
                   n=n+1
                   prof1(iv)=prof1(iv)+profmod1(iprof)
                   prof2(iv)=prof2(iv)+profmod2(iprof)
                enddo
                if(n.eq.0) print*,'!!! n=0 ',iprof,il,ih,iw,iv
                prof1(iv)=prof1(iv)/float(n)
                prof2(iv)=prof2(iv)/float(n)
             enddo
             do iv=1,101
                err1sq=(0.003**2+(0.02*profmeas1(imeas,iv))**2)*21.
                err2sq=(0.03**2+(0.02*profmeas2(imeas,iv))**2)*15.
                chi21=chi21+(prof1(iv)-profmeas1(imeas,iv))**2/err1sq
                chi22=chi22+(prof2(iv)-profmeas2(imeas,iv))**2/err2sq
                chi2=chi21+chi22
                if(iflag.eq.3)then
                  v=float(iv-1)*0.2+vmin
                  write(3,*) imeas,v,prof1(iv),prof2(iv)
               endif
            enddo
         endif
 1    continue
      print 998,chi2/4949/2.,par(1),par(2),par(3),par(4),par(5),
     n     par(6),par(7),par(8),par(9),
     n     par(10),par(11),par(12),par(24),par(13),par(14),
     n     par(20),par(21),par(22),par(23)
      if(iflag.eq.3) print*,'chi2 = ',chi2,chi21/4949.,chi22/4949.
 998  format(26(f7.3,1x))
      return
      end 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   
      subroutine physics(x,y,z,n,icase)
      real par(30)
      real cphi,r,theta,omega
      common/para/npar,par
      common /wind/v,rho
      common/constant/pi,frac1,frac2,fJy,jup1,jup2,eup1,eup2,freq1,freq2
      if(icase.eq.1) then
         jup=jup1
         eup=eup1
         frac=frac1
         freq=freq1
      endif
      if(icase.eq.2) then
         jup=jup2
         eup=eup2
         frac=frac2
         freq=freq2
      endif
         
      r2=x**2+y**2+z**2
      r=sqrt(r2)
      fCO=0.5**((r/11.)**2.5)

c$$$      if(r.gt.2.) t=par(10)*4.**(-par(12))*(r/2.)**(-par(24))
c$$$      if(r.ge.0.5.and.r.lt.2.) then
c$$$         t=par(10)*(r/.5)**(-par(12))
c$$$      endif
c$$$      if(r.lt.0.5) then
c$$$         t=par(10)*(r/.5)**(-par(11))
c$$$      endif

      theta=par(1)+pi/2.
      omega=pi/2.-par(2)
      ct=cos(theta)
      st=sin(theta)
      co=cos(omega)
      so=sin(omega)
      cphi=(x/r)*ct+(y/r)*st*co+(z/r)*st*so 
      calpha=x/r
ccc gaussian distribution
      shape=exp(-(1.-cphi)**2/2./par(3)**2)+
     n     exp(-(1.+cphi)**2/2./par(3)**2)

c$$$      if(r.lt.par(24)) then
c$$$         t=par(10)*(r/par(24))**par(11)
c$$$      endif
c$$$      if(r.ge.par(24)) then
c$$$         t=par(10)*(r/par(24))**par(11)*
c$$$     n        (1.-par(12)*abs(cphi)*(log10(r/par(24)))**2)
c$$$      endif
      t=par(10)*r**(-par(11))
      if(t.lt.2.) t=2.

      prob_emis=(2.*float(jup)+1)*exp(-eup/t)*2.8/t

c      v0=par(5)*shape*(1.+par(22)*cphi)*r**par(7)
c     n     +par(4)*(1.+par(23)*cphi)*r**par(6)
c      v0=par(5)*shape*(1.+par(22)*cphi)*(1.-exp(-r/par(7)))
c     n     +par(4)*(1.+par(23)*cphi)*(1.-exp(-r/par(6)))

      v0=par(5)*shape*(1.+par(22)*cphi)*(1.-par(14)*exp(-r/par(16)))
     n     +par(4)*(1.+par(23)*cphi)*(1.-par(13)*exp(-r/par(15)))
      if(v0.lt.0.001) v0=0.001

c$$$      rc=6.
c$$$      if(r.lt.rc) then
c$$$      v0=par(5)*shape*(1.+par(22)*cphi)*(1.-((r-rc)/par(13))**2)
c$$$     n     +par(4)*(1.+par(23)*cphi)*(1.-((r-rc)/par(14))**2)
c$$$      endif      
c$$$      if(r.ge.rc) then
c$$$      v0=par(5)*shape*(1.+par(22)*cphi)
c$$$     n     +par(4)*(1.+par(23)*cphi)
c$$$      endif

      xmloss=par(9)*(1.+par(20)*cphi)*shape+par(8)*(1.+par(21)*cphi)
c      dens=xmloss/4./pi/r**2/v0*xm2dens*1.E-6
      rho=xmloss/4./pi/r**2/v0*fCO*prob_emis*frac
      wei=rho
      if(n.eq.85) wei=rho
      if(abs(x).lt..5) write(4,*) n,icase,x,y,z,r,wei,xmloss,v0
c      if(icase.eq.1) write(11,*) r,wei
c      if(icase.eq.2) write(12,*) r,wei
      v=v0*calpha
      sigmav1=0.02*sqrt(t)
      sigmav2=0.2
      sigmav=sqrt(sigmav1**2+sigmav2**2)
      fabs=0.68E23/freq**3*(exp(eup/t)-1.) !! c^2/2/h/f^3*(...)
      call fillprofmod(v,sigmav,rho,fabs)
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      subroutine fillprofmod(v,sigmav,rho,fabs)
      integer kv
      real dv
      real profmod1(1001),profmod2(1001),smod(169),garea(1000)
      real ymeas(169),zmeas(169),profmeas1(169,101),profmeas2(169,101)
      common/meas/nmeas,ymeas,zmeas,profmeas1,profmeas2,vmin,vmax,weight
      common/model/kv,dv,xdelta,profmod1,profmod2,smod,garea
      common/profline/pline(1001)
      common/constant/pi,frac1,frac2,fJy,jup1,jup2,eup1,eup2,freq1,freq2

      fkv=float(kv)
      v1=(v-3.*sigmav-vmin)/dv
      v2=(v+3.*sigmav-vmin)/dv
      if(v2.lt.0.) goto 2
      if(v1.gt.fkv) goto 2
c     the array garea(1000) contains the integral of a gaussian of area=1, 
c     sigma=1 and mean=0 in 1000 bins between -3 and 3
c     Fill bin ibin from vlo to vhi
      if(v1.lt.0.)then 
         ibin=1 
         vlo=0.
         vhi=1. 
      endif
      if(v1.ge.0.) then
         ibin=ifix(v1)+1 
         vlo=v1
         vhi=float(ibin) 
      endif

      if(v2-v1.lt.0.0005) then
         if(ibin.ge.1.and.ibin.le.kv) then
            if((1.-fabs*pline(ibin)).gt.0.) then
               pline(ibin)=pline(ibin)+rho*xdelta
     n              *(1.-fabs*pline(ibin))
            endif
         endif
         goto 2
      endif

 1    continue
      glo=1000.*(vlo-v1)/(v2-v1)
      jlo=ifix(glo)+1      
      if(jlo.gt.1) dglo=(garea(jlo)-garea(jlo-1))*(float(jlo)-glo)
      if(jlo.eq.1) dglo=garea(1)*(float(jlo)-glo)
      if(v2.le.vhi) then
         ghi=1000
         jhi=1000
      endif
      if(v2.gt.vhi) then
         ghi=1000.*(vhi-v1)/(v2-v1)
         jhi=ifix(ghi)+1
      endif
      dghi=(garea(jhi)-garea(jhi-1))*(float(jhi)-ghi)
      if((1.-fabs*pline(ibin)).gt.0.) then
         pline(ibin)=pline(ibin)+(garea(jhi)-garea(jlo)
     n        +dglo-dghi)*rho*xdelta*(1.-fabs*pline(ibin))
      endif
c     garea(j) is the integral from -3 to -3+.006*j
      vlo=vhi
      if(vlo.ge.v2) goto 2
      ibin=ibin+1
      if(ibin.gt.kv) goto 2
      vhi=float(ibin) 
      go to 1
 2    continue
      
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function tem(r) !! r in arcsec
!!!! 1 arcsec=0.15*d[in pc] *1.E14 cm
!!! RSCnc
!!! T=T0/r[in cm]**alp  with T0=10**(14*alp+3.05) and alp=0.7  
  
c      tem=650./(r/0.1)**0.7
      tem=580./(r/0.1)**0.7
!!!
      return
      end
