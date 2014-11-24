      program disk3d
      real flux(200,200,200)
      pi=2.*acos(0.)
      print*,'pi: ',pi
c
c a spherical star has a radius r0 = .5", rmin = .2"
c origin of coordinate is at the centre of star
c z rotation axis, to the North; y to the East; x light of sight
c rotation velocity at r0: vr0 = 100 km/s 
c disk radial velocity vd_exp = 10 km/s in (x,y) plane
c pixel size .01 x .01
c each (r,phi,z) has corresponding pixel (y,z)
c for each (ix,iy,iz) I want to project to light of sight direction
c each pixel (iy,iz) has a spectrum integrated over light of sight
      call HLIMAP(20000,'SPEC')
      call hbook1(9,'spec',20.,-200.,200.,0.)
      call hbook1(10,'spec',20.,-100.,100.,0.)
      call hbook1(11,'iv',20.,0.,200.,0.)
c      call hbook1(12,'v',30.,-150.,150.,0.)
c      call hbook2(20,'map',40.,.5,1.5,40.,.5,1.5,0.)
      do ik=1,11
       call hbook1(11+ik,'v',40.,-200.,200.,0.)        
      enddo

      do ik=1,6
c         call hbook2(20+ik,'map',40.,.5,1.5,40.,.5,1.5,0.)
      enddo
      r0=.5
      rmin=.2
      z0=.2
      dv=1.
      rho=1.
      vr0=100.
      vd_exp=10.
      step=.01
         do iy=1,200
            do iz=1,200
               do iv=1,200
               flux(iy,iz,iv)=0.
            enddo
         enddo
      enddo
      ivmax=0
      ivmin=1000
      do iz=1,200
         z=(float(iz)-100.5)*step
         do iy=1,200
            y=(float(iy)-100.5)*step
            do ix=1,200
               x=(float(ix)-100.5)*step
               r2=x**2+y**2+z**2
               r=sqrt(r2)
c calpha is projection of 0 direction
c or velocity on the light of sight
c in the case of spherical model
               calpha=-y*z/r2
!!               v0=vr0*sqrt(x**2+y**2)/r0
c phi is projection of velocity to the 
c light of sigh in case of cylinder model
               phi=atan2(y,x)
               rxy=sqrt(x**2+y**2)
               v0=vr0*sqrt(x**2+y**2)/r0
c disk expansion vd_exp = 10 km/s
        if(abs(rxy).le.r0.and.abs(rxy).ge.rmin.and.abs(z).lt.z0)then
!!                  v=v0*calpha
                  v=v0*sin(phi)+vd_exp*cos(phi)
               do iik=1,11
                  v=v0*sin(phi)+10.*(float(iik)-1.)*cos(phi)
                  call hf1(11+iik,v,1.)
               enddo
                  iv=ifix(v/dv+101.)
                  if(iv.gt.ivmax)ivmax=iv
                  if(iv.lt.ivmin)ivmin=iv
c                  flux(ix,iy,iv)=flux(ix,iy,iv)+rho
                  call hf1(11,float(iv)-.5,1.)
                  call hf1(9,v,rho)
               endif
            enddo
         enddo
      enddo
      print*,'ivmin, ivmax: ',ivmin,ivmax
      stop

      do i=1,200
         do j=1,200
            xx=(float(i)-.5)*.01
            yy=(float(j)-.5)*.01
            do k=1,200
               call hf1(10,float(k)-200.5,flux(i,j,k))
               call hf2(20,xx,yy,flux(i,j,k))
               iik=k/80+1
               call hf2(20+iik,xx,yy,flux(i,j,k))
            enddo
         enddo
      enddo
      stop
      end
