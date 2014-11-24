      program disk
      real flux(200,200,200)
      pi=2.*acos(0.)
      print*,'pi: ',pi
c
c a disk has thickness of z0 = .2, radius r0 = 1., z rotation axis, 
c omega *r0 = 10 km/s
c y east, z north, x light of sight
c pixel size .05 x .05
c each (r,phi,z) has corresponding pixel (y,z)
      call HLIMAP(20000,'SPEC')
c      call hbook1(10,'spec',19.,-10.,9.,0.)
      call hbook1(10,'spec',200.,-100.,100.,0.)
c      call hbook1(10,'spec',40.,-4.,4.,0.)
      call hbook1(11,'iv',200.,0.,200.,0.)
      call hbook1(12,'v',200.,-100.,100.,0.)
      call hbook2(20,'map',40.,-20.,20.,20.,-10.,10.,0.)
      call hbook2(21,'map',40.,-20.,20.,20.,-10.,10.,0.)

      r0=.5
      dv=1.
      z0=.2
      rho=1.
      xdelta=.05
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
         z=(float(iz)-100.5)*.01
         do iy=1,200
            y=(float(iy)-100.5)*.01
            do ix=1,200
               x=(float(ix)-100.5)*.01
               r2=x**2+y**2
               r=sqrt(r2)
               phi=atan2(y,x)
               if(abs(r).le.r0.and.abs(z).le.z0)then
                  v02=(50.*r/r0)
                  v0=50.
                  v=v0*cos(phi)+v02*sin(phi)
                  call hf1(12,v,1.)
                  iv=ifix(v/dv)+100
                  if(iv.gt.ivmax)ivmax=iv
                  if(iv.lt.ivmin)ivmin=iv
                  flux(iy,iz,iv)=flux(iy,iz,iv)+rho*xdelta
                  call hf1(11,float(iv)-.5,1.)
               endif
               enddo
            enddo
         enddo

         print*,ivmin,ivmax

      do i=1,200
         do j=1,200
            do k=1,200
               call hf1(10,float(k)-100,flux(i,j,k))
               call hf2(20,float(i)-50.5,float(j)-50.5,flux(i,j,k))
               call hf2(21,float(i)-50.5,float(k)-10.5,flux(i,j,k))
            enddo
         enddo
      enddo
 
      stop
      end
