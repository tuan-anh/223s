      program disk2d
      real flux(200,200,400)
      pi=2.*acos(0.)
      print*,'pi: ',pi
c
c a 2d disk has thickness has a radius r0 = .5, 
c z rotation axis normal to  (x,y) plane
c vr0 = 100 km/s and radial velocity vexp = 100 km/s
c x light of sight
c pixel size .01 x .01
c each (r,phi,z) has corresponding pixel (y,z)
      call HLIMAP(20000,'SPEC')
      call hbook1(9,'spec',20.,-200.,200.,0.)
      call hbook1(10,'spec',20.,-200.,200.,0.)
      call hbook1(11,'iv',20.,0.,400.,0.)
      call hbook1(12,'v',400.,-200.,200.,0.)
      call hbook2(20,'map',40.,.5,1.5,40.,.5,1.5,0.)
      do ik=1,6
         call hbook2(20+ik,'map',40.,.5,1.5,40.,.5,1.5,0.)
      enddo
      r0=.5
      rmin=.2
      dv=1.
      rho=1.
      vr0=100.
      vexp=100.
         do ix=1,200
            do iy=1,200
               do iv=1,400
               flux(ix,iy,iv)=0.
            enddo
         enddo
      enddo
      ivmax=0
      ivmin=1000
      do iy=1,200
         y=(float(iy)-100.5)*.01
         do ix=1,200
            x=(float(ix)-100.5)*.01
            r2=x**2+y**2
            r=sqrt(r2)
            phi=atan2(y,x)
            if(abs(r).le.r0.and.abs(r).ge.rmin)then
               vr=vr0*r/r0
               v=vexp*cos(phi)+vr*sin(phi)
               call hf1(12,v,1.)
               iv=ifix(v/dv+201.)
               if(iv.gt.ivmax)ivmax=iv
               if(iv.lt.ivmin)ivmin=iv
               flux(ix,iy,iv)=flux(ix,iy,iv)+rho
               call hf1(11,float(iv)-.5,1.)
               call hf1(9,v,rho)
               endif
            enddo
         enddo
         print*,'ivmin, ivmax: ',ivmin,ivmax

      do i=1,200
         do j=1,200
            xx=(float(i)-.5)*.01
            yy=(float(j)-.5)*.01
            do k=1,400
               call hf1(10,float(k)-200.5,flux(i,j,k))
               call hf2(20,xx,yy,flux(i,j,k))
               iik=k/80+1
               call hf2(20+iik,xx,yy,flux(i,j,k))
            enddo
         enddo
      enddo
      stop
      end
