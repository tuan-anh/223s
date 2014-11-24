      program rerange
c
c this version is to resample a regular set of data
c aa is sampled by one-dimension grid size 1.
c bb is a new sample with grid size 1.2
c
      real aa(30), bb(20)

      call HLIMAP(20000,'MAP')
      pi=2.*acos(0.)

      do ii=1,30
         aa(ii)=3.
      enddo

      do i=1,20
c         call transform(-181.5*.12,.12,i,-2.5,.1,floi,ilo)
         call transform(0.,1.,i,0.,1.5,floi,ilo)
         bb(i)=floi*aa(ilo)+(1.-floi)*aa(ilo+1)
         print*,i,bb(i),aa(i), ilo
      enddo

      stop 
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
cc      print*,'inew: ',inew, ilo, xnew
      return
      end
