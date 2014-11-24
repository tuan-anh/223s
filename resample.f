      program resample
      real co32(70,360,360),co65(43,432,432)
      common/meas/profmeas(2,2500,40),smeas(2,2500),fco(2,2500,34)  
      dimension g(70,360,360)

      call HLIMAP(20000,'MAP')

      call hbook1(31,'',60.,-30.,30.,0.)
      call hbook1(32,'',60.,-30.,30.,0.)
      call hbook1(33,'',120.,-60.,60.,0.)
      call hbook1(41,'',60.,-30.,30.,0.)
      call hbook1(42,'',60.,-30.,30.,0.)
      call hbook1(43,'',120.,-60.,60.,0.)
      
      call hbook2(50,'',30.,170.,200.,30.,170.,200.,0.)
      call hbook1(51,'',70.,0.,70.,0.)
      
     
c
c read CO(3-2) data and store in co32
c ii, jj and kk are simply the numbers in the datacubes, incremented by 1
c
      sum3=0.
      sum4=0.

      open(2,file='12CO3-2_channel_image.txt')
      open(3,file='CO6-5_channel_image.txt')
      
      do k=1,70
         do j=1,360
            do i=1,360
               read(2,*) kk, jj, ii, aa
               if(kk.ne.k.or.jj.ne.j.or.ii.ne.i)then
                  print*,i,j,k,ii,jj,kk
                  stop
               endif
               co32(kk,jj,ii)=aa
               g(kk,jj,ii)=1.
               tmp1=abs(float(ii)-182.5)*.12
               tmp2=abs(float(jj)-178.5)*.12
               tmp3=abs(float(kk)-37.5)*.2115
c$$$            if(tmp1.gt.1.8.or.tmp2.gt.1.8.or.tmp3.gt.3.4)then
c$$$               co32(kk,jj,ii)=0.
c$$$               g(kk,jj,ii)=0.
c$$$            endif
            call hf1(31,float(i)-182.5,co32(kk,jj,ii))
            call hf1(32,float(j)-178.5,co32(kk,jj,ii))
            call hf1(33,float(k)-37.5,co32(kk,jj,ii))
            sum3=sum3+co32(kk,jj,ii)               
            enddo
         enddo
      enddo
      close(2)

c      print*,'sum3: ',sum3
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
               if(kk.ne.k.or.jj.ne.j.or.ii.ne.i)then
                  print*,i,j,k,ii,jj,kk
                  stop
               endif
c               co65(44-kk,jj,ii)= aa
               co65(44-kk,jj,ii)= 1.
               tmp1=abs(float(ii)-221.5)*.05
               tmp2=abs(float(jj)-211.5)*.05
               tmp3=abs(float(kk)-21.7)*.4
            if(tmp1.gt.2..or.tmp2.gt.2..or.tmp3.gt.5.4)then
               co65(44-kk,jj,ii)=0.
            endif
            sum4=sum4+co65(44-kk,jj,ii)
            enddo
         enddo
      enddo
      print*,'sum4: ',sum4
      close(3)
c
c convert velocity spectra in 40 bins of 0.4 km/s from -8 km/s to 8 km/s
c map pixels in 50 bins of .1 starting at -2.5 arcsec
c
      sco32=0.
      sco65=0.
      sum1=0.
      sum2=0.
      
      do i=1,50
         do j=1,50
            do k=1,34
               fco(1,i+50*(j-1),k)=0.
               fco(2,i+50*(j-1),k)=0.
            enddo
         enddo
      enddo

      do i=1,50
         do j=1,50
            sco32=0.
            sco65=0.
            do k=1,34
       call transform(-181.5*.12,.12,360,i,-2.5,.1,floi,ilo,fhii,ihi)
       call transform(-177.5*.12,.12,360,j,-2.5,.1,floj,jlo,fhij,jhi)
       call transform(-36.5*.2115,.2115,70,k,-6.8,.4,flok,klo,fhik,khi)
         if(ilo.gt.360.or.jlo.gt.360.or.klo.gt.70) goto 111
         if(ihi.lt.1.or.jhi.lt.1.or.khi.lt.1) goto 111
               do ii=ilo,ihi
                  wi=1.
                  if(ii.eq.ilo) wi=floi
                  if(ii.eq.ihi) wi=fhii
                  if(ilo.eq.ihi) wi=1.-(1.-floi+1.-fhii)
                  do jj=jlo,jhi
                     wj=1.
                     if(jj.eq.jlo) wj=floj
                     if(jj.eq.jhi) wj=fhij
                     if(jlo.eq.jhi) wj=1.-(1.-floj+1.-fhij)
                     do kk=klo,khi                       
                        wk=1.
                        if(kk.eq.klo) wk=flok
                        if(kk.eq.khi) wk=fhik
                        if(klo.eq.khi) wk=1.-(1.-flok+1.-fhik)
                        fco(1,i+50*(j-1),k)=fco(1,i+50*(j-1),k)+
     &                       wi*wj*wk*co32(kk,jj,ii)   
                        if(kk.eq.45.and.klo.ne.43.and.k.ne.21)
     & print 99,i,j,k,
     & kk,ilo,ihi,jlo,jhi,klo,khi,floi,floj,flok,fhii,fhij,fhik
                        g(kk,jj,ii)=g(kk,jj,ii)-wi*wj*wk*co32(kk,jj,ii)         
                     enddo
                  enddo
               enddo
 111           continue
 99            format (10i4,4(1x,e8.2))
               call hf1(41,float(i)-25.5,fco(1,i+50*(j-1),k))
               call hf1(42,float(j)-25.5,fco(1,i+50*(j-1),k))
               call hf1(43,float(k)-17.5,fco(1,i+50*(j-1),k))
               sum1=sum1+fco(1,i+50*(j-1),k)

       call transform(-220.5*.05,.05,432,i,-2.5,.1,floi,ilo,fhii,ihi)
       call transform(-210.5*.05,.05,432,j,-2.5,.1,floj,jlo,fhij,jhi)
       call transform(-20.7*.4,.4,43,k,-8.,.4,flok,klo,fhik,khi)
          if(ilo.gt.432.or.jlo.gt.432.or.klo.gt.43)goto 222
          if(ihi.lt.1.or.jhi.lt.1.or.khi.lt.1) goto 222
               do ii=ilo,ihi
                  wi=1.
                  if(ii.eq.ilo) wi=floi
                  if(ii.eq.ihi) wi=fhii
                  if(ilo.eq.ihi) wi=1.-(1.-floi+1.-fhii)
                  do jj=jlo,jhi
                     wj=1.
                     if(jj.eq.jlo) wj=floj
                     if(jj.eq.jhi) wj=fhij
                     if(jlo.eq.jhi) wj=1.-(1.-floj+1.-fhij)
                     do kk=klo,khi
                        wk=1.
                        if(kk.eq.klo) wk=flok
                        if(kk.eq.khi) wk=fhik
                        if(klo.eq.khi) wk=1.-(1.-flok+1.-fhik)
                        fco(2,i+50*(j-1),k)=fco(2,i+50*(j-1),k)+
     &                       wi*wj*wk*co65(kk,jj,ii)                        
                     enddo
                  enddo
               enddo
 222           continue
               write(11,*)i+50*(j-1),k,fco(1,i+50*(j-1),k)
               write(12,*)i+50*(j-1),k,fco(2,i+50*(j-1),k)
               sco32=sco32+fco(1,i+50*(j-1),k)
               sco65=sco65+fco(2,i+50*(j-1),k)
              
               write(14,*),i,j,k,fco(1,i+50*(j-1),k)


               sum2=sum2+fco(2,i+50*(j-1),k)
            enddo
             smeas(1,i+50*(j-1))=sco32
             smeas(2,i+50*(j-1))=sco65
          enddo
      enddo

      print*,'sum: ',sum1,sum3,sum2,sum4
      print*,'content 25,25,17: ',fco(1,1225,17)

      do ii=1,360
         do jj=1,360
            do kk=1,70
c               if(g(kk,jj,ii).gt.1.e-6.and.kk.ne.45)then
c                  call hf2(50,float(jj)-.5,float(ii)-.5,1.)
c                  call hf1(51,float(kk)-.5,1.)
c                  print*,ii,jj,kk,g(kk,jj,ii)
c               endif
            enddo
         enddo
      enddo



      call plot
      i0=120
      j0=21
c      print*,'',profmeas(1,i0,j0),profmeas(2,i0,j0)
c      print*,smeas(1,2000),smeas(2,2000)
      stop
      end


      subroutine transform(x0,a0,max,i1,x1,a1,flo,ilo,fhi,ihi)
c transform an array starting at x0 in bins of a0 into an array starting at x1 in bins of a1
c under the assumption that a1>=a0. Namely each new bin contains at least 1 old binl
c bin i1 will receive a content that is flo*the content of bin ilo-1, +fhi*the content of bin ihi+1 
c +the contents of bins ilo to ihi
      xlo=x1+(float(i1)-1.)*a1
      xhi=xlo+a1
      ylo=(xlo-x0)/a0
      ilo=ifix1(ylo)
      yhi=(xhi-x0)/a0
      ihi=ifix1(yhi)
      flo=(float(ilo+1)-ylo)
      fhi=(yhi-float(ihi))
      if(ylo.lt.0.)then
         flo=1.
         ilo=0         
      endif
      if(yhi.ge.float(max))then
         fhi=1.
         ihi=max-1
      endif
      if(flo.lt.0.or.flo.gt.1.or.fhi.lt.0.or.fhi.gt.1.)then
         print*,flo,fhi
         print*,x0,a0,max,x1,a1
      endif
      ilo=ilo+1
      ihi=ihi+1
      return
      end

      integer function ifix1(x)
      if(float(ifix(x)).eq.x) then
         ifix1=ifix(x)
         return
      endif
      if(x.ge.0.)ifix1=ifix(x)
      if(x.lt.0.)ifix1=ifix(x)-1
      return
      end

      subroutine plot
      common/meas/profmeas(2,2500,40),smeas(2,2500),fco(2,2500,34)  
      call hbook1(11,'spectrum CO(3-2)meas',40.,0.,40.,0.)
      call hbook1(12,'spectrum CO(6-5)meas',40.,0.,40.,0.)
      call hbook2(21,'map CO(3-2)meas',50.,0.,50.,50.,0.,50.,0.) 
      call hbook2(22,'map CO(6-5)meas',50.,0.,50.,50.,0.,50.,0.)
      do ico=1,2
         do ij=1,2500
            j=ij/50+1
            i=ij-50*(j-1)    
            call hf2(20+ico,float(i)-.5,float(j)-.5,smeas(ico,ij))
            do k=1,40
               call hf1(10+ico,float(k)-.5,fco(ico,ij,k))
            enddo
         enddo
      enddo
      return
      end


   
      subroutine dtransform(dx0,da0,max,i1,dx1,da1,dflo,ilo,dfhi,ihi)
c transform an array starting at x0 in bins of a0 into an array starting at x1 in bins of a1
c under the assumption that a1>=a0. Namely each new bin contains at least 1 old binl
c bin i1 will receive a content that is flo*the content of bin ilo-1, +fhi*the content of bin ihi+1 
c +the contents of bins ilo to ihi
      double precision x1,a1,xlo,xhi,ylo,yhi,x0,a0,flo,fhi
      x0=dx0
      a0=da0
      x1=dx1
      a1=da1
      xlo=x1+(dfloat(i1)-1.)*a1
      xhi=xlo+a1
      ylo=(xlo-x0)/a0
      ilo=ifixd1(ylo)
      yhi=(xhi-x0)/a0
      ihi=ifixd1(yhi)
      flo=(dfloat(ilo+1)-ylo)
      fhi=(yhi-dfloat(ihi))
      if(ylo.lt.0.)then
         flo=1.
         ilo=0         
      endif
      if(yhi.ge.float(max))then
         fhi=1.
         ihi=max-1
      endif
      if(flo.lt.0.)flo=0.
      if(flo.gt.1.)flo=1.
      if(fhi.lt.0.)fhi=0.
      if(fhi.gt.1.)fhi=1.
      ilo=ilo+1
      ihi=ihi+1
      dfhi=fhi
      dflo=flo
      return
      end
      integer function ifixd1(dx)
      double precision dx     
      x=dx
      if(float(ifix(x)).eq.x) then
         ifixd1=ifix(x)
         return
      endif
      if(x.ge.0.)ifixd1=ifix(x)
      if(x.lt.0.)ifixd1=ifix(x)-1
      return
      end
