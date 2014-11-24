ccc
ccc This is to re-sample data to a regular grid
ccc 

      program resample_data
      real co32(70,360,360),co65(43,432,432)
      real profmeas(2,2500,40),smeas(2,2500),fco(2,2500,40)  

      sum=0.

      do i=1,40
         do j=1,2500
            do k=1,2
               fco(k,j,i)=0.
            enddo
         enddo
      enddo

      do i=1,70
         do j=1,360
            do k=1,360
               co32(i,j,k)=0.
            enddo
         enddo
      enddo
      do i=1,43
         do j=1,432
            do k=1,432
               co65(i,j,k)=0.
            enddo
         enddo
      enddo
      co32(35,180,180)=1.
      co32(35,181,180)=1.
      co32(35,180,181)=1.
      co32(35,181,181)=1.

      sco32=0.
      sco65=0.
      
      do i=1,50
         do j=1,50
            sco32=0.
            sco65=0.
            do k=1,40           
         call transform(-181.5*.12,.12,i,-2.5,.1,floi,ilo,fhii,ihi)
         call transform(-177.5*.12,.12,j,-2.5,.1,floj,jlo,fhij,jhi)
         call transform(-36.5*.2115,.2115,k,-8.,.4,flok,klo,fhik,khi)
               do ii=ilo,ihi
                  wi=1.
                  if(ii.eq.ilo) wi=floi
                  if(ii.eq.ihi) wi=fhii
                  if(ilo.eq.ihi) wi=1-(1-floi+1-fhii)
                  do jj=jlo,jhi
                     wj=1.
                     if(jj.eq.jlo) wj=floj
                     if(jj.eq.jhi) wj=fhij
                     if(jlo.eq.jhi) wj=1-(1-floj+1-fhij)
                     do kk=klo,khi
                        wk=1.
                        if(kk.eq.klo) wk=flok
                        if(kk.eq.khi) wk=fhik
                        if(klo.eq.khi) wk=1-(1-flok+1-fhik)
                        fco(1,i+50*(j-1),k)=fco(1,i+50*(j-1),k)+
     &                       wi*wj*wk*co32(kk,jj,ii)

                     enddo
                  enddo
               enddo
          call transform(-220.5*.05,.05,i,-2.5,.1,floi,ilo,fhii,ihi)
          call transform(-210.5*.05,.05,j,-2.5,.1,floj,jlo,fhij,jhi)
          call transform(-20.7*.4,.4,k,-8.,.4,flok,klo,fhik,khi)

               do ii=ilo,ihi
                  wi=1.
                  if(ii.eq.ilo) wi=floi
                  if(ii.eq.ihi) wi=fhii
                  if(ilo.eq.ihi) wi=fhii-floi
                  do jj=jlo,jhi
                     wj=1.
                     if(jj.eq.jlo) wj=floj
                     if(jj.eq.jhi) wj=fhij
                     if(jlo.eq.jhi) wj=fhij-floj
                     do kk=klo,khi
                        wk=1.
                        if(kk.eq.klo) wk=flok
                        if(kk.eq.khi) wk=fhik
                        if(klo.eq.khi) wk=fhik-flok
                        fco(2,i+50*(j-1),k)=fco(2,i+50*(j-1),k)+
     &                       wi*wj*wk*co65(kk,jj,ii)                        
                     enddo
                  enddo
               enddo

               if(fco(1,i+50*(j-1),k).ne.0.)print*,i,j,k,
     &              fco(1,i+50*(j-1),k)
               sum=sum+fco(1,i+50*(j-1),k)
            enddo
             smeas(1,i+50*(j-1))=sco32
             smeas(2,i+50*(j-1))=sco65
          enddo
      enddo

      i0=120
      j0=21
c      print*,'',profmeas(1,i0,j0),profmeas(2,i0,j0)
c      print*,smeas(1,2000),smeas(2,2000)
      print*,'sum: ',sum
      stop
      end


      subroutine transform(x0,a0,i1,x1,a1,flo,ilo,fhi,ihi)
c transform an array starting at x0 in bins of a0 into an array starting at x1 in bins of a1
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

      if(flo.lt.0.or.flo.gt.1.or.fhi.lt.0.or.fhi.gt.1.)print*,flo,fhi
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
