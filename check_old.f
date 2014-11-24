      program check
      real co32(70,360,360),co65(43,432,432)
      real profmeas(2,2500,40),smeas(2,2500),fco(2,2500,40)  

     
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
c$$$      co32(35,181,180)=1.
c$$$      co32(35,180,181)=1.
c$$$      co32(35,181,181)=1.

      sco32=0.
      sco65=0.
      
      do i=1,50
         do j=1,50
            sco32=0.
            sco65=0.
            do k=1,40           
         call transform(-181.5*.12,.12,i,-2.5,.1,floi,ilo,fhii,ihi,ni)
         call transform(-177.5*.12,.12,j,-2.5,.1,floj,jlo,fhij,jhi,nj)
         call transform(-36.5*.2115,.2115,k,-8.,.4,flok,klo,fhik,khi,nk)
c            print*,floi,ilo,fhii,ihi,floj,jlo,fhij,jhi,flok,klo,fhik,khi
               do ii=ilo-1,ihi+1,ni
                  wi=1.
                  if(ii.eq.ilo-1) wi=floi
                  if(ii.eq.ihi+1) wi=fhii
                  do jj=jlo-1,jhi+1,nj
                     wj=1.
                     if(jj.eq.jlo-1) wj=floj
                     if(jj.eq.jhi+1) wj=fhij
                     do kk=klo-1,khi+1,nk
                        wk=1.
                        if(kk.eq.klo-1) wk=flok
                        if(kk.eq.khi+1) wk=fhik
                        fco(1,i+50*(j-1),k)=fco(1,i+50*(j-1),k)+
     &                       wi*wj*wk*co32(kk,jj,ii)

                     enddo
                  enddo
               enddo
          call transform(-220.5*.05,.05,i,-2.5,.1,floi,ilo,fhii,ihi,ni)
          call transform(-210.5*.05,.05,j,-2.5,.1,floj,jlo,fhij,jhi,nj)
          call transform(-20.7*.4,.4,k,-8.,.4,flok,klo,fhik,khi,nk)
               do ii=ilo-1,ihi+1,ni
                  wi=1.
                  if(ii.eq.ilo-1) wi=floi
                  if(ii.eq.ihi+1) wi=fhii
                  do jj=jlo-1,jhi+1,nj
                     wj=1.
                     if(jj.eq.jlo-1) wj=floj
                     if(jj.eq.jhi+1) wj=fhij
                     do kk=klo-1,khi+1,nk
                        wk=1.
                        if(kk.eq.klo-1) wk=flok
                        if(kk.eq.khi+1) wk=fhik

                        fco(2,i+50*(j-1),k)=fco(2,i+50*(j-1),k)+
     &                       wi*wj*wk*co65(kk,jj,ii)                        
                     enddo
                  enddo
               enddo

c$$$               if(fco(1,i+50*(j-1),k).ne.0.)print*,i,j,k,
c$$$     &              fco(1,i+50*(j-1),k)

            enddo
             smeas(1,i+50*(j-1))=sco32
             smeas(2,i+50*(j-1))=sco65
c             print*,i,j,i+50*(j-1),sco32,sco65
          enddo
      enddo

      i0=120
      j0=21
c      print*,'',profmeas(1,i0,j0),profmeas(2,i0,j0)
c      print*,smeas(1,2000),smeas(2,2000)
      stop
      end




      subroutine transform(x0,a0,i1,x1,a1,flo,ilo,fhi,ihi,nstep)
c transform an array starting at x0 in bins of a0 into an array starting at x1 in bins of a1
c under the assumption that a1>=a0. Namely each new bin contains at least 1 old binl
c bin i1 will receive a content that is flo*the content of bin ilo-1, +fhi*the content of bin ihi+1 
c +the contents of bins ilo to ihi
      xlo=x1+(float(i1)-1.)*a1
      xhi=xlo+a1
      ylo=(xlo-x0)/a0
c      ilo=ifix1(ylo)+2
      ilo=ifix1(ylo)
      yhi=(xhi-x0)/a0
c      ihi=ifix1(yhi)
      ihi=ifix1(yhi)
c      flo=(float(ilo-1)-xlo)/a0
c      fhi=(xhi-float(ihi))/a0
      flo=(xlo-float(ilo))/a0
      fhi=(xhi-float(ihi))/a0
      nstep=1
      if(ilo.eq.ihi+1)then
         flo1=fhi
         fhi=flo
         flo=flo1
         nstep=-1
      endif
      if(ilo.eq.ihi+2)then
         ilo=ihi+1
         ihi=ilo
         flo=(flo+fhi-1.)/2.
         fhi=flo
      endif
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
