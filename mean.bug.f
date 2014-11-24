      program red_rectangle
      real co32(70,360,360),co65(43,432,432),fco(2,2500,40)
c
      call HLIMAP(20000,'MAP')
	call hbook1(11,'spectrum CO(3-2)',40.,0.,40.,0.)
	call hbook1(12,'spectrum CO(6-5)',40.,0.,40.,0.)
	call hbook2(21,'map CO(3-2)',50.,0.,50.,50.,0.,50.,0.) 
	call hbook2(22,'map CO(6-5)',50.,0.,50.,50.,0.,50.,0.)
	call hbook2(30,'map CO(6-5)/CO(3-2)',50.,0.,50.,50.,0.,50.,0.) 
        call hbook1(23,'fco32',200.,-.2,1.2,0.)
        call hbook1(24,'fco65',200.,-.2,1.2,0.)
c
c read CO(3-2) data and store in co32
c
      open(2,file='12CO3-2_channel_image.txt')
      open(3,file='CO6-5_channel_image.txt')
      
	do k=1,70
         do j=1,360
            do i=1,360
               read(2,*) kk, jj, ii, aa
               co32(k,j,i) = aa
               x=float(ii)-182.5
               y=float(jj)-178.5
               v = (float(k)-37.5)*.2115
            enddo
         enddo
      enddo
c
c read co65 data and store in co65
c
      do k=1,43
         do j=1,432
            do i=1,432
               read(3,*) kk, jj, ii, aa
               co65(k,j,i)= aa
               x=float(ii)-221.5
               y=float(jj)-211.5
               v = (float(k)-22.7)*.4
            enddo
         enddo
      enddo
c
c convert velocity spectra in 40 bins of 0.4 km/s from -8 km/s to 8 km/s
c map pixels in 50 bins of .1 starting at -2.5 arcsec
c

      tmp =0.

	do i=1,50
	    do j=1,50
		sco32=0.
		sco65=0.
            do k=1,40
		call transform(-181.5*.12,.12,i,-2.5,.1,floi,ilo)
		call transform(-177.5*.12,.12,j,-2.5,.1,floj,jlo)
		call transform(-36.5*.2115,.2115,k,-8.,.4,flok,klo)
		fco(1,i+50*(j-1),k)=
     & floi*floj*(flok*co32(klo,jlo,ilo)+(1.-flok)*co32(klo+1,jlo,ilo))+
     & floi*(1.-floj)*(flok*co32(klo,jlo+1,ilo)+(1.-flok)
     & *co32(klo+1,jlo+1,ilo))+
     & (1.-floi)*floj*(flok*co32(klo,jlo,ilo+1)+(1.-flok)
     &             *co32(klo+1,jlo,ilo+1))+
     & (1.-floi)*(1.-floj)*(flok*co32(klo,jlo+1,ilo+1)+(1.-flok)
     &               *co32(klo+1,jlo+1,ilo+1))
		call transform(-220.5*.05,.05,i,-2.5,.1,floi,ilo)
		call transform(-210.5*.05,.05,j,-2.5,.1,floj,jlo)
		call transform(-21.7*.4,.4,k,-8.,.4,flok,klo)

		fco(2,i+50*(j-1),k)=
     & floi*floj*(flok*co65(klo,jlo,ilo)+(1.-flok)*co65(klo+1,jlo,ilo))+
     & floi*(1.-floj)*(flok*co65(klo,jlo+1,ilo)+(1.-flok)
     &               *co65(klo+1,jlo+1,ilo))+
     & (1.-floi)*floj*(flok*co65(klo,jlo,ilo+1)+(1.-flok)
     &               *co65(klo+1,jlo,ilo+1))+
     & (1.-floi)*(1.-floj)*(flok*co65(klo,jlo+1,ilo+1)+(1.-flok)
     &               *co65(klo+1,jlo+1,ilo+1))
	call hf1(11,float(k)-.5,fco(1,i+50*(j-1),k))
	call hf1(12,float(k)-.5,fco(2,i+50*(j-1),k))
	call hf2(21,float(i)-.5,float(j)-.5,fco(1,i+50*(j-1),k))
	call hf2(22,float(i)-.5,float(j)-.5,fco(2,i+50*(j-1),k))
	sco32=sco32+fco(1,i+50*(j-1),k)
	sco65=sco65+fco(2,i+50*(j-1),k)
        call hf1(23,fco(1,i+50*(j-1),k),1.)
        call hf1(24,fco(2,i+50*(j-1),k),1.)
      enddo
             sco321=.2
             if(sco32.gt..2)sco321=sco32
c             call hf2(30,float(i)-.5,float(j)-.5,sco65/sco321)
             tmp2=sco65/sco321
             if(tmp2.lt.0.)print*, i,j, sco65
             tmp=tmp+sco65/sco321

             if(i.gt.4.and.i.lt.45.and.j.gt.4.and.j.lt.45)then
                if(tmp2.lt.0.)print*, i,j, sco65
                tmp=tmp+sco65/sco321
                call hf2(30,float(i)-.5,float(j)-.5,sco65/sco32)
             endif
        
          enddo
      enddo
      print *,' ', tmp
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
	return
	end


















