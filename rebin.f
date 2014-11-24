      program rebin
c     
c the program rebins a three dimensional array a0(x(3)) into another a1(x(3))
c x(3) runs from x0(3) to x0(3)+n0(3)*d0(3) in a0 
c and from x1(3) to x1(3)+n1(3)*d1(3) in a1
c
      dimension x0(3),n0(3),d0(3),x1(3),n1(3),d1(3)
c      dimension ilo(3,50),ihi(3,50),flo(3,50),fhi(3,50)
      dimension ilo(3,100),ihi(3,100),flo(3,100),fhi(3,100)
ccc      dimension a0(360,360,70),a1(100,100,38)     ! for CO32
      dimension a0(432,432,43),a1(100,100,38)       ! for CO65
c$$$      data x0,n0,d0/-21.6,-21.6,-7.4,
c$$$     &     360,360,70,.12,.12,.2115/                    ! for CO32
      data x0,n0,d0/-10.8,-10.8,-8.6,
     &     432,432,43,.05,.05,.4/                       ! for CO65
ccc      data x1,n1,d1/-2.26,-2.71,-7.3,50,50,38,.1,.1,.4/   ! for CO32 
      data x1,n1,d1/-2.28,-2.76,-6.9,50,50,38,.1,.1,.4/   ! for CO65


      call HLIMAP(20000,'MAP')
      call hbook2(10,'',50.,0.,50.,50.,1.,50.,0.)
      call hbook1(11,'',50.,0.,50.,0.)
      call hbook1(12,'',50.,0.,50.,0.)
      call hbook1(13,'',38.,0.,38.,0.)     
c$$$      call hbook2(10,'',30.,-1.5,1.5,30.,-1.5,1.5,0.)
c$$$      call hbook1(11,'',30.,-1.5,1.5,0.)
c$$$      call hbook1(12,'',30.,-1.5,1.5,0.)
c$$$      call hbook1(13,'',20.,-4.,4.,0.)     

c     zero a1
c     
      do i=1,n1(1)
         do j=1,n1(2)
            do k=1,n1(3)
               a1(i,j,k)=0.
            enddo
         enddo
      enddo
c
c read a0
c 
c      open 2
c      read
c      close 2

      sum1=0.
ccc      open(2,file='12CO3-2_channel_image.txt')
      open(2,file='CO6-5_channel_image.txt')
      do k=1,n0(3)
         do j=1,n0(2)
            do i=1,n0(1)
               read(2,*) kk, jj, ii, aa
               if(kk.ne.k.or.jj.ne.j.or.ii.ne.i)then
                  print*,'read file problem: ',i,j,k,ii,jj,kk
                  stop
               endif
               a0(ii,jj,kk)=aa
c$$$               tmp1=abs(float(ii)-182.5)*.12
c$$$               tmp2=abs(float(jj)-178.5)*.12
c$$$               tmp3=abs(float(kk)-37.5)*.2115
c$$$c            if(tmp1.gt.1.8.or.tmp2.gt.1.8.or.tmp3.gt.3.4)then
c$$$            if(tmp1.gt.2.5.or.tmp2.gt.2.5.or.tmp3.gt.6.4)then
c$$$               a0(ii,jj,kk)=0.
c$$$            endif
c            sum1=sum1+a0(ii,jj,kk)               
            enddo
         enddo
      enddo
      close(2)
c      print*,'sum1: ',sum1
c
c     calculate for each bin i1(3) the bins ilo(3,i1) and ihi(3,i1) 
c     on which to loop in array a0 and the weights flo(3,i1), fhi(3,i1)
c     to give to bins ilo and ihi
c
      do 1 m=1,3
         do 2 i1=1,n1(m)
            ilo(m,i1)=1
            ihi(m,i1)=1
            flo(m,i1)=0.5
            fhi(m,i1)=0.5
            fi1=float(i1)
            xi1=x1(m)+(fi1-1.)*d1(m)
            y0m= x0(m)+float(n0(m))*d0(m)
            if(xi1.ge.y0m) go to 2
            yi1=xi1+d1(m)
            if(yi1.le.x0(m)) go to 2
            if(xi1.le.x0(m).and.yi1.le.x0(m)+d0(m)) then
               fhi(m,i1)=(yi1-x0(m))/d0(m)+0.5
               go to 2
            endif
            if(xi1.ge.y0m-d0(m).and.yi1.ge.y0m) then
               ilo(m,i1)=n0(m)
               ihi(m,i1)=n0(m)
               flo(m,i1)=(y0m-xi1)/d0(m)+0.5
               go to 2
            endif
            yy=(xi1-x0(m))/d0(m)
            iyy=ifix(yy)
            flo(m,i1)=float(iyy)+1.-yy
            ilo(m,i1)=iyy+1
            yy=yy+d1(m)/d0(m)
            iyy=ifix(yy)
            fhi(m,i1)=yy-float(iyy)
            ihi(m,i1)=iyy+1
 2       continue
 1    continue
c
c fill array a1 bin by bin
c
      do i=1,n1(1)
         do j=1,n1(2)
            do k=1,n1(3)
c
               do ii=ilo(1,i),ihi(1,i)
                  wi=1.
                  if(ii.eq.ilo(1,i)) wi=flo(1,i)
                  if(ii.eq.ihi(1,i)) wi=fhi(1,i)
                  if(ilo(1,i).eq.ihi(1,i)) wi=flo(1,i)+fhi(1,i)-1.
                  do jj=ilo(2,j),ihi(2,j)
                     wj=1.
                     if(jj.eq.ilo(2,j)) wj=flo(2,j)
                     if(jj.eq.ihi(2,j)) wj=fhi(2,j)
                     if(ilo(2,j).eq.ihi(2,j)) wj=flo(2,j)+fhi(2,j)-1.
                     do kk=ilo(3,k),ihi(3,k)
                        wk=1.
                        if(kk.eq.ilo(3,k)) wk=flo(3,k)
                        if(kk.eq.ihi(3,k)) wk=fhi(3,k)
                        if(ilo(3,k).eq.ihi(3,k)) wk=flo(3,k)+fhi(3,k)-1.
                        a1(i,j,k)=a1(i,j,k)+wi*wj*wk*a0(ii,jj,kk)  

c$$$               if(k.eq.17.and.i.eq.25.and.j.eq.25)print 99,i,j,k,
c$$$     &        ilo(1,i),ihi(1,i),ilo(2,j),ihi(2,j),ilo(3,k),ihi(3,k)
c$$$     &        ,flo(1,i),fhi(1,i),flo(2,j),fhi(2,j),flo(3,k),fhi(3,k)
c$$$ 99                     format (9i4,6(1x,e8.2))

                     enddo
                  enddo
               enddo
c               sum2=sum2+a1(i,j,k)
            enddo
         enddo
      enddo    


c      print*,'sum2: ',sum2,(sum1-sum2)/sum1*1000.
      print*,'check content 25, 25, 17:', a1(25,25,17)
c
c write a1
c
      print*,'here:'
      do ii=1,n1(1)
         do jj=1,n1(2)
            do kk=1,n1(3)
               write(26,*)ii,jj,kk,a1(ii,jj,kk)
               call hf2(10,float(ii)-.5,float(jj)-.5,a1(ii,jj,kk))
               call hf1(11,float(ii)-.5,a1(ii,jj,kk))
               call hf1(12,float(jj)-.5,a1(ii,jj,kk))              
               if(ii.gt.20.and.ii.lt.30.and.jj.gt.20.and.jj.lt.30)
     &              call hf1(13,float(kk)-.5,a1(ii,jj,kk))
            enddo
         enddo
      enddo

      stop 
      end
