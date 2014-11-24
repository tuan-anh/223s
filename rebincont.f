      program rebincont
c     
c the program rebins a three dimensional array a0(x(3)) into another a1(x(3))
c x(3) runs from x0(3) to x0(3)+n0(3)*d0(3) in a0 
c and from x1(3) to x1(3)+n1(3)*d1(3) in a1
c
      dimension x0(2),n0(2),d0(2),x1(2),n1(2),d1(2)
c      dimension ilo(3,50),ihi(3,50),flo(3,50),fhi(3,50)
      dimension ilo(2,100),ihi(2,100),flo(2,100),fhi(2,100)
      dimension a0(360,360),a1(100,100)   ! for CO32
ccc      dimension a0(432,432),a1(100,100)  ! for CO65
      data x0,n0,d0/-21.6,-21.6,
     &     360,360,.12,.12/                    ! for CO32
c$$$      data x0,n0,d0/-10.8,-10.8,
c$$$     &     432,432,.05,.05/                       ! for CO65
      data x1,n1,d1/-2.26,-2.71,50,50,.1,.1/   ! for CO32 
ccc      data x1,n1,d1/-2.28,-2.76,50,50,.1,.1/   ! for CO65


      call HLIMAP(20000,'MAP')
      call hbook2(10,'',50.,0.,50.,50.,1.,50.,0.)
      call hbook1(11,'',50.,0.,50.,0.)
      call hbook1(12,'',50.,0.,50.,0.)
      call hbook1(13,'',34.,0.,34.,0.)     
c$$$      call hbook2(10,'',30.,-1.5,1.5,30.,-1.5,1.5,0.)
c$$$      call hbook1(11,'',30.,-1.5,1.5,0.)
c$$$      call hbook1(12,'',30.,-1.5,1.5,0.)
c$$$      call hbook1(13,'',20.,-4.,4.,0.)     

      print*,'array: ',n1(1),n1(2),n0(1),n0(2)

c     zero a1
c     
      do i=1,n1(1)
         do j=1,n1(2)
            a1(i,j)=0.
         enddo
      enddo
c
c read a0
c 
      sum1=0.
      sum2=0.
      open(2,file='cont_co32_av_image.txt')
ccc      open(2,file='cont_co65_av_image.txt')
      do j=1,n0(2)
         do i=1,n0(1)
            read(2,*) jj, ii, aa
            if(jj.ne.j.or.ii.ne.i)then
               print*,'read file problem: ',i,j,ii,jj
               stop
            endif
            a0(ii,jj)=aa
c$$$               tmp1=abs(float(ii)-221.5)*.05
c$$$               tmp2=abs(float(jj)-211.5)*.05
c$$$            if(tmp1.gt.2..or.tmp2.gt.2.)then
c$$$               a0(ii,jj)=0.
c$$$            endif
c$$$            sum1=sum1+a0(ii,jj)
         enddo
      enddo
      close(2)
c
c     calculate for each bin i1(2) the bins ilo(2,i1) and ihi(2,i1) 
c     on which to loop in array a0 and the weights flo(2,i1), fhi(2,i1)
c     to give to bins ilo and ihi
c
      do 1 m=1,2
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
                  a1(i,j)=a1(i,j)+wi*wj*a0(ii,jj)  
               enddo
            enddo
            sum2=sum2+a1(i,j)
         enddo
      enddo    
c
c write a1
c
      print*,'here:',sum1,sum2,n1(1),n1(2)
      do ii=1,n1(1)
         do jj=1,n1(2)
            write(21,*)ii,jj,a1(ii,jj)    
         enddo
      enddo
c     
      stop 
      end
