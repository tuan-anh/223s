      program fluxcomp
c      dimension s1(85000),s2(85000)
      call HLIMAP(20000,'COMP')
      call hbook1(10,'',1000.,-1.e-4,1.e-4,0.)
      open(2,file='fort.13')
      open(3,file='fort.14')

      do ii=1,50
         do jj=1,50
            do kk=1,34
               read(2,*) i1,j1,k1,aa
               read(3,*) i2,j2,k2,bb
               if(i1.ne.i2.or.i1.ne.ii.or.j1.ne.j2.or.j1.ne.jj
     & .or.k1.ne.k2.or.k1.ne.kk) then
                  print*,ii,jj,kk,i1,j1,k1,i2,j2,k2
                  stop
               endif
               tmp=aa-bb
               call hf1(10,tmp,1.)
            enddo
         enddo
      enddo

      stop
      end
