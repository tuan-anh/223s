      program colines
      dimension co32(50,50,38),co65(50,50,38),s32(50,50),s65(50,50)
     &     ,t32(38),t65(38),cont32(50,50),cont65(50,50)
      logical center, between
      call HLIMAP(20000,'COMP')
      call hbook2(201,'',50.,0.,50.,50.,0.,50.,0.)
      call hbook2(202,'',50.,0.,50.,50.,0.,50.,0.)
      call hbook1(203,'',50.,0.,50.,0.)
      call hbook1(204,'',50.,0.,50.,0.)
      call hbook1(205,'',50.,0.,50.,0.)
      call hbook1(206,'',50.,0.,50.,0.)
      call hbook1(301,'',38.,0.,38.,0.)     
      call hbook1(302,'',38.,0.,38.,0.)    
      call hbook2(305,'',50.,0.,50.,50.,0.,50.,0.)
      
      do ii=1,50
         call hbook1(400+ii,'',38.,0.,38.,0.)
      enddo

c$$$      call hbook2(101,'',50.,0.,50.,50.,0.,50.,0.)
c$$$      call hbook2(102,'',50.,0.,50.,50.,0.,50.,0.)
      call hbook1(107,'',100.,-.003,.003,0.)
      call hbook1(108,'',100.,-.2,.2,0.)
c$$$
      call hbook1(303,'',100.,-2.,4.,0.)
      call hbook1(304,'',100.,-10.,20.,0.)



      do ii=1,4
         call hbook1(102+ii,'',50.,0.,50.,0.)
      enddo

      call hbook1(207,'',100.,-.2,.5,0.)
      call hbook1(208,'',100.,-2.,3.,0.)

c
c     read continuum data
c
      open(2,file='co32-rebin-cont2.txt')
      open(3,file='co65-rebin-cont2.txt')
      do ii=1,50
         do jj=1,50
            read(2,*)i1,j1,aa
            read(3,*)i2,j2,bb
            if(ii.ne.i1.or.ii.ne.i2.or.jj.ne.j1.or.jj.ne.j2)then
               print*,'read file problem: ',ii,jj,i1,j1,i2,j2
               stop
            endif
            cont32(ii,jj)=aa
            cont65(ii,jj)=bb
            call hf2(101,float(ii)-.5,float(jj)-.5,aa)
            call hf2(102,float(ii)-.5,float(jj)-.5,bb)
            call hf1(103,float(ii)-.5,cont32(ii,jj))            
            call hf1(104,float(jj)-.5,cont32(ii,jj))
            call hf1(105,float(ii)-.5,cont65(ii,jj))            
            call hf1(106,float(jj)-.5,cont65(ii,jj))
            call hf1(107,cont32(ii,jj),1.)
            call hf1(108,cont65(ii,jj),1.)
         enddo
      enddo
      close(2)
      close(3)

c
c     read lines data
c
      open(2,file='co32-rebin2.txt')
      open(3,file='co65-rebin2.txt')

      do ii=1,50
         do jj=1,50
            do kk=1,38
               read(2,*) i1,j1,k1,aa
               daa=.5*(1.+.032*(float(ii)-25.5))
     &              +.02*(1.+.28*(float(jj)-25.5))
               aa=aa+daa/38.
               read(3,*) i2,j2,k2,bb
               dbb=2.*(1.+.04*(float(ii)-25.5))
     &              +1.6*(1.+.02*(float(jj)-25.5))
               bb=bb+dbb/38.
               if(i1.ne.i2.or.i1.ne.ii.or.j1.ne.j2.or.j1.ne.jj
     &              .or.k1.ne.k2.or.k1.ne.kk) then
                  print*,'read problem: ',ii,jj,kk,i1,j1,k1,i2,j2,k2
                  stop
               endif
               call hf1(207,aa,1.)
               call hf1(208,bb,1.)        
               if(aa.lt..052)aa=0.
               if(bb.lt..616)bb=0.
               co32(ii,jj,39-kk)=aa
               co65(ii,jj,kk)=bb
                                       
               call hf2(201,float(ii)-.5,float(jj)-.5,aa)
               call hf2(202,float(ii)-.5,float(jj)-.5,bb)
               call hf1(203,float(ii)-.5,aa)
               call hf1(204,float(jj)-.5,aa)
               call hf1(205,float(ii)-.5,bb)
               call hf1(206,float(jj)-.5,bb)
         
               if(ii.ge.20.and.ii.le.30.and.jj.ge.20.and.jj.le.30)then
                  call hf1(301,float(39-kk)-.5,aa)
                  call hf1(302,float(kk)-.5,bb)
               endif
            enddo
         enddo
      enddo
      close(2)
      close(3)
c
ccc     look spectrum by spectrum
c     
      do ii=1,50
         do jj=1,50
            s32(ii,jj)=0.
            s65(ii,jj)=0.
            do kk=1,38
               write(27,*)ii,jj,kk,co32(ii,jj,kk)
               write(28,*)ii,jj,kk,co65(ii,jj,kk)
               s32(ii,jj)=s32(ii,jj)+co32(ii,jj,kk)
               s65(ii,jj)=s65(ii,jj)+co65(ii,jj,kk)
            enddo            
            call hf1(303,s32(ii,jj),1.)
            call hf1(304,s65(ii,jj),1.)                    
            if(s32(ii,jj).gt..84*2..and.s65(ii,jj).gt.4.4*2.)then
               ratio=s65(ii,jj)/s32(ii,jj)
               call hf2(305,float(ii)-.5,float(jj)-.5,ratio)               
            endif
         enddo
      enddo

      center=.false.
      between=.true.
      do kch=1,25
         mm=(kch-1)/5
         ii1=10*mm+1
         if(center)ii1=2*mm+21
         if(between)ii1=6*mm+11
         ii2=ii1+9
         if(center)ii2=ii1+1
         if(between)ii2=ii1+5
         nn=kch-5*mm-1
         jj1=10*nn+1
         if(center)jj1=2*nn+21
         if(between)jj1=6*nn+11
         jj2=jj1+9               
         if(center)jj2=jj1+1
         if(between)jj2=jj1+5
         print*,kch,jj1,jj2,ii1,ii2
      do kk=1,38
         t32(kk)=0.
         t65(kk)=0.
         do ii=ii1,ii2
            do jj=jj1,jj2
               t32(kk)=t32(kk)+co32(jj,ii,kk)
               t65(kk)=t65(kk)+co65(jj,ii,kk)               
            enddo
         enddo
         call hf1(400+kch,float(kk)-.5,t32(kk))
         call hf1(425+kch,float(kk)-.5,t65(kk))
      enddo
      enddo
   
      stop
      end
