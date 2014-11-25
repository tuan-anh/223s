      program asym_co
      dimension profmeas(2,2500,38),smeas(2,2500),ameas(2,36)
      dimension rmeas(2,8,25), bmeas(2,2,25)
      dimension alim(8)
      data alim,rlim/35.,75.,105.,145.,220.,260.,280.,320.,15./  
      call HLIMAP(20000,'MAP')
      pi=2.*acos(0.)
      open(2,file='co32-rebin3.txt')
      open(3,file='co65-rebin3.txt')
      do m=1,36
         ameas(1,m)=0.
         ameas(2,m)=0.
      enddo
      do n=1,8
         call hbook1(110+n,'rmeas1',25.,0.,25.,0.)
         call hbook1(120+n,'rmeas2',25.,0.,25.,0.)
         call hbook1(130+n,'rmeas2/1',25.,0.,25.,0.)
         do m=1,25
            rmeas(1,n,m)=0.
            rmeas(2,n,m)=0.
            bmeas(1,1,m)=0.
            bmeas(1,2,m)=0.
            bmeas(2,1,m)=0.
            bmeas(2,2,m)=0.
         enddo
      enddo
      call hbook1(101,'ameas1',36.,0.,36.,0.)
      call hbook1(102,'ameas2',36.,0.,36.,0.)
      call hbook1(103,'ameas2/1',36.,0.,36.,0.)
      call hbook2(21,'map1',20.,15.,35.,20.,15.,35.,0.)
      call hbook2(22,'map2',20.,15.,35.,20.,15.,35.,0.)
      call hbook2(23,'map2/1',20.,15.,35.,20.,15.,35.,0.)
      call hbook1(31,'',100.,0.,1000.,0.)
      call hbook1(32,'',100.,0.,2000.,0.)
      call hbook1(141,'',25.,0.,25.,0.)
      call hbook1(142,'',25.,0.,25.,0.)
      call hbook1(151,'',25.,0.,25.,0.)
      call hbook1(152,'',25.,0.,25.,0.)
      call hbook1(161,'',25.,0.,25.,0.)
      call hbook1(162,'',25.,0.,25.,0.)
      do i=1,50
         do j=1,50
            ij=i+50*(j-1)
            xij=float(i)-25.5
            yij=float(j)-25.5
            rij=sqrt(xij**2+yij**2)
            nrij=ifix(rij)+1
            aij=atan2(yij,xij)*180./pi
            aij=aij-13.
            if(aij.lt.0.) aij=aij+360.
            lij=0
            if(aij.lt.45..or.aij.gt.315.) lij=1
            if(aij.gt.135..and.aij.lt.225.) lij=1
            if(aij.gt.40..and.aij.lt.70.) lij=2
            if(aij.gt.110..and.aij.lt.150.) lij=2
            if(aij.gt.220..and.aij.lt.260.) lij=2
            if(aij.gt.280..and.aij.lt.320.) lij=2
            nij=1
            do na=1,7
               if(aij.ge.alim(na).and.aij.lt.alim(na+1)) nij=na+1
            enddo
            mij=0
            if(rij.le.rlim)mij=ifix(aij/10.)+1
            smeas(1,ij)=0.
            smeas(2,ij)=0.
            if(i.eq.3.and.j.eq.48)print*,'north-east',aij,mij
            if(i.eq.3.and.j.eq.3)print*,'south-east',aij,mij
            if(i.eq.48.and.j.eq.48)print*,'north-west',aij,mij
            if(i.eq.48.and.j.eq.3)print*,'south-west',aij,mij
            do k=1,38              
               read(2,*) i1,j1,k1,aa
               read(3,*) i2,j2,k2,bb
               if(i.ne.i1.or.i.ne.i2.or.j.ne.j1.or.j.ne.j2
     &              .or.k.ne.k1.or.k.ne.k2)then
                  print*,'read file problem: ',i,j,k,i1,j1,k1,i2,j2,k2
                  stop
               endif
               profmeas(1,ij,k)=aa
               profmeas(2,ij,k)=bb
               smeas(1,ij)=smeas(1,ij)+aa
               smeas(2,ij)=smeas(2,ij)+bb
               if(mij.ne.0) ameas(1,mij)=ameas(1,mij)+aa
               if(mij.ne.0) ameas(2,mij)=ameas(2,mij)+bb


c$$$  if(ameas(1,mij).gt.300..or.ameas(2,mij).gt.1000.)
c$$$  &             print*,aa,bb,mij,ij,i,j,ameas(1,mij),ameas(2,mij)
               if(nrij.le.25) then
                  rmeas(1,nij,nrij)=rmeas(1,nij,nrij)+aa 
                  rmeas(2,nij,nrij)=rmeas(2,nij,nrij)+bb
               endif
               if(lij.ne.0)then
                  bmeas(1,lij,nrij)=bmeas(1,lij,nrij)+aa
                  bmeas(2,lij,nrij)=bmeas(2,lij,nrij)+bb
               endif                  
            enddo
		if(i.gt.15.and.i.le.35..and.j.gt.15.and.j.le.35) then
                   call hf2(21,float(i)-.5,float(j)-.5,smeas(1,ij))
                   call hf2(22,float(i)-.5,float(j)-.5,smeas(2,ij))
                   ratio=0.
                   if(smeas(1,ij).gt.0.) ratio=smeas(2,ij)/smeas(1,ij)
                   call hf2(23,float(i)-.5,float(j)-.5,ratio)
		endif               
             enddo
          enddo
          do mij=1,36
             call hf1(31,ameas(1,mij),1.)
             call hf1(32,ameas(2,mij),1.)
             call hf1(101,float(mij)-.5,ameas(1,mij))
             call hf1(102,float(mij)-.5,ameas(2,mij))
             ratio=0.
             if(ameas(1,mij).gt.0.) ratio=ameas(2,mij)/ameas(1,mij)
             call hf1(103,float(mij)-.5,ratio)
          enddo
          do nij=1,8
             do nrij=1,25
                call hf1(110+nij,float(nrij)-.5,rmeas(1,nij,nrij))
                call hf1(120+nij,float(nrij)-.5,rmeas(2,nij,nrij))                
                ratio=0.
                if(rmeas(2,nij,nrij).gt.0.) 
     &               ratio=rmeas(2,nij,nrij)/rmeas(1,nij,nrij)
                call hf1(130+nij,float(nrij)-.5,ratio)
             enddo
          enddo
          do lij=1,2
             do nrij=1,25
                call hf1(140+lij,float(nrij)-.5,bmeas(1,lij,nrij))
                call hf1(150+lij,float(nrij)-.5,bmeas(2,lij,nrij))                
                ratio=0.
                if(rmeas(1,nij,nrij).gt.0.) 
     &               ratio=rmeas(2,nij,nrij)/rmeas(1,nij,nrij)
                ratio=0.
                if(bmeas(1,lij,nrij).gt.0.) 
     &               ratio=bmeas(2,lij,nrij)/bmeas(1,lij,nrij)
                call hf1(160+lij,float(nrij)-.5,ratio)
             enddo
          enddo

          stop
          end
