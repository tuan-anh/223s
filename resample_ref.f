      program resample
      do k=1,40
         call transform(-181.5*.12,.12,i,-2.5,.1,floi,ilo)
         call transform(-177.5*.12,.12,j,-2.5,.1,floj,jlo)
         call transform(-36.5*.2115,.2115,k,-8.,.4,flok,klo)               
         profmeas(1,i+50*(j-1),k)=
     & floi*floj*(flok*co32(klo,jlo,ilo)+(1.-flok)*co32(klo+1,jlo,ilo))+
     & floi*(1.-floj)*(flok*co32(klo,jlo+1,ilo)+(1.-flok)
     & *co32(klo+1,jlo+1,ilo))+
     & (1.-floi)*floj*(flok*co32(klo,jlo,ilo+1)+(1.-flok)
     &             *co32(klo+1,jlo,ilo+1))+
     & (1.-floi)*(1.-floj)*(flok*co32(klo,jlo+1,ilo+1)+(1.-flok)
     &               *co32(klo+1,jlo+1,ilo+1))
         call transform(-220.5*.05,.05,i,-2.5,.1,floi,ilo)
         call transform(-210.5*.05,.05,j,-2.5,.1,floj,jlo)
         call transform(-20.7*.4,.4,k,-8.,.4,flok,klo)
               profmeas(2,i+50*(j-1),k)=
     & floi*floj*(flok*co65(klo,jlo,ilo)+(1.-flok)*co65(klo+1,jlo,ilo))+
     & floi*(1.-floj)*(flok*co65(klo,jlo+1,ilo)+(1.-flok)
     &               *co65(klo+1,jlo+1,ilo))+
     & (1.-floi)*floj*(flok*co65(klo,jlo,ilo+1)+(1.-flok)
     &               *co65(klo+1,jlo,ilo+1))+
     & (1.-floi)*(1.-floj)*(flok*co65(klo,jlo+1,ilo+1)+(1.-flok)
     &               *co65(klo+1,jlo+1,ilo+1))
         sco32=sco32+profmeas(1,i+50*(j-1),k)
         sco65=sco65+profmeas(2,i+50*(j-1),k)
               
      enddo
