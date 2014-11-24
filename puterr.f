      program puterr
 

      call HLIMAP(20000,'TMP')
      call hbook1(10,'tmp',10.,0.,10.,0.)

      do i=1,5
         x=float(i)*2
         call hf1(10,float(i),x)
      enddo
      
      stop
      end
