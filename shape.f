      program shape
      
      do i=1,21
         cphi=-1.1+float(i)*.1
         shapeb=exp(-0.5*(abs(cphi)-1.+.1)**2/.03**2)
         shaped=exp(-0.5*(cphi)**2/.08**2)
         print*,i,cphi,shapeb,shaped
      enddo
      stop
      end
