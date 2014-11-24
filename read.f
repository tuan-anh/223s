      program red_rectangle
      real dat_cub(43,432,432)
      logical cond
      call HLIMAP(20000,'MAP')
      call hbook1(10,'map',100.,-1.,4.5,0.)

ccc interested region: rectangle [170,160] [270,260]
ccc
      open(3,file='CO6-5_channel_image_rgn.txt')
      open(2,file='CO6-5_channel_image.txt')
      do k=1,43
         do j=1,432
            do i=1,432
               read(2,*) kk, jj, ii, aa
               dat_cub(k,j,i) = aa
               call hf1(10,aa,1.)
               cond=(j.ge.160.and.j.lt.260.and.i.ge.170.and.i.lt.270)
               if (cond) then
                  write(3,*) k, j, i, aa
               endif
            enddo
         enddo
      enddo


      stop
      end
