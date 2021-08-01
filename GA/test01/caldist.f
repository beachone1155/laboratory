      integer ix(50),iy(50)
      open(11,file="tmp_mod.dat")
      do i=1,50
         read(11,*) idummy,ix(i),iy(i)
c        write(*,*) ix(i),iy(i)
      enddo
      idist=0
      do i=1,49
      idist=idist+int(0.5 +
     & sqrt(real((ix(i+1)-ix(i))**2+(iy(i+1)-iy(i))**2)))
      enddo
      idist=idist+int(0.5 +
     & sqrt(real((ix(1)-ix(50))**2+(iy(1)-iy(50))**2))) 
      write(*,*) "適応度=",-idist
      stop
      end
