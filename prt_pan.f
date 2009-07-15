c===========================================================================
      subroutine prt_pan(var,il,jl,npan,lab)

      real var(il,jl)
      real temp(il,il)
      character*(*) lab
      character*3 cpan

      write(6,*) "prt_pan il,jl,npan,lab=",il,jl,npan,lab

      write(cpan,'("p",i1,"-")')npan

      j1 = 1+il*(npan-1)
      j2 = j1+il-1
      write(6,*) "j1,j2=",j1,j2,il

      do j=j1,j2
       do i=1,il
         !write(6,*)i,j,var(i,j)
         temp(i,j-j1+1) = var(i,j)
       end do
      end do
 
      call amap ( temp, il, il, cpan//lab, 0., 0. )

      return
      end
