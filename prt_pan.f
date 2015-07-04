! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------

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
