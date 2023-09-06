! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2023 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
      
      subroutine fill(a,il,jl,value)

! now  assumes actual spval < value

c     routine fills in interior of an array which has undefined points
      real a(il,jl)         ! input and output array
      real value            ! array value denoting undefined
      real b(il,jl)
      dimension in(8), jn(8)   ! specifies neighbours
      data in/-1,-1,-1,0,1,1, 1, 0/
      data jn/-1, 0, 1,1,1,0,-1,-1/

      if ( all(a(:,:).lt.value) ) then
        write(6,*) "WARN: No valid data for fill"
        return
      end if
      
      write(6,*)"fill il,jl,value=",il,jl,value

      do while ( any(a(2:il-1,2:jl-1)<value) )
       do j=2,jl-1
        do i=2,il-1
         b(i,j)=a(i,j)
         if(a(i,j).lt.value)then
          neighb=0
          av=0.
          do nbs=1,8
           if(a(i+in(nbs),j+jn(nbs)).gt.value)then
            neighb=neighb+1
            av=av+a(i+in(nbs),j+jn(nbs))
           endif
          end do
          if(neighb.gt.0)then
           b(i,j)=av/neighb
          endif
         endif
        end do
       end do 
       do j=2,jl-1
        do i=2,il-1
         a(i,j)=b(i,j)
        enddo
       enddo
      end do ! while (any(a<value))

!     fix up any boundary points
      do i=2,il-1
       if(a(i,1).lt.value)a(i,1)=a(i,2)
       if(a(i,jl).lt.value)a(i,jl)=a(i,jl-1)
      end do
      do j=2,jl-1
       if(a(1,j).lt.value)a(1,j)=a(2,j)
       if(a(il,j).lt.value)a(il,j)=a(il-1,j)
      end do
      if(a(1,1).lt.value)a(1,1)=a(2,2)
      if(a(il,1).lt.value)a(il,1)=a(il-1,2)
      if(a(1,jl).lt.value)a(1,jl)=a(2,jl-1)
      if(a(il,jl).lt.value)a(il,jl)=a(il-1,jl-1)
      return
      end
c******************************************************************************
