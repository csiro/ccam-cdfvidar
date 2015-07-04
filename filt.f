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
      
      subroutine filt ( a, b, n1, n2, we, wc, spp )

c Applies a simple filter over data in 2 dimensions
c   the first pass will take a and put into b
c   the second pass will take b and put back into a

c******************************************************************

      parameter ( c0=0., c2=2. )

c input/output variables
      real a(n1,n2)  ! input/output array
      real b(n1,n2)  ! work array
      real we(1)     ! end weight
      real wc        ! center weight
      integer n1     ! 1st dimension of a,b
      integer n2     ! 2nd dimension of a,b
      real spp       ! minimum acceptable value

c local variables
      real tw        ! sum of all weights
      real cw        ! corrected center weight
      real ew        ! corrected edge weight
      integer n1m    ! n1-1
      integer n2m    ! n2-1

c calculate sum total of weights
      tw = c0
      do k = 1, 1
        tw = tw + c2 * we(k)
      end do ! k=1,1
      tw = tw + wc

c correct weights to sum to 1
      cw = wc/tw
      ew = we(1)/tw

      print *,'filt called with cw,ew=',cw,ew

c set up do loop limits
      n1m=n1-1
      n2m=n2-1

c initialize b array to a
      do j = 1, n2
        do i = 1,n1
          b(i,j) = a(i,j)
        end do ! i=1,n1
      end do ! j=1,n2

c 1st pass- e/w filter of a put into b
      do j = 1, n2
        do i = 2, n1m
          if(a(i-1,j).gt.spp.and.a(i+1,j).gt.spp.and.a(i,j).gt.spp)then
            b(i,j) = ew*(a(i-1,j)+a(i+1,j))+cw*a(i,j)
          endif
        end do ! i=2,n1m
      end do ! j=1,n2

c set boundary values not done in 2nd pass
      do i = 1, n1
        a(i,1) = b(i,1)
        a(i,n2) = b(i,n2)
      end do ! i=1,n1

c 2nd pass- n/s filter of b put into a
      do j = 2, n2m
        do i = 1, n1
          if(b(i,j-1).gt.spp.and.b(i,j+1).gt.spp.and.b(i,j).gt.spp)then
            a(i,j) = ew*(b(i,j-1)+b(i,j+1))+cw*b(i,j)
          else
            a(i,j) = b(i,j)
          endif
        end do ! i=1,n1
      end do ! j=2,n2m

c fix corners
      a(1,1)=.5*(a(1,2)+a(2,1))
      a(n1,1)=.5*(a(n1,2)+a(n1m,1))
      a(1,n2)=.5*(a(1,n2m)+a(2,n2))
      a(n1,n2)=.5*(a(n1,n2m)+a(n1m,n2))

      return
      end
