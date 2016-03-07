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
      
      subroutine mslp ( ps, pmsl, zs, t, sig, ndim, npts, kl )
      
      implicit none
c------------------------------------------------------------------------------
c
c routine to calculate the mean sea level pressure
c
c   input:
c     ps = sfc.pres.
c     zs = sfc geop in m2/s2
c     t  = temps ( need level 2 temps )
c     sig= sigma values
c     ndim = 1st dim.
c     npts = number of points to do
c     kl = 2nd dim.
c   output:
c     ps = mslp (same units as ps)
c
c------------------------------------------------------------------------------

      real, parameter :: c=9.806/6.5e-3
      real, parameter :: conr=c/287.
      real, parameter :: rconr=287./c

      integer :: ndim, kl, npts
      integer :: i
      real :: con
      
      real, dimension(ndim,kl) :: t
      real, dimension(ndim) :: ps,pmsl,zs
      real, dimension(kl) :: sig
      
      integer ktemp

      ktemp = 2

      con=sig(ktemp)**rconr/c

!     write(6,*)c,conr,rconr,con,ktemp

      do i=1,npts
!       write(6,*)i,ps(i),zs(i),t(i,ktemp)
        pmsl(i) = ps(i)*min(1.+con*zs(i)/t(i,ktemp),1.e5)**conr ! MJT suggestion for single precision
      end do ! i=1,npts

      return
      end
