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
      
      subroutine findxn ( ar, nval, spp, amax, k1, amin, k2 )
c
c determines maximum (amax) of ar located at k1
c        and minimum (amin) of ar located at k2
c
c ar dimensioned nval
c
c ignores any values less than spp
c
      dimension ar(nval)
      integer, dimension(1) :: pos
c
c preset amax and amin
c
      pos = maxloc(ar,ar>=spp)
      k1 = pos(1)
      pos = minloc(ar,ar>=spp)
      k2 = pos(1)
      amax = ar(k1)
      amin = ar(k2)

      write(6,300) amax, k1, amin, k2
 300  format(" x=",f12.6," i=",i6," n=",f12.4," i=",i6)

      return
      end
