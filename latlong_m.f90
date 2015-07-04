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
    
module latlong_m

private
public rlat,rlong,slat,slon,dlat,dlon,latlongalloc,latlongdealloc

real, save :: slat,slon,dlat,dlon
real, dimension(:), allocatable, save :: rlat,rlong

contains

subroutine latlongalloc(il)

implicit none

integer, intent(in) :: il
integer jl,ifull

jl=6*il
ifull=il*jl

allocate(rlat(ifull),rlong(ifull))

return
end subroutine latlongalloc

subroutine latlongdealloc

implicit none

deallocate(rlat,rlong)

return
end subroutine latlongdealloc

end module latlong_m