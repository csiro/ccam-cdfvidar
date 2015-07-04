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
    
module cll_m

private
public clon,clat,cllalloc,clldealloc

real, dimension(:), allocatable, save :: clon,clat

contains

subroutine cllalloc(il)

implicit none

integer, intent(in) :: il
integer jl,ifull

jl=6*il
ifull=il*jl

allocate(clat(ifull),clon(ifull))

return
end subroutine cllalloc

subroutine clldealloc

implicit none

deallocate(clat,clon)

return
end subroutine clldealloc

end module cll_m