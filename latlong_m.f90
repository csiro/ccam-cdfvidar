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