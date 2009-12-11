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