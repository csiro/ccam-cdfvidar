module comsig_m

private
public dsgx,sgmlx,sgx,comsigalloc,comsigdealloc

real, dimension(:), allocatable, save :: dsgx,sgmlx,sgx

contains

subroutine comsigalloc(kl)

implicit none

integer, intent(in) :: kl

allocate(dsgx(kl),sgmlx(kl),sgx(kl+1))

return
end subroutine comsigalloc

subroutine comsigdealloc

implicit none

deallocate(dsgx,sgmlx,sgx)

return
end subroutine comsigdealloc

end module comsig_m
