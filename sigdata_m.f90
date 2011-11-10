module sigdata_m

private
public pmsl,sfct,zs,ps,us,vs,ts,rs,hs,psg_m,zsi_m,lsm_m,sigdataalloc,sigdatadealloc
public fracice

real, dimension(:), allocatable, save :: pmsl,sfct,zs,ps,fracice
real, dimension(:), allocatable, save :: psg_m,zsi_m,lsm_m
real, dimension(:,:), allocatable, save :: us,vs,ts,rs,hs

contains

subroutine sigdataalloc(il,kl)

implicit none

integer, intent(in) :: il,kl
integer jl,ifull

jl=6*il
ifull=il*jl

allocate(pmsl(ifull),sfct(ifull),zs(ifull),ps(ifull),fracice(ifull))
allocate(psg_m(ifull),zsi_m(ifull),lsm_m(ifull))
allocate(us(ifull,kl),vs(ifull,kl),ts(ifull,kl),rs(ifull,kl),hs(ifull,kl))

return
end subroutine sigdataalloc

subroutine sigdatadealloc

implicit none

deallocate(pmsl,sfct,zs,ps,fracice)
deallocate(psg_m,zsi_m,lsm_m)
deallocate(us,vs,ts,rs,hs)

return
end subroutine sigdatadealloc

end module sigdata_m
