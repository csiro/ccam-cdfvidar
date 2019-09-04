! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2019 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

module outcdf_m
    
private
public outcdf

integer ixp, iyp, idlev, idnt
integer idproc, idgpnode, idgpoff
integer vnode_nproc, nxp, nyp
integer il_l, jl_l, npan_l, ipan, jpan
integer dlen
integer, dimension(:), allocatable :: ipoff, jpoff, npoff

contains

subroutine outcdf(ihr,idy,imon,iyr,iout,nt,time,mtimer,sig,cdffile_in,ddss,il,kl, &
                  minlon,maxlon,minlat,maxlat,llrng,procformat_nproc)

use netcdf_m
      
implicit none
      
integer il,jl,kl,ifull
integer ihr,idy,imon,iyr,iout
integer nt
real ddss
real du,tanl,rnml,stl1,stl2

common/mapproj/du,tanl,rnml,stl1,stl2

character(len=10) rundate
      
integer, intent(in) :: procformat_nproc

integer, save :: idnc0, idnc1, idncm1
integer kdate, ktime, iarch
integer idnc, ier, imode
integer ktau, icy, icm, icd
integer ich, icmi, ics, mtimer, idv
integer, parameter :: nextout=0
integer, parameter :: itype=1
integer, parameter :: nihead=54
integer nahead(nihead)
integer, dimension(1) :: dimids
real, intent(in) :: minlon, maxlon, minlat, maxlat, llrng
logical ptest

integer, parameter :: nrhead = 14
real ahead(nrhead)

real, dimension(kl) :: sig
real time,dt,ds

character(len=1024) cdffile_in
character(len=1032) cdffile

!common/cdfind/ixp,iyp,idlev,idnt

integer, dimension(5) :: dim
integer xdim,ydim,zdim,tdim
integer pdim,gpdim
integer oldmode
character(len=20) :: timorg
character(len=33) :: grdtim
character(len=3), dimension(12) :: month
data month/'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'/

data idnc1/0/, idnc0/0/, idncm1/0/
data rundate/"ncepavnanl"/

jl=6*il
ifull=il*jl
il_l=il ! default
jl_l=jl ! default

nahead = 0
ahead = 0.
ddss = 0.

write(6,*)"outcdf ihr,idy,imon,iyr,iout=",ihr,idy,imon,iyr,iout
write(6,*)"time=",time

dt=0
kdate=iyr*10000+imon*100+idy
ktime=ihr*100
write(6,*)"kdate,ktime=",kdate,ktime

! itype=1 outfile
iarch=nt
idnc=idnc1

write(6,'("outcdf itype,idnc,iarch,cdffile=",3i5," ",a80)') itype,idnc,iarch,cdffile_in

if ( procformat_nproc>0 ) then
  ! assume face decomposition
  vnode_nproc = procformat_nproc
  ptest = .true.
  do while ( ptest )
    if ( mod(vnode_nproc,6)==0 .or. mod(6,vnode_nproc)==0 ) then
      nxp = max( 1, nint(sqrt(real(vnode_nproc)/6.)) ) ! number of proc in x-direction
      nyp = vnode_nproc/nxp                            ! number of proc in y-direction
      do while ( (mod(il,max(nxp,1))/=0 .or. mod(vnode_nproc/6,max(nxp,1))/=0 .or. &
                  mod(jl,max(nyp,1))/=0) .and. nxp>0 )
        nxp = nxp - 1
        nyp = vnode_nproc/max(nxp,1)
      end do
      if ( nxp>0 ) then
        ptest = .false.
      else
        vnode_nproc = vnode_nproc - 1  
      end if    
    else 
      vnode_nproc = vnode_nproc - 1  
    end if
  end do
  write(6,*) "Procformat nproc,nxp,nyp ",vnode_nproc,nxp,nyp
  dlen = 5
else
  vnode_nproc = 0
  nxp = 1
  nyp = 1
  dlen = 4
end if
      
!#######################################################################
! netcdf output
!#######################################################################

if ( iarch.lt.1 ) then
  write(6,*) "wrong iarch in outcdf"
  call finishbanner
  stop
end if
if ( iarch.eq.1 ) then
        
  ! expect to create new file
  cdffile = ""
  if ( vnode_nproc>0 ) then
    write(cdffile,"(a,'.',i6.6)") trim(cdffile_in), 0  
  else  
    cdffile = trim(cdffile_in)  
  end if    
  write(6,*)'nccre of ',cdffile
        
#ifdef usenc3
#ifdef no64bit_offset
  ier = nf_create(cdffile,NF_CLOBBER,idnc)
#else
  ier = nf_create(cdffile,NF_64BIT_OFFSET,idnc)
#endif
#else
  ier = nf_create(cdffile,NF_NETCDF4,idnc)
#endif
  write(6,*)'idnc,ier=',idnc,ier
  ! Turn off the data filling
  ier = nf_set_fill(idnc,nf_nofill,oldmode)
  ! Create dimensions, lon, lat
  if ( vnode_nproc>0 ) then
    il_l = il/nxp  
    ier = nf_def_dim(idnc,'longitude', il_l,         xdim)
    jl_l = jl/nyp
    ier = nf_def_dim(idnc,'latitude',  jl_l,         ydim)       
  else
    ier = nf_def_dim(idnc,'longitude', il,           xdim)
    ier = nf_def_dim(idnc,'latitude',  jl,           ydim)
  end if   
  ier = nf_def_dim(idnc,'lev',       kl,           zdim)
  if ( vnode_nproc>0 ) then
    ier = nf_def_dim(idnc,'processor',vnode_nproc,pdim)
    ier = nf_def_dim(idnc,'gprocessor',vnode_nproc,gpdim)
  else
    pdim = 0
    gpdim = 0
  end if
  ier = nf_def_dim(idnc,'time',      nf_unlimited, tdim)
  write(6,*) "xdim=",xdim," ydim=",ydim," zdim=",zdim," tdim=",tdim
          
  ! define coords.
          
  if ( vnode_nproc>0 ) then
    dim(1) = xdim
    dim(2) = pdim
    ier = nf_def_var(idnc,'longitude',nf_float,2,dim(1:2),ixp)
    ier = nf_put_att_text(idnc,ixp,'point_spacing',4,'even')
    ier = nf_put_att_text(idnc,ixp,'units',12,'degrees_east')
    dim(1) = ydim
    dim(2) = pdim
    ier = nf_def_var(idnc,'latitude',nf_float,2,dim(1:2),iyp)
    ier = nf_put_att_text(idnc,iyp,'point_spacing',4,'even')
    ier = nf_put_att_text(idnc,iyp,'units',13,'degrees_north')   
  else
    dimids = xdim
    ier = nf_def_var(idnc,'longitude',nf_float,1,dimids,ixp)
    ier = nf_put_att_text(idnc,ixp,'point_spacing',4,'even')
    ier = nf_put_att_text(idnc,ixp,'units',12,'degrees_east')
    dimids = ydim
    ier = nf_def_var(idnc,'latitude',nf_float,1,dimids,iyp)
    ier = nf_put_att_text(idnc,iyp,'point_spacing',4,'even')
    ier = nf_put_att_text(idnc,iyp,'units',13,'degrees_north')
  end if  
  write(6,*)'ixp,iyp=',ixp,iyp

  dimids = zdim
  ier = nf_def_var(idnc,'lev',nf_float,1,dimids,idlev)
  write(6,*)'idlev,ier=',idlev,ier
  ier = nf_put_att_text(idnc,idlev,'positive',4,'down')
  ier = nf_put_att_text(idnc,idlev,'point_spacing',6,'uneven')
  ier = nf_put_att_text(idnc,idlev,'units',11,'sigma_level')
  ier = nf_put_att_text(idnc,idlev,'long_name',11,'sigma_level')
  write(6,*)'idlev=',idlev

  if ( vnode_nproc>0 ) then
    dim(1) = pdim
    ier = nf_def_var(idnc,'processor',nf_int,1,dim(1),idproc)
    ier = nf_put_att_text(idnc,idproc,'long_name',16,'processor number')
    dim(1) = gpdim
    ier = nf_def_var(idnc,'gprocnode',nf_int,1,dim(1),idgpnode)
    ier = nf_put_att_text(idnc,idgpnode,'long_name',25,'global processor node map')
    ier = nf_def_var(idnc,'gprocoffset',nf_int,1,dim(1),idgpoff)
    ier = nf_put_att_text(idnc,idgpnode,'long_name',27,'global processor offset map')
  end if

  write(6,*)'tdim,idnc=',tdim,idnc
  dimids = tdim
  ier = nf_def_var(idnc,'time',nf_int,1,dimids,idnt)
  write(6,*)'idnt=',idnt
  ier = nf_put_att_text(idnc,idnt,'point_spacing',4,'even')

  write(6,*)'kdate,ktime=',kdate,ktime

  icy=kdate/10000
  icm=max(1,min(12,(kdate-icy*10000)/100))
  icd=max(1,min(31,(kdate-icy*10000-icm*100)))
  ich=ktime/100
  icmi=(ktime-ich*100)
  ics=0
  write(6,*) icy,icm,icd,ich,icmi,ics
  write(timorg,'(i2.2,"-",a3,"-",i4.4,3(":",i2.2))') icd,month(icm),icy,ich,icmi,ics
  write(6,*)'timorg=',timorg
  ier = nf_put_att_text(idnc,idnt,'time_origin',20,timorg)

  write(grdtim,'("minutes since ",i4.4,"-",i2.2,"-",i2.2," ",2(i2.2,":"),i2.2)') icy,icm,icd,ich,icmi,ics
  write(6,*)'grdtim=',grdtim
  ier = nf_put_att_text(idnc,idnt,'units',33,grdtim)

  if ( vnode_nproc>0 ) then
    dim(1) = xdim
    dim(2) = ydim
    dim(3) = zdim
    dim(4) = pdim
    dim(5) = tdim
  else
    dim(1) = xdim
    dim(2) = ydim
    dim(3) = zdim
    dim(4) = tdim
    dim(5) = 0
  end if
  write(6,*) "dim=",dim

  ! create the attributes of the header record of the file
  nahead(1)=il         ! needed by cc2hist
  nahead(2)=jl         ! needed by cc2hist
  nahead(3)=kl         ! needed by cc2hist
  nahead(4)=0 !m
  nahead(5)=0 !nsd        ! not needed now
  nahead(6)=0 !io_in
  nahead(7)=0 !nbd
  nahead(8)=0 !nps
  nahead(9)=0 !mex
  nahead(10)=0 !mup
  nahead(11)=0 !nem
  nahead(12)=mtimer ! should be mtimer
  nahead(13)=0 ! nmi
  nahead(14)=0 !ndt       ! needed by cc2hist
  nahead(15)=0 !npsav
  nahead(16)=0 !nhor
  nahead(17)=0 !nkuo
  nahead(18)=0 !khdif
  nahead(19)=0 !kwt       ! needed by cc2hist
  nahead(20)=0 !iaa
  nahead(21)=0 !jaa
  nahead(22)=0 !nvad
  nahead(23)=0 !nqg       ! not needed now      
  nahead(24)=0 !lbd
  nahead(25)=0 !nrun
  nahead(26)=0 !nrunx
  nahead(27)=0 !khor
  nahead(28)=0 !ksc
  nahead(29)=0 !kountr
  nahead(30)=0 !ndiur
  nahead(31)=0 !nspare
  nahead(32)=0 !nhorps
  nahead(33)=0 !nsoil
  nahead(34)=0 !ms        ! needed by cc2hist
  nahead(35)=0 !ntsur
  nahead(36)=0 !nrad
  nahead(37)=0 !kuocb
  nahead(38)=0 !nvmix
  nahead(39)=0 !ntsea
  nahead(40)=0
  nahead(41)=0 !nextout
  nahead(42)=0 !ilt
  nahead(43)=0 !ntrac     ! needed by cc2hist
  nahead(44)=0 !nsib
  nahead(45)=0 !nrungcm
  nahead(46)=0 !ncvmix
  nahead(47)=0 !ngwd
  nahead(48)=0 !lgwd
  nahead(49)=0 !mup
  nahead(50)=0 !nritch
  nahead(51)=0 !ifullw
  nahead(52)=0 !nevapls
  nahead(53)=0 !nevapcc
  nahead(54)=0.!nhadq
  write(6,'("nahead=",(20i4))') nahead
  ds=ddss
  ahead(1)=ds
  ahead(2)=0.  !difknbd
  ahead(3)=0.  ! was rhkuo for kuo scheme
  ahead(4)=0.  !du
  ahead(5)=du ! rlong0     ! needed by cc2hist
  ahead(6)=tanl ! rlat0      ! needed by cc2hist
  ahead(7)=rnml ! schmidt    ! needed by cc2hist
  ahead(8)=0.  !stl2
  ahead(9)=0.  !relaxt
  ahead(10)=0.  !hourbd
  ahead(11)=0. !tss_sh
  ahead(12)=0 !vmodmin
  ahead(13)=0 !av_vmod
  ahead(14)=0 !epsp
  write(6,*) "ahead=",ahead
  ier = nf_put_att_int(idnc,nf_global,'int_header',nf_int,nihead,nahead)
  if(ier.ne.0)write(6,*)"ncapt int idnc,ier=",idnc,ier
  ier = nf_put_att_real(idnc,nf_global,'real_header',nf_float,nrhead,ahead)
  if(ier.ne.0)write(6,*)"ncapt real idnc,ier=",idnc,ier
  ier = nf_put_att_text(idnc,nf_global,'date_header',10,rundate)
  if(ier.ne.0)write(6,*)"ncaptc date idnc,ier=",idnc,ier
        
endif ! ( iarch=1 ) then

write(6,*)'call openhist for itype= ',itype
call openhist(idnc,iarch,itype,dim,sig,kdate,ktime,time,mtimer,il,kl, &
              minlon,maxlon,minlat,maxlat,llrng)

! MJT notes - nf_sync is used for multiple time-steps in output file
ier = nf_sync(idnc)
if(ier.ne.0)write(6,*)"ncsnc idnc,ier=",idnc,ier
      
!ier = nf_close(idnc)

if ( itype.eq.1 ) then
! itype=1 outfile
  idnc1=idnc
elseif ( itype.eq.0 ) then
! itype=0 climcdf
  idnc0=idnc
elseif ( itype.eq.-1 ) then
! itype=-1 restfile
  idncm1=idnc
endif ! ( itype.eq.1 ) then

return ! outcdf
end subroutine outcdf
!=======================================================================
subroutine openhist(idnc,iarch,itype,dim,sig,kdate,ktime,time,mtimer,il,kl, &
    minlon,maxlon,minlat,maxlat,llrng)

use cll_m
use netcdf_m
use sigdata_m
!use xyzinfo_m, only : em,f

!     this routine creates attributes and writes output

integer il,jl,kl,ifull
integer ioff, joff, noff, np, myface, mtmp

character(len=50) :: lname, expdesc
integer dim(5)
integer idim2(4)
integer id2len
integer, dimension(1) :: ivals
integer, dimension(1) :: start, ncount
integer, dimension(4) :: ostart, ocount
integer, dimension(:), allocatable :: procnode, procoffset, vnode_dat
integer nrun
real, dimension(:), allocatable :: xpnt,ypnt
real, dimension(kl), intent(in) :: sig
real, intent(in) :: minlon, maxlon, minlat, maxlat, llrng

!common/cdfind/ixp,iyp,idlev,idnt
real, dimension(:), allocatable :: tst,tsb
!       *** qscrn_ave not presently written     
real, dimension(:), allocatable :: aa,bb,cc
real, dimension(:,:), allocatable :: cfrac
real, dimension(1) :: rvals
      
integer iq, varid
integer i, j, n
real x_wgt, y_wgt
real, dimension(:), allocatable :: cc_wgt, origdata2d
real, dimension(:,:), allocatable :: origdata3d

jl=6*il
ifull=il*jl
      
if ( vnode_nproc>=6 ) then
  npan_l = 1
  allocate( ipoff(vnode_nproc), jpoff(vnode_nproc), npoff(vnode_nproc) )
else if ( vnode_nproc>0 ) then
  npan_l = 6/vnode_nproc
  allocate( ipoff(vnode_nproc), jpoff(vnode_nproc), npoff(vnode_nproc) )
else
  npan_l = 1
end if

id2len = dlen - 1
il_l = il/nxp
jl_l = jl/nyp
ipan = il/nxp
jpan = jl/nyp
if ( vnode_nproc>=6 ) then
  do np = 0,vnode_nproc-1
    myface = np*6/(nxp*nyp)
    mtmp = np - myface*nxp*nyp/6
    npoff(np+1) = 1 - myface
    jpoff(np+1) = (mtmp/nxp)*jpan
    ipoff(np+1) = modulo( mtmp, nxp )*ipan
  end do  
  write(6,*) "ipan,jpan,npan_l=",ipan,jpan,npan_l
  write(6,*) "ipoff=",ipoff
  write(6,*) "jpoff=",jpoff
  write(6,*) "npoff=",npoff
else if ( vnode_nproc>0 ) then
  do np = 0,vnode_nproc-1
    npoff(np+1) = 1 - np*npan_l
    ipoff(np+1) = 0
    jpoff(np+1) = 0
  end do    
  write(6,*) "ipan,jpan,npan_l=",ipan,jpan,npan_l
  write(6,*) "ipoff=",ipoff
  write(6,*) "jpoff=",jpoff
  write(6,*) "npoff=",npoff
end if
      
allocate( xpnt(il),ypnt(6*il) )
allocate( tst(6*il*il),tsb(6*il*il) )
allocate( aa(6*il*il),bb(6*il*il),cc(6*il*il) )
allocate( cfrac(6*il*il,kl) )
allocate( cc_wgt(6*il*il),origdata2d(6*il*il) )
allocate( origdata3d(6*il*il,kl) )

write(6,*)'openhist iarch,idnc,itype=',iarch,idnc,itype

!if this is the first archive, set up some global attributes
if(iarch.eq.1) then
  write(6,*)'dim=',dim(1:dlen)
  if ( vnode_nproc>0 ) then
    idim2(1)=dim(1)    ! x
    idim2(2)=dim(2)    ! y
    idim2(3)=dim(4)    ! proc
    idim2(4)=dim(dlen) ! time
  else    
    idim2(1)=dim(1)    ! x
    idim2(2)=dim(2)    ! y
    idim2(3)=dim(dlen) ! time
  end if  
  write(6,*)'idim2=',idim2
  write(6,*)'id2len=',id2len

! Create global attributes
! Model run number
  nrun = 1
  ivals = nrun
  ier = nf_put_att_int(idnc,nf_global,'nrun',nf_int,1,ivals)
  write(6,*)"nrun=",nrun," ier=",ier

! Experiment description
  expdesc = 'CCAM model run'
  ier = nf_put_att_text(idnc,nf_global,'expdesc',50,expdesc)
  write(6,*)"expdesc=",expdesc," ier=",ier

  if ( vnode_nproc>0 ) then
    ier = nf_put_att_int(idnc,nf_global,'nproc',nf_int,1,vnode_nproc)  
    ier = nf_put_att_int(idnc,nf_global,'procmode',nf_int,1,vnode_nproc)
    ier = nf_put_att_text(idnc,nf_global,'decomp',4,'face')
  end if
        
! Sigma levels
! write(6,*)'sig=',sig
  ier = nf_put_att_real(idnc,nf_global,'sigma',nf_float,kl,sig)

  lname = 'year-month-day at start of run'
  ier = nf_def_var(idnc,'kdate',nf_int,1,dim(4:4),idkdate)
  ier = nf_put_att_text(idnc,idkdate,'long_name',len_trim(lname),lname)

  lname = 'hour-minute at start of run'
  ier = nf_def_var(idnc,'ktime',nf_int,1,dim(4:4),idktime)
  ier = nf_put_att_text(idnc,idktime,'long_name',len_trim(lname),lname)

  lname = 'timer (hrs)'
  ier = nf_def_var(idnc,'timer',nf_float,1,dim(4:4),idnter)
  ier = nf_put_att_text(idnc,idnter,'long_name',len_trim(lname),lname)

  lname = 'mtimer (mins)'
  ier = nf_def_var(idnc,'mtimer',nf_int,1,dim(4:4),idmtimer)
  ier = nf_put_att_text(idnc,idmtimer,'long_name',len_trim(lname),lname)

  lname = 'timeg (UTC)'
  ier = nf_def_var(idnc,'timeg',nf_float,1,dim(4:4),idnteg)
  ier = nf_put_att_text(idnc,idnteg,'long_name',len_trim(lname),lname)

  lname = 'number of time steps from start'
  ier = nf_def_var(idnc,'ktau',nf_int,1,dim(4:4),idktau)
  ier = nf_put_att_text(idnc,idktau,'long_name',len_trim(lname),lname)

  ier = nf_def_var(idnc,'sigma',nf_float,1,dim(3:3),idv)
  ier = nf_put_att_text(idnc,idv,'positive',4,'down')

  write(6,*)'define attributes of variables'

  lname ='Scaled Log Surface pressure'
  call attrib(idnc,idim2,id2len,'psf',lname,'none',-1.3,0.2)

  lname ='Mean sea level pressure'
  call attrib(idnc,idim2,id2len,'pmsl',lname,'hPa',800.,1200.)
  lname = 'Surface geopotential'
  call attrib(idnc,idim2,id2len-1,'zht',lname,'m2/s2',-2.e3,128.e3) ! MJT lsmask
  lname = 'Soil type'                                               ! MJT lsmask
  call attrib(idnc,idim2,id2len-1,'soilt',lname,'none',0.,65.e3)    ! MJT lsmask

! For time invariant surface fields
  lname = 'clon'
  call attrib(idnc,idim2,id2len-1,'clon',lname,'none',-360.,360.)
  lname = 'clat'
  call attrib(idnc,idim2,id2len-1,'clat',lname,'none',-90.,90.)

! For time varying surface fields
  lname = 'Surface temperature'
  call attrib(idnc,idim2,id2len,'tsu',lname,'K',0.,350.)

  if ( all(soiltemp<0.) ) then
    lname = 'Soil temperature top'
    call attrib(idnc,idim2,id2len,'tb3',lname,'K',100.,400.)
    lname = 'Soil temperature bottom'
    call attrib(idnc,idim2,id2len,'tb2',lname,'K',100.,400.)
  else
    lname = 'Soil temperature lev 1'
    call attrib(idnc,idim2,id2len,'tgg1',lname,'K',100.,400.)
    lname = 'Soil temperature lev 2'
    call attrib(idnc,idim2,id2len,'tgg2',lname,'K',100.,400.)
    lname = 'Soil temperature lev 3'
    call attrib(idnc,idim2,id2len,'tgg3',lname,'K',100.,400.)
    lname = 'Soil temperature lev 4'
    call attrib(idnc,idim2,id2len,'tgg4',lname,'K',100.,400.)
    lname = 'Soil temperature lev 5'
    call attrib(idnc,idim2,id2len,'tgg5',lname,'K',100.,400.)
    lname = 'Soil temperature lev 6'
    call attrib(idnc,idim2,id2len,'tgg6',lname,'K',100.,400.)
  end if
      
  if ( all(soilmoist<0.) ) then
    lname = 'Soil moisture top'
    call attrib(idnc,idim2,id2len,'wfg',lname,'none',0.,.4)
    lname = 'Soil moisture bottom'
    call attrib(idnc,idim2,id2len,'wfb',lname,'none',0.,.4)
  else
    lname = 'Soil moisture 1'
    call attrib(idnc,idim2,id2len,'wb1',lname,'m3/m3',0.,2.)
    lname = 'Soil moisture 2'
    call attrib(idnc,idim2,id2len,'wb2',lname,'m3/m3',0.,2.)
     lname = 'Soil moisture 3'
    call attrib(idnc,idim2,id2len,'wb3',lname,'m3/m3',0.,2.)
    lname = 'Soil moisture 4'
    call attrib(idnc,idim2,id2len,'wb4',lname,'m3/m3',0.,2.)
    lname = 'Soil moisture 5'
    call attrib(idnc,idim2,id2len,'wb5',lname,'m3/m3',0.,2.)
    lname = 'Soil moisture 6'
    call attrib(idnc,idim2,id2len,'wb6',lname,'m3/m3',0.,2.)
  end if

  if (any(fracice>=0.)) then
    lname = 'Sea ice fraction'
    call attrib(idnc,idim2,id2len,'fracice',lname,'none',0.,6.5)
  end if
        
  if (any(snod>1.e-8)) then
    lname = 'Snow depth (liquid water)'
    call attrib(idnc,idim2,id2len,'snd',lname,'none',0.,6500.)
  end if

  write(6,*)'3d variables'
  call attrib(idnc,dim,dlen,'temp','Air temperature','K',100.,350.)
  call attrib(idnc,dim,dlen,'u','x-component wind','m/s',-150.,150.)
  call attrib(idnc,dim,dlen,'v','y-component wind','m/s',-150.,150.)
  lname= 'Water mixing ratio'
  call attrib(idnc,dim,dlen,'mixr',lname,'kg/kg',0.,.05)

  write(6,*)'finished defining attributes'
! Leave define mode
  ier = nf_enddef(idnc)
  write(6,*)'leave define mode: ier=',ier

  if ( vnode_nproc>0 ) then
    il_l = il/nxp
    jl_l = jl/nyp
    do np = 0,vnode_nproc-1
      ioff = ipoff(np+1)
      joff = jpoff(np+1)
      noff = npoff(np+1)
      do i = 1,ipan
        xpnt(i) = float(i + ioff)
      end do
      ostart(1) = 1
      ostart(2) = np + 1
      ocount(1) = il_l
      ocount(2) = 1
      ier = nf_put_vara_real(idnc,ixp,ostart(1:2),ocount(1:2),xpnt(1:il_l))
      do n = 1,npan_l
        do j = 1,jpan
          i = j + (n-noff)*jpan 
          ypnt(i) = float(j + joff + (n-noff)*il)
        end do  
      end do
      ostart(1) = 1
      ostart(2) = np + 1
      ocount(1) = jl_l
      ocount(2) = 1
      ier = nf_put_vara_real(idnc,iyp,ostart(1:2),ocount(1:2),ypnt(1:jl_l))
    end do  
  else
    do i=1,il
      xpnt(i) = float(i)
    end do
    start = 1
    ncount = il
    ier = nf_put_vara_real(idnc,ixp,start,ncount,xpnt)
    do j=1,jl
      ypnt(j) = float(j)
    end do
    start = 1
    ncount = jl
    ier = nf_put_vara_real(idnc,iyp,start,ncount,ypnt)
  end if      

  start = 1
  ncount = kl
  ier = nf_put_vara_real(idnc,idlev,start,ncount,sig)
  write(6,*) "sig=",sig

  ier = nf_inq_varid(idnc,'sigma',idv)
  start = 1
  ncount = kl
  ier = nf_put_vara_real(idnc,idv,start,ncount,sig)

  if ( vnode_nproc>0 ) then
    ! store local processor id in output file  
    allocate( vnode_dat(vnode_nproc) )
    do np = 0,vnode_nproc-1
      vnode_dat(np+1) = np
    end do  
    ier = nf_put_vara_int(idnc,idproc,(/1/),(/vnode_nproc/),vnode_dat)
    deallocate( vnode_dat )
    ! store file id for a given processor number in output file number 000000
    allocate( procnode(vnode_nproc) )
    do np = 0,vnode_nproc-1
      procnode(np+1) = 0
    end do  
    ier = nf_put_vara_int(idnc,idgpnode,(/1/),(/vnode_nproc/),procnode)
    deallocate( procnode )
    ! store offset within a file for a given processor number in output file number 000000
    allocate( procoffset(vnode_nproc) )
    do np = 0,vnode_nproc-1
      procoffset(np+1) = np
    end do  
    ier = nf_put_vara_int(idnc,idgpoff,(/1/),(/vnode_nproc/),procoffset)
    deallocate( procoffset )
  end if

endif ! iarch.eq.1
     
!------------------------------------------------------------------      
ktau=0
!timer=0
timer=time/60. ! MJT quick fix
!timeg=0
timeg=mod(time/60.,24.) ! MJT quick fix
write(6,*)'outcdf processing kdate,ktime,ktau,time,mtimer: ',kdate,ktime,ktau,time,mtimer

!set time to number of minutes since start
ier = nf_inq_varid(idnc,'time',idv)
start = iarch
ier = nf_put_var1_int(idnc,idv,start,int(time))
write(6,*)"int(time)=",int(time)

ier = nf_inq_varid(idnc,'kdate',idv)
ier = nf_put_var1_int(idnc,idv,start,kdate)
ier = nf_inq_varid(idnc,'ktime',idv)
ier = nf_put_var1_int(idnc,idv,start,ktime)
ier = nf_inq_varid(idnc,'timer',idv)
ier = nf_put_var1_real(idnc,idv,start,timer)   
ier = nf_inq_varid(idnc,'mtimer',idv)
ier = nf_put_var1_int(idnc,idv,start,nint(time)) ! MJT quick fix
ier = nf_inq_varid(idnc,'timeg',idv)
ier = nf_put_var1_real(idnc,idv,start,timeg)
ier = nf_inq_varid(idnc,'ktau',idv)
ier = nf_put_var1_int(idnc,idv,start,ktau)      

write(6,*)'kdate,ktime,ktau=',kdate,ktime,ktau
write(6,*)'timer,timeg=',timer,timeg

write(6,*)'now write out variables'

if(ktau.eq.0.or.itype.eq.-1)then  ! also for restart file
  call histwrt3(clon,'clon',idnc,iarch,il)
  call histwrt3(clat,'clat',idnc,iarch,il)
endif ! (ktau.eq.0) 

do iq=1,ifull
  aa(iq)=log(ps(iq)/1.e5)
enddo
call histwrt3(aa,'psf',idnc,iarch,il)

is=il/2+(il+il/2-1)*il
xsfct=-1.
xpmsl=-1.
do iq=1,ifull
  ps(iq)=ps(iq)/100.
  zsi_m(iq)=zs(iq)*9.80616
  xsfct=max(xsfct,sfct(iq))
  xpmsl=max(xpmsl,pmsl(iq))
enddo ! iq=1,ifull

write(6,*)"ps(hPa)=",(ps(is+i),i=1,5)
write(6,*)"zs(m2/s2)=",(zsi_m(is+i),i=1,5)
write(6,*)"temp(k=2)=",(ts(is+i,2),i=1,5)
call prt_pan(zs,il,jl,2,'zs(m)')
call prt_pan(zsi_m,il,jl,2,'zs*g(m2/s2)')

call histwrt3(zsi_m,'zht',idnc,iarch,il)   ! always from 13/9/02
call histwrt3(lsm_m*65.e3,'soilt',idnc,iarch,il)

if(xpmsl.lt.500.)then
  write(6,*)"call mslp(ps,pmsl,zs,ts,sig,ifull,ifull,kl)"
  call mslp(ps,pmsl,zsi_m,ts,sig,ifull,ifull,kl)
endif!(xpmsl.lt.500.e2)then

write(6,*)"pmsl=",(pmsl(is+i),i=1,5)
call histwrt3(pmsl,'pmsl',idnc,iarch,il)

if ( xsfct .lt. 200. ) then
  write(6,*)"###################################################"
  write(6,*)"setting sfct temp lowest model level temperature!!!"
  do iq=1,ifull
    sfct(iq)= ts(iq,1)
  enddo ! iq=1,ifull
endif ! ( xsfct .lt. 200. ) then

if ( any( sfct<100. .or. sfct>400. ) ) then
  write(6,*) "ERROR: Invalid sfct"
  write(6,*) "minval,maxval ",minval(sfct),maxval(sfct)
  call finishbanner
  stop -1
end if 
      
call histwrt3(sfct,'tsu',idnc,iarch,il)

if ( all(soiltemp<0.) ) then
  soiltemp(:,1) = ts(:,2)
  soiltemp(:,2) = ts(:,2)
  call histwrt3(soiltemp(:,2),'tb3',idnc,iarch,il) ! top
  call histwrt3(soiltemp(:,1),'tb2',idnc,iarch,il) ! bottom
else
  call histwrt3(soiltemp(:,1),'tgg1',idnc,iarch,il)
  call histwrt3(soiltemp(:,2),'tgg2',idnc,iarch,il)
  call histwrt3(soiltemp(:,3),'tgg3',idnc,iarch,il)
  call histwrt3(soiltemp(:,4),'tgg4',idnc,iarch,il)
  call histwrt3(soiltemp(:,5),'tgg5',idnc,iarch,il)
  call histwrt3(soiltemp(:,6),'tgg6',idnc,iarch,il)
end if
      

if ( all(soilmoist<0.) ) then
  soilmoist(:,1) = 0.14
  soilmoist(:,2) = 0.14
  call histwrt3(soilmoist(:,1),'wfg',idnc,iarch,il)
  call histwrt3(soilmoist(:,2),'wfb',idnc,iarch,il)
else
  call histwrt3(soilmoist(:,1),'wb1',idnc,iarch,il)
  call histwrt3(soilmoist(:,2),'wb2',idnc,iarch,il)
  call histwrt3(soilmoist(:,3),'wb3',idnc,iarch,il)
  call histwrt3(soilmoist(:,4),'wb4',idnc,iarch,il)
  call histwrt3(soilmoist(:,5),'wb5',idnc,iarch,il)
  call histwrt3(soilmoist(:,6),'wb6',idnc,iarch,il)
end if
      
if ( any( fracice >= 0. ) ) then
  call histwrt3(fracice,'fracice',idnc,iarch,il)
end if
      
if ( any( snod > 1.e-8 ) ) then
  call histwrt3(snod,'snd',idnc,iarch,il)
end if
      
write(6,*)'netcdf save of 3d variables'
call histwrt4(ts,'temp',idnc,iarch,il,kl)
call histwrt4(us,'u',idnc,iarch,il,kl)
call histwrt4(vs,'v',idnc,iarch,il,kl)
call histwrt4(rs,'mixr',idnc,iarch,il,kl)

deallocate( xpnt,ypnt )
deallocate( tst,tsb )
deallocate( aa,bb,cc )
deallocate( cfrac )
deallocate( cc_wgt,origdata2d )
deallocate( origdata3d )
if ( vnode_nproc>0 ) then
  deallocate( ipoff, jpoff, npoff )  
end if
      
return ! subroutine openhist(idnc,iarch,itype,dim,sig
end subroutine openhist
!=======================================================================
subroutine attrib(cdfid,dim,ndim,name,lname,units,xmin,xmax)

use netcdf_m

integer(kind=2) minv, maxv, missval   ! was integer*2
parameter(minv = -32500, maxv = 32500, missval = -32501)
integer ndim
integer cdfid, idv, dim(ndim)
integer, dimension(ndim) :: chunks
character(len=*) name, lname, units
real xmin, xmax
real, dimension(1) :: rvals
integer(kind=2), dimension(1) :: i2vals

ier = nf_def_var(cdfid, name, nf_int2, ndim, dim, idv)
if ( ier.ne.0 ) then
  write(6,*)ier,' Error in variable declaration ', name
  write(6,*) nf_strerror(ier)
  call finishbanner
  stop
end if

if ( vnode_nproc>0 ) then
  select case(ndim)
    case(5)
      chunks = (/ dim(1), dim(2), 1, vnode_nproc, 1 /)
    case(4)
      chunks = (/ dim(1), dim(2), vnode_nproc, 1 /)  
    case(3)
      chunks = (/ dim(1), dim(2), vnode_nproc /)
    case default
      write(6,*) "ERROR: Invlid ndim in attrib ",ndim
      call finishbanner
      stop
  end select   
  ier = nf_def_var_chunking(cdfid,idv,nf_chunked,chunks)
  if ( ier.ne.0 ) then
    write(6,*)ier,' Error in variable chunking ', name
    write(6,*) nf_strerror(ier)
    call finishbanner
    stop
  end if
end if

ier = nf_put_att_text(cdfid,idv,'long_name',len_trim(lname),lname)
if(len_trim(units)/=0)then
  ier = nf_put_att_text(cdfid,idv,'units',len_trim(units),units)
end if
i2vals = minv
ier = nf_put_att_int2(cdfid,idv,'valid_min',nf_int2,1,i2vals)
i2vals = maxv
ier = nf_put_att_int2(cdfid,idv,'valid_max',nf_int2,1,i2vals)
i2vals = missval
ier = nf_put_att_int2(cdfid,idv,'missing_value',nf_int2,1,i2vals)
!scalef=(xmax-xmin)/float(maxv - minv)
scalef=(xmax-xmin)/(real(maxv)-real(minv)) ! jlm fix for precision problems
addoff=xmin-scalef*minv
rvals = addoff
ier = nf_put_att_real(cdfid,idv,'add_offset',nf_float,1,rvals)
rvals = scalef
ier = nf_put_att_real(cdfid,idv,'scale_factor',nf_float,1,rvals)
ier = nf_put_att_text(cdfid,idv,'FORTRAN_format',5,'G11.4')
return
end subroutine attrib
!=======================================================================
subroutine histwrt3(var,sname,idnc,iarch,il)
! Write 2d+t fields from the savegrid array.

use netcdf_m
      
implicit none
      
integer il,jl,ifull
integer iarch

!include 'newmpar.h'
!include 'parm.h'

integer mid, start(4), count(4)
integer imn, imx, jmn, jmx
integer i, j, ndims
integer np, ioff, joff, noff, n
!integer(kind=2), dimension(:,:), allocatable :: ipack ! was integer*2 
character(len=*), intent(in) :: sname
!character*8 sname
integer(kind=2) minv, maxv, missval ! was integer*2 
parameter(minv = -32500, maxv = 32500, missval = -32501)
real addoff, scale_f
real varn, varx, xmin, xmax, pvar
integer, intent(in) :: idnc
integer ier

real, dimension(il,6*il), intent(in) :: var
integer(kind=2), dimension(il,6*il) :: ipack
integer(kind=2), dimension(il_l,jl_l,vnode_nproc) :: pf_ipack

jl=6*il
ifull=il*jl
      
!allocate( ipack(il,6*il) )

write(6,*)"histwrt3 sname=",sname," iarch=",iarch," idnc=",idnc

! find variable index
ier = nf_inq_varid(idnc,sname,mid)
ier = nf_get_att_real(idnc,mid,'add_offset',addoff)
ier = nf_get_att_real(idnc,mid,'scale_factor',scale_f)

xmin=addoff+scale_f*minv
!xmax=xmin+scale_f*float(maxv-minv)
xmax=xmin+scale_f*(real(maxv)-real(minv)) ! jlm fix for precision problems

varn= 1.e29
varx=-1.e29
do j=1,jl
  do i=1,il
    if(var(i,j).lt.varn)then
      varn=var(i,j)
      imn=i
      jmn=j
    endif
    if(var(i,j).gt.varx)then
      varx=var(i,j)
      imx=i
      jmx=j
    endif
    pvar = max(xmin,min(xmax,var(i,j))) ! limited output variable
    ipack(i,j)=nint((pvar-addoff)/scale_f)
    ipack(i,j)=max(min(ipack(i,j),maxv),minv)
  end do
end do

ier = nf_inq_varndims(idnc,mid,ndims)
      
if ( vnode_nproc>0 ) then
  ! reorder data for procfomat
  start(1) = 1
  start(2) = 1
  start(3) = 1
  start(4) = iarch
  count(1) = il_l
  count(2) = jl_l
  count(3) = vnode_nproc
  count(4) = 1
  do np = 0,vnode_nproc-1
    ioff = ipoff(np+1)
    joff = jpoff(np+1)
    noff = npoff(np+1)
    do n = 1,npan_l
      do j = 1,jpan
        do i = 1,ipan
          pf_ipack(i,j,np+1) = ipack(i+ioff,j+joff+(n-noff)*jpan)
        end do  
      end do
    end do  
  end do
  ier = nf_put_vara_int2(idnc, mid, start(1:ndims+1), count(1:ndims+1), pf_ipack)
  if(ier.ne.0) then
    write(6,*) "in histwrt3 ier not zero",ier,sname
    call finishbanner
    stop
  end if  
else
  start(1) = 1
  start(2) = 1
  start(3) = iarch
  count(1) = il
  count(2) = jl
  count(3) = 1  
  ier = nf_put_vara_int2(idnc, mid, start(1:ndims), count(1:ndims), ipack)
  if(ier.ne.0) then
    write(6,*) "in histwrt3 ier not zero",ier,sname
    call finishbanner
    stop
  end if  
end if

write(6,'("histwrt3:",a7," nt=",i4," n=",f12.4," ij=",2i4," x=",f12.4," ij=",2i4)') sname,iarch,varn,imn,jmn,varx,imx,jmx

!deallocate( ipack )
      
return
end subroutine histwrt3
!=======================================================================
subroutine histwrt4(var,sname,idnc,iarch,il,kl)
! Write 3d+t fields from the savegrid array.

use netcdf_m
      
integer il,jl,kl,ifull

!include 'newmpar.h'
!include 'parm.h'

integer mid, start(5), count(5)
character(len=*), intent(in) :: sname
!character*8 sname
integer*2 minv, maxv, missval ! was integer*2 
parameter(minv = -32500, maxv = 32500, missval = -32501)
real addoff, scale_f
integer np, ioff, joff, noff, n

real, dimension(il,6*il,kl), intent(in) :: var
integer(kind=2), dimension(il,6*il,kl) :: ipack
integer(kind=2), dimension(il_l,jl_l,kl,vnode_nproc) :: pf_ipack

jl=6*il
ifull=il*jl

!allocate( ipack(il,6*il,kl) )
      
write(6,*)"histwrt4 sname=",sname," iarch=",iarch," idnc=",idnc


! find variable index
ier = nf_inq_varid(idnc,sname,mid)
ier = nf_get_att_real(idnc,mid,'add_offset',addoff)
ier = nf_get_att_real(idnc,mid,'scale_factor',scale_f)

xmin=addoff+scale_f*minv
!xmax=xmin+scale_f*float(maxv-minv)
xmax=xmin+scale_f*(real(maxv)-real(minv)) ! jlm fix for precision problems

varn= 1.e29
varx=-1.e29
imx=0
jmx=0
kmx=0
imn=0
jmn=0
kmn=0
do k=1,kl
  do j=1,jl
    do i=1,il
      pvar = max(xmin,min(xmax,var(i,j,k))) ! limited output variable
      ipack(i,j,k)=nint((pvar-addoff)/scale_f)
      ipack(i,j,k)=max(min(ipack(i,j,k),maxv),minv)
      if(var(i,j,k).gt.varx)then
        varx=var(i,j,k)
        imx=i
        jmx=j
        kmx=k
      endif
      if(var(i,j,k).lt.varn)then
        varn=var(i,j,k)
        imn=i
        jmn=j
        kmn=k
      endif
    end do
  end do
end do

if ( vnode_nproc>0 ) then
  start(1) = 1
  start(2) = 1
  start(3) = 1
  start(4) = 1
  start(5) = iarch
  count(1) = il_l
  count(2) = jl_l
  count(3) = kl
  count(4) = vnode_nproc
  count(5) = 1
  do np = 0,vnode_nproc-1
    ioff = ipoff(np+1)
    joff = jpoff(np+1)
    noff = npoff(np+1)
    do k = 1,kl
      do n = 1,npan_l
        do j = 1,jpan
          do i = 1,ipan
            pf_ipack(i,j,k,np+1) = ipack(i+ioff,j+joff+(n-noff)*jpan,k)
          end do
        end do  
      end do
    end do  
  end do
  ier = nf_put_vara_int2(idnc, mid, start, count, pf_ipack)
else    
  start(1) = 1
  start(2) = 1
  start(3) = 1
  start(4) = iarch
  count(1) = il
  count(2) = jl
  count(3) = kl
  count(4) = 1
  ier = nf_put_vara_int2(idnc, mid, start(1:4), count(1:4), ipack)
end if  

write(6,'("histwrt4:",a7," nt=",i4," n=",f12.4," ijk=",3i4," x=",f12.4," ijk=",3i4)') &
    sname,iarch,varn,imn,jmn,kmn,varx,imx,jmx,kmx

!deallocate( ipack )
      
return
end subroutine histwrt4
     
subroutine mtimerget(mtimer,kdate1,ktime1,kdate2,ktime2) ! jlm
!     returns mtimer in minutes, corr. to (kdate2,ktime2) -  (kdate1,ktime1)    
dimension ndoy(12)   ! days from beginning of year (1st Jan is 0)
data ndoy/ 0,31,59,90,120,151,181,212,243,273,304,334/
common/leap_yr/leap  ! 1 to allow leap years
 
if(leap.ne.0) then
  write(6,*) 'leap years not catered for in mtimerget'
  call finishbanner
  stop
end if
!     Set up number of minutes from beginning of year
!     For GCM runs assume year is <1980 (e.g. ~321-460 for 140 year run)
jyear1=kdate1/10000
jmonth=(kdate1-jyear1*10000)/100
jday=kdate1-jyear1*10000-jmonth*100
jhour=ktime1/100
jmin=ktime1-jhour*100
mstart1=1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of y

jyear2=kdate2/10000
jmonth=(kdate2-jyear2*10000)/100
jday=kdate2-jyear2*10000-jmonth*100
jhour=ktime2/100
jmin=ktime2-jhour*100
mstart2=1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of y

mtimer=mstart2-mstart1+(jyear2-jyear1)*365*24*60
return
end subroutine mtimerget

end module outcdf_m
