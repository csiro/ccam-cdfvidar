! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2020 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
public readvar, readsst, readsoil
public netcdferror, readpress, datefix

integer ixp, iyp, idlev, idnt

interface readvar
  module procedure readvar3d, readvar2d, readvarinv
end interface

contains

subroutine outcdf(ihr,idy,imon,iyr,iout,nt,time,mtimer,sig,cdffile_in,ddss,il,kl, &
                  minlon,maxlon,minlat,maxlat,llrng,calendar,rlong0,rlat0,schmidt)

use netcdf_m
      
implicit none
      
integer il,jl,kl,ifull
integer ihr,idy,imon,iyr,iout
integer nt
real ddss

character(len=*), intent(in) :: calendar
character(len=10) rundate

integer, save :: idnc0, idnc1, idncm1
integer, save :: kdate, ktime
integer iarch
integer idnc, ier, imode
integer ktau, icy, icm, icd
integer ich, icmi, ics, mtimer, idv
integer, parameter :: nextout=0
integer, parameter :: itype=1
integer, parameter :: nihead=54
integer nahead(nihead)
integer, dimension(1) :: dimids
real, intent(in) :: minlon, maxlon, minlat, maxlat, llrng
real, intent(in) :: rlong0, rlat0, schmidt
real(kind=8), intent(inout) :: time
real(kind=8), save :: first_time

integer, parameter :: nrhead = 14
real ahead(nrhead)

real, dimension(kl) :: sig
real dt,ds

character(len=1024) cdffile_in
character(len=1032) cdffile

integer, dimension(4) :: dim
integer xdim,ydim,zdim,tdim
integer oldmode
character(len=20) :: timorg
character(len=33) :: grdtim
character(len=3), dimension(12) :: month
data month/'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'/

data idnc1/0/, idnc0/0/, idncm1/0/
data rundate/"ncepavnanl"/

logical, save :: first = .true.

jl=6*il
ifull=il*jl

nahead = 0
ahead = 0.
ddss = 0.

write(6,*)"outcdf ihr,idy,imon,iyr,iout=",ihr,idy,imon,iyr,iout
write(6,*)"time=",time

dt=0

if ( first ) then
  kdate=iyr*10000+imon*100+idy
  ktime=ihr*100
  write(6,*)"kdate,ktime=",kdate,ktime
  first_time = time
  call datefix(kdate,ktime,time,calendar)
  first = .false.
else  
  time = time - first_time  
end if  
write(6,*)"time=",time
write(6,*)"kdate,ktime=",kdate,ktime


! itype=1 outfile
iarch=nt
idnc=idnc1

write(6,'("outcdf itype,idnc,iarch,cdffile=",3i5," ",a80)') itype,idnc,iarch,cdffile_in
      
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
  cdffile = trim(cdffile_in)  
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
  ier = nf_def_dim(idnc,'longitude', il,           xdim)
  ier = nf_def_dim(idnc,'latitude',  jl,           ydim)
  ier = nf_def_dim(idnc,'lev',       kl,           zdim)
  ier = nf_def_dim(idnc,'time',      nf_unlimited, tdim)
  write(6,*) "xdim=",xdim," ydim=",ydim," zdim=",zdim," tdim=",tdim
          
  ! define coords.
          
  dimids = xdim
  ier = nf_def_var(idnc,'longitude',nf_float,1,dimids,ixp)
  ier = nf_put_att_text(idnc,ixp,'point_spacing',4,'even')
  ier = nf_put_att_text(idnc,ixp,'units',12,'degrees_east')
  dimids = ydim
  ier = nf_def_var(idnc,'latitude',nf_float,1,dimids,iyp)
  ier = nf_put_att_text(idnc,iyp,'point_spacing',4,'even')
  ier = nf_put_att_text(idnc,iyp,'units',13,'degrees_north')
  write(6,*)'ixp,iyp=',ixp,iyp

  dimids = zdim
  ier = nf_def_var(idnc,'lev',nf_float,1,dimids,idlev)
  write(6,*)'idlev,ier=',idlev,ier
  ier = nf_put_att_text(idnc,idlev,'positive',4,'down')
  ier = nf_put_att_text(idnc,idlev,'point_spacing',6,'uneven')
  ier = nf_put_att_text(idnc,idlev,'units',11,'sigma_level')
  ier = nf_put_att_text(idnc,idlev,'long_name',11,'sigma_level')
  write(6,*)'idlev=',idlev

  write(6,*)'tdim,idnc=',tdim,idnc
  dimids = tdim
  ier = nf_def_var(idnc,'time',nf_double,1,dimids,idnt)
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

  if ( calendar/="" ) then
    ier = nf_put_att_text(idnc,idnt,'calendar',len_trim(calendar),calendar)
  end if  
  
  dim(1) = xdim
  dim(2) = ydim
  dim(3) = zdim
  dim(4) = tdim
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
  ahead(5)=rlong0     ! needed by cc2hist
  ahead(6)=rlat0      ! needed by cc2hist
  ahead(7)=schmidt    ! needed by cc2hist
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
integer dim(4)
integer idim2(3)
integer, dimension(1) :: ivals
integer, dimension(1) :: start, ncount
integer, dimension(4) :: ostart, ocount
integer nrun
real, dimension(il) :: xpnt
real, dimension(6*il) :: ypnt
real, dimension(kl), intent(in) :: sig
real, intent(in) :: minlon, maxlon, minlat, maxlat, llrng
real(kind=8), intent(in) :: time

!common/cdfind/ixp,iyp,idlev,idnt
!       *** qscrn_ave not presently written     
real, dimension(6*il*il) :: aa
real, dimension(1) :: rvals
      
integer iq, varid
integer i, j, n

jl=6*il
ifull=il*jl

write(6,*)'openhist iarch,idnc,itype=',iarch,idnc,itype

!if this is the first archive, set up some global attributes
if(iarch.eq.1) then
  write(6,*)'dim=',dim(1:4)
  idim2(1)=dim(1)    ! x
  idim2(2)=dim(2)    ! y
  idim2(3)=dim(4)    ! time
  write(6,*)'idim2=',idim2

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
        
! Sigma levels
! write(6,*)'sig=',sig
  ier = nf_put_att_real(idnc,nf_global,'sigma',nf_float,kl,sig)

  !lname = 'year-month-day at start of run'
  !ier = nf_def_var(idnc,'kdate',nf_int,1,dim(4:4),idkdate)
  !ier = nf_put_att_text(idnc,idkdate,'long_name',len_trim(lname),lname)

  !lname = 'hour-minute at start of run'
  !ier = nf_def_var(idnc,'ktime',nf_int,1,dim(4:4),idktime)
  !ier = nf_put_att_text(idnc,idktime,'long_name',len_trim(lname),lname)

  !lname = 'timer (hrs)'
  !ier = nf_def_var(idnc,'timer',nf_float,1,dim(4:4),idnter)
  !ier = nf_put_att_text(idnc,idnter,'long_name',len_trim(lname),lname)

  !lname = 'mtimer (mins)'
  !ier = nf_def_var(idnc,'mtimer',nf_int,1,dim(4:4),idmtimer)
  !ier = nf_put_att_text(idnc,idmtimer,'long_name',len_trim(lname),lname)

  !lname = 'timeg (UTC)'
  !ier = nf_def_var(idnc,'timeg',nf_float,1,dim(4:4),idnteg)
  !ier = nf_put_att_text(idnc,idnteg,'long_name',len_trim(lname),lname)

  !lname = 'number of time steps from start'
  !ier = nf_def_var(idnc,'ktau',nf_int,1,dim(4:4),idktau)
  !ier = nf_put_att_text(idnc,idktau,'long_name',len_trim(lname),lname)

  ier = nf_def_var(idnc,'sigma',nf_float,1,dim(3:3),idv)
  ier = nf_put_att_text(idnc,idv,'positive',4,'down')

  write(6,*)'define attributes of variables'

  lname ='Scaled Log Surface pressure'
  call attrib(idnc,idim2,3,'psf',lname,'none',-1.3,0.2)

  lname ='Mean sea level pressure'
  call attrib(idnc,idim2,3,'pmsl',lname,'hPa',800.,1200.)
  lname = 'Surface geopotential'
  call attrib(idnc,idim2,2,'zht',lname,'m2/s2',-2.e3,128.e3) ! MJT lsmask
  lname = 'Soil type'                                               ! MJT lsmask
  call attrib(idnc,idim2,2,'soilt',lname,'none',0.,65.e3)    ! MJT lsmask

! For time invariant surface fields
  lname = 'clon'
  call attrib(idnc,idim2,2,'clon',lname,'none',-360.,360.)
  lname = 'clat'
  call attrib(idnc,idim2,2,'clat',lname,'none',-90.,90.)

! For time varying surface fields
  lname = 'Surface temperature'
  call attrib(idnc,idim2,3,'tsu',lname,'K',0.,350.)

  if ( all(soiltemp<0.) ) then
    lname = 'Soil temperature top'
    call attrib(idnc,idim2,3,'tb3',lname,'K',100.,400.)
    lname = 'Soil temperature bottom'
    call attrib(idnc,idim2,3,'tb2',lname,'K',100.,400.)
  else
    lname = 'Soil temperature lev 1'
    call attrib(idnc,idim2,3,'tgg1',lname,'K',100.,400.)
    lname = 'Soil temperature lev 2'
    call attrib(idnc,idim2,3,'tgg2',lname,'K',100.,400.)
    lname = 'Soil temperature lev 3'
    call attrib(idnc,idim2,3,'tgg3',lname,'K',100.,400.)
    lname = 'Soil temperature lev 4'
    call attrib(idnc,idim2,3,'tgg4',lname,'K',100.,400.)
    lname = 'Soil temperature lev 5'
    call attrib(idnc,idim2,3,'tgg5',lname,'K',100.,400.)
    lname = 'Soil temperature lev 6'
    call attrib(idnc,idim2,3,'tgg6',lname,'K',100.,400.)
  end if
      
  if ( all(soilmoist<0.) ) then
    lname = 'Soil moisture top'
    call attrib(idnc,idim2,3,'wfg',lname,'none',0.,.4)
    lname = 'Soil moisture bottom'
    call attrib(idnc,idim2,3,'wfb',lname,'none',0.,.4)
  else
    lname = 'Soil moisture 1'
    call attrib(idnc,idim2,3,'wb1',lname,'m3/m3',0.,2.)
    lname = 'Soil moisture 2'
    call attrib(idnc,idim2,3,'wb2',lname,'m3/m3',0.,2.)
     lname = 'Soil moisture 3'
    call attrib(idnc,idim2,3,'wb3',lname,'m3/m3',0.,2.)
    lname = 'Soil moisture 4'
    call attrib(idnc,idim2,3,'wb4',lname,'m3/m3',0.,2.)
    lname = 'Soil moisture 5'
    call attrib(idnc,idim2,3,'wb5',lname,'m3/m3',0.,2.)
    lname = 'Soil moisture 6'
    call attrib(idnc,idim2,3,'wb6',lname,'m3/m3',0.,2.)
  end if

  if (any(fracice>=0.)) then
    lname = 'Sea ice fraction'
    call attrib(idnc,idim2,3,'fracice',lname,'none',0.,6.5)
  end if
        
  if (any(snod>1.e-8)) then
    lname = 'Snow depth (liquid water)'
    call attrib(idnc,idim2,3,'snd',lname,'none',0.,6500.)
  end if

  write(6,*)'3d variables'
  call attrib(idnc,dim,4,'temp','Air temperature','K',100.,350.)
  call attrib(idnc,dim,4,'u','x-component wind','m/s',-150.,150.)
  call attrib(idnc,dim,4,'v','y-component wind','m/s',-150.,150.)
  lname= 'Water mixing ratio'
  call attrib(idnc,dim,4,'mixr',lname,'kg/kg',0.,.05)

  write(6,*)'finished defining attributes'
! Leave define mode
  ier = nf_enddef(idnc)
  write(6,*)'leave define mode: ier=',ier

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

  start = 1
  ncount = kl
  ier = nf_put_vara_real(idnc,idlev,start,ncount,sig)
  write(6,*) "sig=",sig

  ier = nf_inq_varid(idnc,'sigma',idv)
  start = 1
  ncount = kl
  ier = nf_put_vara_real(idnc,idv,start,ncount,sig)

endif ! iarch.eq.1
     
!------------------------------------------------------------------      
ktau=0
!timer=0
timer=real(time/60._8) ! MJT quick fix
!timeg=0
timeg=mod(real(time/60._8),24.) ! MJT quick fix
write(6,*)'outcdf processing kdate,ktime,ktau,time,mtimer: ',kdate,ktime,ktau,time,mtimer

!set time to number of minutes since start
ier = nf_inq_varid(idnc,'time',idv)
start = iarch
ier = nf_put_var1_double(idnc,idv,start,time)
write(6,*)"int(time)=",int(time,8)

!ier = nf_inq_varid(idnc,'kdate',idv)
!ier = nf_put_var1_int(idnc,idv,start,kdate)
!ier = nf_inq_varid(idnc,'ktime',idv)
!ier = nf_put_var1_int(idnc,idv,start,ktime)
!ier = nf_inq_varid(idnc,'timer',idv)
!ier = nf_put_var1_real(idnc,idv,start,timer)   
!ier = nf_inq_varid(idnc,'mtimer',idv)
!ier = nf_put_var1_int(idnc,idv,start,nint(time)) ! MJT quick fix
!ier = nf_inq_varid(idnc,'timeg',idv)
!ier = nf_put_var1_real(idnc,idv,start,timeg)
!ier = nf_inq_varid(idnc,'ktau',idv)
!ier = nf_put_var1_int(idnc,idv,start,ktau)      

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
      
return ! subroutine openhist(idnc,iarch,itype,dim,sig
end subroutine openhist
!=======================================================================
subroutine attrib(cdfid,dim,ndim,name,lname,units,xmin,xmax)

use netcdf_m

integer(kind=2) minv, maxv, missval   ! was integer*2
parameter(minv = -32500, maxv = 32500, missval = -32501)
integer ndim
integer cdfid, idv, dim(ndim)
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

integer mid, start(3), count(3)
integer imn, imx, jmn, jmx
integer i, j, ndims
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

jl=6*il
ifull=il*jl
      
!allocate( ipack(il,6*il) )

write(6,*)"histwrt3 sname=",sname," iarch=",iarch," idnc=",idnc

! find variable index
ier = nf_inq_varid(idnc,sname,mid)
ier = nf_get_att_real(idnc,mid,'add_offset',addoff)
ier = nf_get_att_real(idnc,mid,'scale_factor',scale_f)

xmin=addoff+scale_f*minv
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

integer mid, start(4), count(4)
character(len=*), intent(in) :: sname
!character*8 sname
integer*2 minv, maxv, missval ! was integer*2 
parameter(minv = -32500, maxv = 32500, missval = -32501)
real addoff, scale_f

real, dimension(il,6*il,kl), intent(in) :: var
integer(kind=2), dimension(il,6*il,kl) :: ipack

jl=6*il
ifull=il*jl

!allocate( ipack(il,6*il,kl) )
      
write(6,*)"histwrt4 sname=",sname," iarch=",iarch," idnc=",idnc


! find variable index
ier = nf_inq_varid(idnc,sname,mid)
ier = nf_get_att_real(idnc,mid,'add_offset',addoff)
ier = nf_get_att_real(idnc,mid,'scale_factor',scale_f)

xmin=addoff+scale_f*minv
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

start(1) = 1
start(2) = 1
start(3) = 1
start(4) = iarch
count(1) = il
count(2) = jl
count(3) = kl
count(4) = 1
ier = nf_put_vara_int2(idnc, mid, start(1:4), count(1:4), ipack)

write(6,'("histwrt4:",a7," nt=",i4," n=",f12.4," ijk=",3i4," x=",f12.4," ijk=",3i4)') &
    sname,iarch,varn,imn,jmn,kmn,varx,imx,jmx,kmx

!deallocate( ipack )
      
return
end subroutine histwrt4
     
! Updated version for reading file data

subroutine readvar3d(ncid,varname,kdate,ktime,iarchi,sdiag, &
                     ptype,in_plev,dataout,ier,units)

use netcdf_m

implicit none

integer, intent(in) :: kdate, ktime, iarchi
integer, intent(in) :: ncid
integer, intent(out) :: ier
integer idv, lonid, latid
integer ix, iy, il, khout, khin
integer idpres, ivpres, ivtim, nplev, idvar
integer in_iyr, in_imn, in_idy, in_ihr, in_imi
integer itype
integer in_kdate, in_ktime, i
integer, dimension(4) :: start, ncount
integer, dimension(:), allocatable :: ivar
real, dimension(:,:,:), intent(out) :: dataout
real, dimension(:), intent(in) :: in_plev
real, dimension(:), allocatable :: glon, glat
real, dimension(:), allocatable :: datan
real, dimension(size(in_plev)) :: plev, plev_b
real fill_float, addoff, sf
real(kind=8), dimension(:), allocatable :: dvar
real(kind=8) in_time
logical orev, osig_in
logical, intent(in) :: sdiag
character(len=*), intent(in) :: varname
character(len=*), intent(in) :: ptype
character(len=*), intent(inout), optional :: units
character(len=1) in_type
character(len=60) timorg, cu
character(len=20) in_calendar

ier = nf_inq_varid(ncid,varname,idvar)  
if ( ier/=nf_noerr ) then
  return 
end if

if ( present(units) ) then
  ier = nf_get_att_text(ncid,idvar,'units',units)  
end if

! check date
ier = nf_inq_varid(ncid,'time',ivtim)
call netcdferror(ier)
ier = nf_get_att_text(ncid,ivtim,'units',timorg)
call netcdferror(ier)
ier = nf_get_var1_double(ncid,ivtim,iarchi,in_time)
call netcdferror(ier)
in_calendar=""
ier = nf_get_att_text(ncid,ivtim,'calendar',in_calendar)
call processdatestring(timorg,in_iyr,in_imn,in_idy,in_ihr,in_imi)
in_kdate = in_iyr*10000 + in_imn*100 + in_idy
in_ktime = in_ihr*100 + in_imi
i=scan(timorg,' ')-1
cu=''  ! clear string to ensure blank
cu(1:i)=timorg(1:i)
if ( cu(1:i) == "since" ) then
  cu="hours"
endif
select case(cu) ! MJT quick fix
  case('days')
    in_time=in_time*1440._8 
  case('hours')
    in_time=in_time*60._8 
  case('minutes')
    ! no change	
  case('seconds')
    in_time=in_time/60._8
  case DEFAULT
    write(6,*) "cannot convert unknown time unit ",trim(cu)
    call finishbanner
    stop -1
end select
call datefix(in_kdate,in_ktime,in_time,in_calendar)
if ( kdate/=in_kdate .or. in_ktime/=ktime ) then
  write(6,*) "ERROR: Inconsistent time units with ",trim(varname)
  call finishbanner
  stop -1
end if

! read lat/lon
ier = nf_inq_dimid(ncid,'lon',lonid)
if ( ier==nf_noerr ) then
  ier = nf_inq_varid(ncid,'lon',idv)
else
  ier = nf_inq_dimid(ncid,'longitude',lonid)
  call netcdferror(ier)
  ier = nf_inq_varid(ncid,'longitude',idv)
end if
ier = nf_inq_dimlen(ncid,lonid,ix)
allocate( glon(ix) )
ier = nf_get_var_real(ncid,idv,glon)

ier = nf_inq_dimid(ncid,'lat',latid)
if ( ier==nf_noerr ) then
  ier = nf_inq_varid(ncid,'lat',idv)    
else
  ier = nf_inq_dimid(ncid,'latitude',latid)
  call netcdferror(ier)
  ier = nf_inq_varid(ncid,'latitude',idv)    
end if
ier = nf_inq_dimlen(ncid,latid,iy)
allocate( glat(iy) )
ier = nf_get_var_real(ncid,idv,glat)

! read lev
call readpress(ncid,in_type,plev,plev_b,nplev,osig_in,orev)
if ( nplev/=size(dataout,3) ) then
  write(6,*) "ERROR: Inconsistent number of vertical levels with ",trim(varname)
  call finishbanner
  stop -1
end if
if ( in_type/=ptype ) then
  write(6,*) "ERROR: Inconsistent type of vertical levels with ",trim(varname)
  call finishbanner
  stop -1
end if
if ( any( abs(plev(1:nplev)-in_plev)>1.e-4 ) ) then
  write(6,*) "ERROR: Inconsistent values for vertical levels with ",trim(varname)
  write(6,*) "plev ",plev
  write(6,*) "in_plev ",in_plev
  call finishbanner
  stop -1
end if

! read data
allocate( datan(ix*iy*nplev) )
start(1) = 1
start(2) = 1
start(3) = 1
start(4) = iarchi
ncount(1) = ix
ncount(2) = iy
ncount(3) = nplev
ncount(4) = 1
ier = nf_get_att_real(ncid,idvar,'add_offset',addoff)
if ( ier/=nf_noerr ) addoff = 0.
ier = nf_get_att_real(ncid,idvar,'scale_factor',sf)
if ( ier/=nf_noerr ) sf = 1.
!ier = nf_inq_vartype(ncid,idvar,itype)
!call netcdferror(ier)
!select case(itype)
!  case ( nf_short )
!    allocate( ivar(ix*iy*nplev) )
!    ier = nf_get_vara_int(ncid,idvar,start,ncount,ivar)
!    call netcdferror(ier)
!    datan(1:ix*iy*nplev) = sf*real(ivar(1:ix*iy*nplev)) + addoff
!    deallocate( ivar )
!  case ( nf_float )
    ier = nf_get_vara_real(ncid,idvar,start,ncount,datan)
    call netcdferror(ier)
    datan(1:ix*iy*nplev) = sf*real(datan(1:ix*iy*nplev)) + addoff      
!  case ( nf_double )
!    allocate( dvar(ix*iy*nplev) )
!    ier = nf_get_vara_double(ncid,idvar,start,ncount,dvar)
!    call netcdferror(ier)
!    datan(1:ix*iy*nplev) = sf*real(dvar(1:ix*iy*nplev)) + addoff
!    deallocate( dvar )      
!  case default
!    write(6,*) "Variable is unkown"
!    call finishbanner
!    stop -1
!end select

ier = nf_get_att_real(ncid,idvar,'_FillValue',fill_float)
if ( ier/=nf_noerr ) then
  ier = nf_get_att_real(ncid,idvar,'missing_value',fill_float)    
end if
if ( ier==nf_noerr ) then
  where ( datan==fill_float )
    datan = 1.e10
  end where
  !call getvalidlev(validlevhost,datan,ix,iy,nplev)
  call filldat(datan,ix,iy,nplev)
end if
ier =  nf_noerr ! reset error even if fillvalue is not located
  
! interpolate data to cc grid
il = size(dataout,1)
call amap(datan(1:ix*iy),ix,iy,varname,0.,0.)
call amap(datan(1+ix*iy*(nplev-1):ix*iy*nplev),ix,iy,varname,0.,0.)
! vertical pressure levels are reversed just before cdfvidar
if ( orev ) then
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(nplev,datan,ix,iy,dataout,glon,glat,sdiag,il) PRIVATE(khin,khout)
  do khin = 1,nplev    
    khout = khin  
    call sintp16(datan(1+ix*iy*(khin-1):ix*iy*khin),ix,iy,dataout(:,:,khout),glon,glat,sdiag,il)  
  end do
!$OMP END PARALLEL DO 
else
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(nplev,datan,ix,iy,dataout,glon,glat,sdiag,il) PRIVATE(khin,khout)
  do khin = 1,nplev    
    khout = nplev + 1 - khin  
    call sintp16(datan(1+ix*iy*(khin-1):ix*iy*khin),ix,iy,dataout(:,:,khout),glon,glat,sdiag,il)  
  end do
!$OMP END PARALLEL DO  
end if

deallocate( glon, glat )
deallocate( datan )

return
end subroutine readvar3d

subroutine readvar2d(ncid,varname,kdate,ktime,iarchi,sdiag, &
                     dataout,ier,units)

use netcdf_m

implicit none

integer, intent(in) :: kdate, ktime, iarchi, ncid
integer, intent(out) :: ier
integer, dimension(:), allocatable :: ivar
integer, dimension(3) :: start, ncount
integer ix, iy, lonid, latid
integer idv, idvar, itype, il, ivtim
integer in_iyr, in_imn, in_idy, in_ihr, in_imi
integer in_kdate, in_ktime, i
real, dimension(:,:), intent(out) :: dataout
real, dimension(:), allocatable, save :: glon, glat
real, dimension(:), allocatable :: datan
real sf, addoff, fill_float
real(kind=8), dimension(:), allocatable :: dvar
real(kind=8) in_time
logical, intent(in) :: sdiag
character(len=*), intent(in) :: varname
character(len=*), intent(inout), optional :: units
character(len=60) timorg, cu
character(len=20) in_calendar

ier = nf_inq_varid(ncid,varname,idvar)  
if ( ier/=nf_noerr ) then
  return 
end if

if ( present(units) ) then
  ier = nf_get_att_text(ncid,idvar,'units',units) 
end if

! check date
ier = nf_inq_varid(ncid,'time',ivtim)
call netcdferror(ier)
ier = nf_get_att_text(ncid,ivtim,'units',timorg)
call netcdferror(ier)
ier = nf_get_var1_double(ncid,ivtim,iarchi,in_time)
call netcdferror(ier)
in_calendar=""
ier = nf_get_att_text(ncid,ivtim,'calendar',in_calendar)
call processdatestring(timorg,in_iyr,in_imn,in_idy,in_ihr,in_imi)
in_kdate = in_iyr*10000 + in_imn*100 + in_idy
in_ktime = in_ihr*100 + in_imi
i=scan(timorg,' ')-1
cu=''  ! clear string to ensure blank
cu(1:i)=timorg(1:i)
if ( cu(1:i) == "since" ) then
  cu="hours"
endif
select case(cu) ! MJT quick fix
  case('days')
    in_time=in_time*1440._8 
  case('hours')
    in_time=in_time*60._8 
  case('minutes')
    ! no change	
  case('seconds')
    in_time=in_time/60._8
  case DEFAULT
    write(6,*) "cannot convert unknown time unit ",trim(cu)
    call finishbanner
    stop -1
end select
call datefix(in_kdate,in_ktime,in_time,in_calendar)
if ( kdate/=in_kdate .or. in_ktime/=ktime ) then
  write(6,*) "ERROR: Inconsistent time units with ",trim(varname)
  call finishbanner
  stop -1
end if

! read lat/lon
ier = nf_inq_dimid(ncid,'lon',lonid)
if ( ier==nf_noerr ) then
  ier = nf_inq_varid(ncid,'lon',idv)
else
  ier = nf_inq_dimid(ncid,'longitude',lonid)
  call netcdferror(ier)
  ier = nf_inq_varid(ncid,'longitude',idv)
end if
ier = nf_inq_dimlen(ncid,lonid,ix)
allocate( glon(ix) )
ier = nf_get_var_real(ncid,idv,glon)

ier = nf_inq_dimid(ncid,'lat',latid)
if ( ier==nf_noerr ) then
  ier = nf_inq_varid(ncid,'lat',idv)    
else
  ier = nf_inq_dimid(ncid,'latitude',latid)
  call netcdferror(ier)
  ier = nf_inq_varid(ncid,'latitude',idv)    
end if
ier = nf_inq_dimlen(ncid,latid,iy)
allocate( glat(iy) )
ier = nf_get_var_real(ncid,idv,glat)

! read data
allocate( datan(ix*iy) )
start(1) = 1
start(2) = 1
start(3) = iarchi
ncount(1) = ix
ncount(2) = iy
ncount(3) = 1
ier = nf_get_att_real(ncid,idvar,'add_offset',addoff)
if ( ier/=nf_noerr ) addoff = 0.
ier = nf_get_att_real(ncid,idvar,'scale_factor',sf)
if ( ier/=nf_noerr ) sf = 1.
!ier = nf_inq_vartype(ncid,idvar,itype)
!call netcdferror(ier)
!select case(itype)
!  case ( nf_short )
!    allocate( ivar(ix*iy) )
!    ier = nf_get_vara_int(ncid,idvar,start,ncount,ivar)
!    call netcdferror(ier)
!    datan(1:ix*iy) = sf*real(ivar(1:ix*iy)) + addoff
!    deallocate( ivar )
!  case ( nf_float )
    ier = nf_get_vara_real(ncid,idvar,start,ncount,datan)
    call netcdferror(ier)
    datan(1:ix*iy) = sf*real(datan(1:ix*iy)) + addoff      
!  case ( nf_double )
!    allocate( dvar(ix*iy) )
!    ier = nf_get_vara_double(ncid,idvar,start,ncount,dvar)
!    call netcdferror(ier)
!    datan(1:ix*iy) = sf*real(dvar(1:ix*iy)) + addoff
!    deallocate( dvar )      
!  case default
!    write(6,*) "Variable is unkown"
!    call finishbanner
!    stop -1
!end select
  
ier = nf_get_att_real(ncid,idvar,'_FillValue',fill_float)
if ( ier/=nf_noerr ) then
  ier = nf_get_att_real(ncid,idvar,'missing_value',fill_float)    
end if
if ( ier==nf_noerr ) then
  where ( datan==fill_float )
    datan = 1.e10
  end where
  call filldat(datan,ix,iy,1)
end if
ier =  nf_noerr ! reset error even if fillvalue is not located  
  
! interpolation
il = size(dataout,1)
call amap( datan(1:ix*iy), ix, iy, varname, 0., 0. )
call sintp16(datan(1:ix*iy),ix,iy,dataout,glon,glat,sdiag,il)

deallocate( glon, glat )
deallocate( datan )

return
end subroutine readvar2d

subroutine readvarinv(ncid,varname,sdiag,dataout,ier,units,island)

use netcdf_m

implicit none

integer, intent(in) :: ncid
integer, intent(out) :: ier
integer, dimension(:), allocatable :: ivar
integer, dimension(3) :: start, ncount
integer ix, iy, lonid, latid
integer idv, idvar, itype, il, ivtim
real, dimension(:,:), intent(out) :: dataout
real, dimension(:), allocatable, save :: glon, glat
real, dimension(:), allocatable :: datan
real sf, addoff
real(kind=8), dimension(:), allocatable :: dvar
logical, intent(in) :: sdiag
logical, intent(in), optional :: island
logical is_land
character(len=*), intent(in) :: varname
character(len=*), intent(inout), optional :: units

is_land = .false.
if ( present(island) ) then
  is_land = island
end if

ier = nf_inq_varid(ncid,varname,idvar)  
if ( ier/=nf_noerr ) then
  return 
end if

if ( present(units) ) then
  ier = nf_get_att_text(ncid,idvar,'units',units) 
end if

! read lat/lon
ier = nf_inq_dimid(ncid,'lon',lonid)
if ( ier==nf_noerr ) then
  ier = nf_inq_varid(ncid,'lon',idv)
else
  ier = nf_inq_dimid(ncid,'longitude',lonid)
  call netcdferror(ier)
  ier = nf_inq_varid(ncid,'longitude',idv)
end if
ier = nf_inq_dimlen(ncid,lonid,ix)
allocate( glon(ix) )
ier = nf_get_var_real(ncid,idv,glon)

ier = nf_inq_dimid(ncid,'lat',latid)
if ( ier==nf_noerr ) then
  ier = nf_inq_varid(ncid,'lat',idv)    
else
  ier = nf_inq_dimid(ncid,'latitude',latid)
  call netcdferror(ier)
  ier = nf_inq_varid(ncid,'latitude',idv)    
end if
ier = nf_inq_dimlen(ncid,latid,iy)
allocate( glat(iy) )
ier = nf_get_var_real(ncid,idv,glat)

! read data
allocate( datan(ix*iy) )
start(1) = 1
start(2) = 1
start(3) = 1
ncount(1) = ix
ncount(2) = iy
ncount(3) = 1
ier = nf_get_att_real(ncid,idvar,'add_offset',addoff)
if ( ier/=nf_noerr ) addoff = 0.
ier = nf_get_att_real(ncid,idvar,'scale_factor',sf)
if ( ier/=nf_noerr ) sf = 1.
!ier = nf_inq_vartype(ncid,idvar,itype)
!call netcdferror(ier)
!select case(itype)
!  case ( nf_short )
!    allocate( ivar(ix*iy) )
!    ier = nf_get_vara_int(ncid,idvar,start,ncount,ivar)
!    call netcdferror(ier)
!    datan(1:ix*iy) = sf*real(ivar(1:ix*iy)) + addoff
!    deallocate( ivar )
!  case ( nf_float )
    ier = nf_get_vara_real(ncid,idvar,start,ncount,datan)
    call netcdferror(ier)
    datan(1:ix*iy) = sf*real(datan(1:ix*iy)) + addoff      
!  case ( nf_double )
!    allocate( dvar(ix*iy) )
!    ier = nf_get_vara_double(ncid,idvar,start,ncount,dvar)
!    call netcdferror(ier)
!    datan(1:ix*iy) = sf*real(dvar(1:ix*iy)) + addoff
!    deallocate( dvar )      
!  case default
!    write(6,*) "Variable is unkown"
!    call finishbanner
!    stop -1
!end select
  
! interpolation
il = size(dataout,1)
if ( is_land ) then
   where ( abs(datan)>=1.e10 ) ! quick fix for land_mask
      datan = 0.
   end where
end if   
call amap( datan(1:ix*iy), ix, iy, varname, 0., 0. )
call sintp16(datan(1:ix*iy),ix,iy,dataout,glon,glat,sdiag,il)

deallocate( glon, glat )
deallocate( datan )

return
end subroutine readvarinv

subroutine readsoil(ncid,varname,kdate,ktime,iarchi,sdiag, &
                    soildepth_ccam,dataout,ier,units)

use netcdf_m

implicit none

integer, intent(in) :: kdate, ktime, iarchi
integer, intent(in) :: ncid
integer, intent(out) :: ier
integer lonid, latid, idvar, idv
integer ix, iy, il, khout, khin
integer idsoillvl, ivsoillvl, ivtim, new_nsoillvl
integer k, ksearch, ktest
integer in_iyr, in_imn, in_idy, in_ihr, in_imi
integer itype
integer in_kdate, in_ktime, i
integer, save :: nsoillvl
integer, dimension(4) :: start, ncount
integer, dimension(:), allocatable :: ivar
real, dimension(:), intent(in) :: soildepth_ccam
real, dimension(:,:,:), intent(out) :: dataout
real, dimension(:), allocatable :: glon, glat
real, dimension(:), allocatable :: datan, datatemp
real, dimension(:), allocatable, save :: soildepth_in
real fill_float, xfrac
real addoff, sf
real(kind=8), dimension(:), allocatable :: dvar
real(kind=8) in_time
logical, intent(in) :: sdiag
character(len=*), intent(in) :: varname
character(len=*), intent(inout), optional :: units
character(len=60) timorg, cu
character(len=20) in_calendar

ier = nf_inq_varid(ncid,varname,idvar)  
if ( ier/=nf_noerr ) then
  return 
end if

if ( present(units) ) then
  ier = nf_get_att_text(ncid,idvar,'units',units)  
end if

! check date
ier = nf_inq_varid(ncid,'time',ivtim)
call netcdferror(ier)
ier = nf_get_att_text(ncid,ivtim,'units',timorg)
call netcdferror(ier)
ier = nf_get_var1_double(ncid,ivtim,iarchi,in_time)
call netcdferror(ier)
in_calendar=""
ier = nf_get_att_text(ncid,ivtim,'calendar',in_calendar)
call processdatestring(timorg,in_iyr,in_imn,in_idy,in_ihr,in_imi)
in_kdate = in_iyr*10000 + in_imn*100 + in_idy
in_ktime = in_ihr*100 + in_imi
i=scan(timorg,' ')-1
cu=''  ! clear string to ensure blank
cu(1:i)=timorg(1:i)
if ( cu(1:i) == "since" ) then
  cu="hours"
endif
select case(cu) ! MJT quick fix
  case('days')
    in_time=in_time*1440._8 
  case('hours')
    in_time=in_time*60._8 
  case('minutes')
    ! no change	
  case('seconds')
    in_time=in_time/60._8
  case DEFAULT
    write(6,*) "cannot convert unknown time unit ",trim(cu)
    call finishbanner
    stop -1
end select
call datefix(in_kdate,in_ktime,in_time,in_calendar)
if ( kdate/=in_kdate .or. in_ktime/=ktime ) then
  write(6,*) "ERROR: Inconsistent time units with ",trim(varname)
  call finishbanner
  stop -1
end if

! read lat/lon
ier = nf_inq_dimid(ncid,'lon',lonid)
if ( ier==nf_noerr ) then
  ier = nf_inq_varid(ncid,'lon',idv)
else
  ier = nf_inq_dimid(ncid,'longitude',lonid)
  call netcdferror(ier)
  ier = nf_inq_varid(ncid,'longitude',idv)
end if
ier = nf_inq_dimlen(ncid,lonid,ix)
allocate( glon(ix) )
ier = nf_get_var_real(ncid,idv,glon)

ier = nf_inq_dimid(ncid,'lat',latid)
if ( ier==nf_noerr ) then
  ier = nf_inq_varid(ncid,'lat',idv)    
else
  ier = nf_inq_dimid(ncid,'latitude',latid)
  call netcdferror(ier)
  ier = nf_inq_varid(ncid,'latitude',idv)    
end if
ier = nf_inq_dimlen(ncid,latid,iy)
allocate( glat(iy) )
ier = nf_get_var_real(ncid,idv,glat)

! read lev
ier = nf_inq_dimid(ncid,'soil_lvl',idsoillvl)
if ( ier/=nf_noerr ) then
  write(6,*) "ERROR: Cannot real dimension soil_lvl when reading ",trim(varname)
  call finishbanner
  stop -1 
end if
ier = nf_inq_dimlen(ncid,idsoillvl,new_nsoillvl)
if ( .not.allocated(soildepth_in) ) then
  nsoillvl = new_nsoillvl
  ier = nf_inq_varid(ncid,'soil_lvl',ivsoillvl)
  if ( ier==nf_noerr ) then
    ier = nf_inq_dimid(ncid,'soil_lvl',idsoillvl)
    ier = nf_inq_dimlen(ncid,idsoillvl,nsoillvl)
    allocate( soildepth_in(nsoillvl) )
    ier = nf_get_var_real(ncid,ivsoillvl,soildepth_in)
  end if   
else
  if ( new_nsoillvl/=nsoillvl ) then
    write(6,*) "ERROR: Change in number of soil levels"
    call finishbanner
    stop -1
  end if
end if

! read data
allocate( datan(ix*iy*nsoillvl) )
start(1) = 1
start(2) = 1
start(3) = 1
start(4) = iarchi
ncount(1) = ix
ncount(2) = iy
ncount(3) = nsoillvl
ncount(4) = 1
ier = nf_get_att_real(ncid,idvar,'add_offset',addoff)
if ( ier/=nf_noerr ) addoff = 0.
ier = nf_get_att_real(ncid,idvar,'scale_factor',sf)
if ( ier/=nf_noerr ) sf = 1.
!ier = nf_inq_vartype(ncid,idvar,itype)
!call netcdferror(ier)
!select case(itype)
!  case ( nf_short )
!    allocate( ivar(ix*iy*nsoillvl) )
!    ier = nf_get_vara_int(ncid,idvar,start,ncount,ivar)
!    call netcdferror(ier)
!    datan(1:ix*iy*nsoillvl) = sf*real(ivar(1:ix*iy*nsoillvl)) + addoff
!    deallocate( ivar )
!  case ( nf_float )
    ier = nf_get_vara_real(ncid,idvar,start,ncount,datan)
    call netcdferror(ier)
    datan(1:ix*iy*nsoillvl) = sf*real(datan(1:ix*iy*nsoillvl)) + addoff      
!  case ( nf_double )
!    allocate( dvar(ix*iy*nsoillvl) )
!    ier = nf_get_vara_double(ncid,idvar,start,ncount,dvar)
!    call netcdferror(ier)
!    datan(1:ix*iy*nsoillvl) = sf*real(dvar(1:ix*iy*nsoillvl)) + addoff
!    deallocate( dvar )      
!  case default
!    write(6,*) "Variable is unkown"
!    call finishbanner
!    stop -1
!end select

ier = nf_get_att_real(ncid,idvar,'_FillValue',fill_float)
if ( ier/=nf_noerr ) then
  ier = nf_get_att_real(ncid,idvar,'missing_value',fill_float)    
end if
if ( ier==nf_noerr ) then
  where ( datan==fill_float )
    datan = 1.e10
  end where
  call filldat(datan,ix,iy,nsoillvl)
end if
ier = nf_noerr

! interpolate to CCAM soil levels
il = size(dataout,1)
allocate( datatemp(ix*iy) )
ksearch = 2
do k = 1,6
  write(6,*)"************************************************k=",k  
  do ktest = ksearch,nsoillvl
    if ( soildepth_in(ktest)>soildepth_ccam(k) ) exit
  end do
  ksearch = max( min( ktest, nsoillvl ), 2 )
  if ( soildepth_ccam(k)<soildepth_in(ksearch-1) ) then
    ! extrapolate
    datatemp(:) = datan(1:ix*iy) ! 1st level
  else if ( soildepth_ccam(k)>soildepth_in(ksearch) ) then
    ! extrapolate
    datatemp(:) = datan((nsoillvl-1)*ix*iy+1:nsoillvl*ix*iy) ! last level
  else
    ! interpolate
    xfrac = ( soildepth_ccam(k) - soildepth_in(ksearch-1) ) / ( soildepth_in(ksearch) - soildepth_in(ksearch-1) )
    datatemp(:) = (1.-xfrac)*datan((ksearch-2)*ix*iy+1:(ksearch-1)*ix*iy) &
                +      xfrac*datan((ksearch-1)*ix*iy+1:ksearch*ix*iy)
  end if
  call sintp16(datatemp(:),ix,iy,dataout(:,:,k),glon,glat,sdiag,il)
end do
deallocate( datatemp )

deallocate( glon, glat )
deallocate( datan )

return
end subroutine readsoil

subroutine readsst(ncid,lsm_ncid,varname,kdate,ktime,iarchi,sdiag, &
                   lsm_m,sfct,sstmode,ier,units)

use netcdf_m

implicit none

integer, intent(in) :: kdate, ktime, iarchi, sstmode, ncid, lsm_ncid
integer, intent(out) :: ier
integer, dimension(:), allocatable :: ivar
integer, dimension(3) :: start, ncount
integer ix, iy, lonid, latid
integer idv, idvar, itype, il, ivtim
integer in_iyr, in_imn, in_idy, in_ihr, in_imi
integer nlpnts, nopnts
integer i, j, iq
integer in_kdate, in_ktime
real, dimension(:), intent(in) :: lsm_m
real, dimension(:,:), intent(out) :: sfct
real, dimension(size(sfct,1),size(sfct,2)) :: sfcto_m
real, dimension(:), allocatable, save :: glon, glat, glon_lsm, glat_lsm
real, dimension(:), allocatable :: datan, datan_tmp, datan_ocn
real, dimension(:), allocatable, save :: lsm_gbl
real sf, addoff, spval
real fill_float
real(kind=8), dimension(:), allocatable :: dvar
real(kind=8) in_time
logical, save :: olsm_gbl
logical, intent(in) :: sdiag
character(len=*), intent(in) :: varname
character(len=*), intent(inout), optional :: units
character(len=60) timorg, cu
character(len=20) in_calendar


ier = nf_inq_varid(ncid,varname,idvar)
if ( ier/=nf_noerr ) then
  return 
end if

if ( present(units) ) then
  ier = nf_get_att_text(ncid,idvar,'units',units)  
end if

! check date
ier = nf_inq_varid(ncid,'time',ivtim)
call netcdferror(ier)
ier = nf_get_att_text(ncid,ivtim,'units',timorg)
call netcdferror(ier)
ier = nf_get_var1_double(ncid,ivtim,iarchi,in_time)
call netcdferror(ier)
in_calendar=""
ier = nf_get_att_text(ncid,ivtim,'calendar',in_calendar)
call processdatestring(timorg,in_iyr,in_imn,in_idy,in_ihr,in_imi)
in_kdate = in_iyr*10000 + in_imn*100 + in_idy
in_ktime = in_ihr*100 + in_imi
i=scan(timorg,' ')-1
cu=''  ! clear string to ensure blank
cu(1:i)=timorg(1:i)
if ( cu(1:i) == "since" ) then
  cu="hours"
endif
select case(cu) ! MJT quick fix
  case('days')
    in_time=in_time*1440._8 
  case('hours')
    in_time=in_time*60._8 
  case('minutes')
    ! no change	
  case('seconds')
    in_time=in_time/60._8
  case DEFAULT
    write(6,*) "cannot convert unknown time unit ",trim(cu)
    call finishbanner
    stop -1
end select
call datefix(in_kdate,in_ktime,in_time,in_calendar)
if ( kdate/=in_kdate .or. in_ktime/=ktime ) then
  write(6,*) "ERROR: Inconsistent time units with ",trim(varname)
  call finishbanner
  stop -1
end if

! read lat/lon
ier = nf_inq_dimid(ncid,'lon',lonid)
if ( ier==nf_noerr ) then
  ier = nf_inq_varid(ncid,'lon',idv)
else
  ier = nf_inq_dimid(ncid,'longitude',lonid)
  call netcdferror(ier)
  ier = nf_inq_varid(ncid,'longitude',idv)
end if
ier = nf_inq_dimlen(ncid,lonid,ix)
allocate( glon(ix) )
ier = nf_get_var_real(ncid,idv,glon)

ier = nf_inq_dimid(ncid,'lat',latid)
if ( ier==nf_noerr ) then
  ier = nf_inq_varid(ncid,'lat',idv)    
else
  ier = nf_inq_dimid(ncid,'latitude',latid)
  call netcdferror(ier)
  ier = nf_inq_varid(ncid,'latitude',idv)    
end if
ier = nf_inq_dimlen(ncid,latid,iy)
allocate( glat(iy) )
ier = nf_get_var_real(ncid,idv,glat)

! read data
write(6,*)"input data has sfc data, now read in"
allocate( datan(ix*iy), datan_ocn(ix*iy), datan_tmp(ix*iy) )
start(1) = 1
start(2) = 1
start(3) = iarchi
ncount(1) = ix
ncount(2) = iy
ncount(3) = 1
ier = nf_get_att_real(ncid,idvar,'add_offset',addoff)
if ( ier/=nf_noerr ) addoff = 0.
ier = nf_get_att_real(ncid,idvar,'scale_factor',sf)
if ( ier/=nf_noerr ) sf = 1.
!ier = nf_inq_vartype(ncid,idvar,itype)
!call netcdferror(ier)
!select case(itype)
!  case ( nf_short )
!    allocate( ivar(ix*iy) )
!    ier = nf_get_vara_int(ncid,idvar,start,ncount,ivar)
!    call netcdferror(ier)
!    datan(1:ix*iy) = sf*real(ivar(1:ix*iy)) + addoff
!    deallocate( ivar )
!  case ( nf_float )
    ier = nf_get_vara_real(ncid,idvar,start,ncount,datan)
    call netcdferror(ier)
    datan(1:ix*iy) = sf*real(datan(1:ix*iy)) + addoff      
!  case ( nf_double )
!    allocate( dvar(ix*iy) )
!    ier = nf_get_vara_double(ncid,idvar,start,ncount,dvar)
!    call netcdferror(ier)
!    datan(1:ix*iy) = sf*real(dvar(1:ix*iy)) + addoff
!    deallocate( dvar )      
!  case default
!    write(6,*) "Variable is unkown"
!    call finishbanner
!    stop -1
!end select

spval = -1.e10
if ( sstmode==0 ) then
  write(6,*)"spval=",spval
  if (any(datan(1:ix*iy)<100..or.datan(1:ix*iy)>400.)) then
    write(6,*) "Missing data found in sfc temp"
    where (datan(1:ix*iy)<100..or.datan(1:ix*iy)>400.)
      datan(1:ix*iy)=spval
    end where
    call fill(datan(1:ix*iy),ix,iy,.1*spval,datan_tmp)
  end if
end if  

ier = nf_get_att_real(ncid,idvar,'_FillValue',fill_float)
if ( ier/=nf_noerr ) then
  ier = nf_get_att_real(ncid,idvar,'missing_value',fill_float)    
end if
if ( ier==nf_noerr ) then
  where ( datan==fill_float )
    datan = spval
  end where
  call fill(datan(1:ix*iy),ix,iy,.1*spval,datan_tmp)
end if
ier =  nf_noerr ! reset error even if fillvalue is not located

! read land-sea mask
if ( .not.allocated(lsm_gbl) ) then
  ier = nf_inq_varid(lsm_ncid,"land",idvar) 
  if ( ier/=nf_noerr ) then
    ier = nf_inq_varid(lsm_ncid,"lsm",idvar)  
  end if
  if ( ier/=nf_noerr ) then
    ier = nf_inq_varid(lsm_ncid,"sftlf",idvar)  
  end if
  if ( ier/=nf_noerr ) then
    ier = nf_inq_varid(lsm_ncid,"sfc_lsm",idvar)  
  end if
  if ( ier/=nf_noerr ) then
    ier = nf_inq_varid(lsm_ncid,"land_mask",idvar)  
  end if
  if ( ier/=nf_noerr ) then
    ier = nf_inq_varid(lsm_ncid,"lnd_mask",idvar)  
  end if
  if ( ier==nf_noerr ) then
    olsm_gbl = .true.
    ! read lat/lon
    ier = nf_inq_dimid(lsm_ncid,'lon',lonid)
    if ( ier==nf_noerr ) then
      ier = nf_inq_varid(lsm_ncid,'lon',idv)
    else
      ier = nf_inq_dimid(lsm_ncid,'longitude',lonid)
      call netcdferror(ier)
      ier = nf_inq_varid(lsm_ncid,'longitude',idv)
    end if
    ier = nf_inq_dimlen(lsm_ncid,lonid,ix)
    allocate( glon_lsm(ix) )
    ier = nf_get_var_real(lsm_ncid,idv,glon_lsm)

    ier = nf_inq_dimid(lsm_ncid,'lat',latid)
    if ( ier==nf_noerr ) then
      ier = nf_inq_varid(lsm_ncid,'lat',idv)    
    else
      ier = nf_inq_dimid(lsm_ncid,'latitude',latid)
      call netcdferror(ier)
      ier = nf_inq_varid(lsm_ncid,'latitude',idv)    
    end if
    ier = nf_inq_dimlen(lsm_ncid,latid,iy)
    allocate( glat_lsm(iy) )
    ier = nf_get_var_real(lsm_ncid,idv,glat_lsm)

    if ( size(glon)/=size(glon_lsm) ) then
      write(6,*) "ERROR: Inconsistent longitude when reading SST and land mask"
      call finishbanner
      stop -1
    end if
    if ( size(glat)/=size(glat_lsm) ) then
      write(6,*) "ERROR: Inconsistent latitude when reading SST and land mask"
      call finishbanner
      stop -1
    end if
    if ( any( abs( glon-glon_lsm )>1.e-4 ) ) then
      write(6,*) "ERROR: Inconsistent longitude when reading SST and land mask"
      call finishbanner
      stop -1
    end if
    if ( any( abs( glat-glat_lsm )>1.e-4 ) ) then
      write(6,*) "ERROR: Inconsistent longitude when reading SST and land mask"
      call finishbanner
      stop -1
    end if
    deallocate( glon_lsm, glat_lsm )

    allocate( lsm_gbl(ix*iy) )
    start(1) = 1
    start(2) = 1
    ncount(1) = ix
    ncount(2) = iy
    ier = nf_get_att_real(lsm_ncid,idvar,'add_offset',addoff)
    if ( ier/=nf_noerr ) addoff = 0.
    ier = nf_get_att_real(lsm_ncid,idvar,'scale_factor',sf)
    if ( ier/=nf_noerr ) sf = 1.
    !ier = nf_inq_vartype(lsm_ncid,idvar,itype)
    !call netcdferror(ier)
    !select case(itype)
    !  case ( nf_short )
    !    allocate( ivar(ix*iy) )
    !    ier = nf_get_vara_int(lsm_ncid,idvar,start,ncount,ivar)
    !    call netcdferror(ier)
    !    lsm_gbl(1:ix*iy) = sf*real(ivar(1:ix*iy)) + addoff
    !    deallocate( ivar )
    !  case ( nf_float )
        ier = nf_get_vara_real(lsm_ncid,idvar,start,ncount,datan_tmp)
        call netcdferror(ier)
        lsm_gbl(1:ix*iy) = sf*real(datan_tmp(1:ix*iy)) + addoff     
    !  case ( nf_double )
    !    allocate( dvar(ix*iy) )
    !    ier = nf_get_vara_double(lsm_ncid,idvar,start,ncount,dvar)
    !    call netcdferror(ier)
    !    lsm_gbl(1:ix*iy) = sf*real(dvar(1:ix*iy)) + addoff
    !    deallocate( dvar )      
    !  case default
    !    write(6,*) "Variable is unkown"
    !    call finishbanner
    !    stop -1
    !end select
    ! MJT quick fix 
    where (abs(lsm_gbl(1:ix*iy))>=1.e10)
      lsm_gbl(1:ix*iy)=0.
    end where
    where (lsm_gbl(1:ix*iy)==-1.)
      lsm_gbl(1:ix*iy)=1.
    end where
  else  
    write(6,*) "WARN: Cannot locate land-sea mask"  
    olsm_gbl = .false.
    allocate( lsm_gbl(ix*iy) )
    lsm_gbl = 1.  
  end if
end if
  
write(6,*)"###################### do we have olsm_gbl=",olsm_gbl
write(6,*)"prepare to interp. for sea and land separately"
write(6,*)"putting only land values into datan"
write(6,*)"putting only ocean values into datan_ocn"

nlpnts = 0
nopnts = 0
datan_ocn = datan
if ( olsm_gbl ) then
  do j = 1,iy
    do i = 1,ix
      iq = i+(j-1)*ix
      if ( lsm_gbl(iq) < 0.5 ) then
         datan(iq) = spval          ! land, fill in ocean pts
         nlpnts = nlpnts + 1
      else !!!  ( lsm_gbl(iq) .lt. .5 ) then
         datan_ocn(iq) = spval      ! ocean, fill in land pts
         nopnts = nopnts + 1
      endif ! ( lsm_gbl(iq) .gt. .5 ) then
    enddo ! ix
  enddo ! iy
endif!(olsm_gbl)then

write(6,*)"two global arrays with spval=", spval
write(6,*)"fill in missing values nlpnts,nopnts=",nlpnts,nopnts

write(6,*)"=======> for land array, fill in ocean values"
call fill(datan(1:ix*iy),ix,iy,.1*spval,datan_tmp)

if ( olsm_gbl .and. nopnts>0 ) then
   write(6,*)"=======> for ocean array, fill in land values"
   call fill(datan_ocn,ix,iy,.1*spval,datan_tmp)
endif!(olsm_gbl)then

! interpolation
il = size(sfct,1)
write(6,*)"=========================> now interp. land data"
call sintp16(datan(1:ix*iy),ix,iy,sfct,glon,glat,sdiag,il) ! land
!sfct = min( max( sfct, 100. ), 425. )

if(olsm_gbl .and. nopnts.gt.0)then
   write(6,*)"=========================> now interp. ocean data"
   call sintp16(datan_ocn,ix,iy,sfcto_m,glon,glat,sdiag,il)   ! ocean
   !sfcto_m = min( max( sfcto_m, 100. ), 425. )
endif!(olsm_gbl)then

! recombine
if (olsm_gbl) then
  write(6,*)"now recombine two (land/ocean) fields"    
  do j=1,6*il
    do i=1,il
      iq=i+(j-1)*il
      if ( lsm_m(iq)<0.5 ) sfct(i,j)=sfcto_m(i,j)  ! set to ocean interp pnt
    enddo ! i
  enddo ! j
end if    

deallocate( glon, glat )
deallocate( datan, datan_ocn, datan_tmp )

return
end subroutine readsst

subroutine netcdferror(ier)

use netcdf_m

implicit none

integer, intent(in) :: ier

if ( ier==nf_noerr ) return

write(6,*) "ERROR: ",nf_strerror(ier)
call finishbanner
stop -1

end subroutine netcdferror

subroutine datefix(kdate_r,ktime_r,time,calendar)

implicit none

integer, intent(inout) :: kdate_r,ktime_r
real(kind=8), intent(inout) :: time
integer(kind=8), dimension(12) :: mdays = (/31_8,28_8,31_8,30_8,31_8,30_8,31_8,31_8,30_8,31_8,30_8,31_8/)
integer leap_l
integer(kind=8) mtimer_r
integer(kind=8) iyr,imo,iday,ihr,imins
integer(kind=8) mtimerh,mtimerm
integer(kind=8) mdays_save
integer(kind=8), parameter :: minsday = 1440
character(len=*), intent(in) :: calendar

leap_l = 1
if ( calendar(1:7)=="365_day" ) then
  leap_l = 0
end if

mtimer_r = int(time,8)

iyr=int(kdate_r,8)/10000_8
imo=(int(kdate_r,8)-10000_8*iyr)/100_8
iday=int(kdate_r,8)-10000_8*iyr-100_8*imo
ihr=int(ktime_r,8)/100_8
imins=int(ktime_r,8)-100_8*ihr
write(6,*) 'entering datefix'
write(6,*) 'iyr,imo,iday:       ',iyr,imo,iday
write(6,*) 'ihr,imins,mtimer_r: ',ihr,imins,mtimer_r

mdays(2)=28_8
if ( leap_l==1 ) then
  if ( mod(iyr,4_8)==0   ) mdays(2)=29_8
  if ( mod(iyr,100_8)==0 ) mdays(2)=28_8
  if ( mod(iyr,400_8)==0 ) mdays(2)=29_8
end if
do while ( mtimer_r>minsday*mdays(imo) )
  mtimer_r=mtimer_r-minsday*mdays(imo)
  imo=imo+1_8
  if ( imo>12_8 ) then
    imo=1_8
    iyr=iyr+1_8
    if ( leap_l==1 ) then
      mdays(2)=28_8      
      if ( mod(iyr,4_8)==0   ) mdays(2)=29_8
      if ( mod(iyr,100_8)==0 ) mdays(2)=28_8
      if ( mod(iyr,400_8)==0 ) mdays(2)=29_8
    end if
  end if
end do
write(6,*)'b datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                     iyr,imo,iday,ihr,imins,mtimer_r
  
iday=iday+mtimer_r/minsday
mtimer_r=mod(mtimer_r,minsday)
write(6,*)'c datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                     iyr,imo,iday,ihr,imins,mtimer_r
  
! at this point mtimer_r has been reduced to fraction of a day
mtimerh=mtimer_r/60_8
mtimerm=mtimer_r-mtimerh*60_8  ! minutes left over
ihr=ihr+mtimerh
imins=imins+mtimerm

ihr=ihr+imins/60_8
imins=mod(imins,60_8)
write(6,*)'d datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                     iyr,imo,iday,ihr,imins,mtimer_r
  
iday=iday+ihr/24_8
ihr=mod(ihr,24_8)
write(6,*)'e datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                     iyr,imo,iday,ihr,imins,mtimer_r
  
mdays_save=mdays(imo)
imo=imo+(iday-1_8)/mdays(imo)
iday=mod(iday-1_8,mdays_save)+1_8

iyr=iyr+(imo-1_8)/12_8
imo=mod(imo-1_8,12_8)+1_8

kdate_r=int(iday+100_8*(imo+100_8*iyr),8)
ktime_r=int(ihr*100_8+imins,8)
write(6,*)'end datefix iyr,imo,iday,ihr,imins,mtimer_r: ', &
                       iyr,imo,iday,ihr,imins,mtimer_r
  
write(6,*)'leaving datefix kdate_r,ktime_r: ',kdate_r,ktime_r

time=real(mtimer_r,8)

return
end subroutine datefix

subroutine readpress(ncid,in_type,plev,plev_b,nplev,osig_in,orev)

use netcdf_m

implicit none

integer, intent(in) :: ncid
integer, intent(out) :: nplev
integer ivpres, idpres, ier, ivpres_a, ivpres_b, ivpres_0
integer k, maxplev
real, dimension(:), intent(out) :: plev, plev_b
real, dimension(size(plev)) :: datan
real xplev, plev0
logical, intent(out) :: osig_in, orev
character(len=*), intent(out) :: in_type
character(len=80) presname, typename
character(len=10) presunits

ier = nf_inq_varid(ncid,'pres',ivpres)
ier = nf_inq_dimid(ncid,'pres',idpres)
in_type="p"
if ( ier/=nf_noerr ) then
  ier = nf_inq_varid(ncid,'plev',ivpres)
  ier = nf_inq_dimid(ncid,'plev',idpres)
  in_type="p"
end if
if ( ier/=nf_noerr ) then
  ier = nf_inq_varid(ncid,'lev',ivpres)
  ier = nf_inq_dimid(ncid,'lev',idpres)
  in_type="p"
end if
if ( ier/=nf_noerr ) then
  ier = nf_inq_varid(ncid,'lvl',ivpres)
  ier = nf_inq_dimid(ncid,'lvl',idpres)
  in_type="s"
end if
ier = nf_get_att_text(ncid,ivpres,'type',typename)
if ( ier==nf_noerr ) then
  if ( typename(1:8)=="pressure" ) then
    in_type="p"
  end if
end if
ier = nf_get_att_text(ncid,ivpres,'long_name',presname)
call netcdferror(ier)
if ( presname(1:32) == "hybrid sigma pressure coordinate" ) then
  in_type="h"  
end if
ier = nf_inq_dimlen(ncid,idpres,nplev)
call netcdferror(ier)
ier = nf_get_att_text(ncid,ivpres,'units',presunits)
call netcdferror(ier)
maxplev = size(plev)
if ( maxplev<nplev ) then
  write(6,*) "ERROR: maxplev is less then nplev"
  write(6,*) "This could indicate inconsistant vertical levels"
  write(6,*) "nplev,maxplev =",nplev,maxplev
  call finishbanner
  stop -1
end if
ier = nf_get_var_real(ncid,ivpres,plev(1:nplev))
call netcdferror(ier)
xplev = maxval( plev(1:nplev) )
osig_in = .false.
plev_b = 0.
if ( in_type == "h" ) then
  write(6,*)"^^^^^^^^^hybrid sigma levels^^^^^^^^"
  osig_in = .true.
  ier = nf_inq_varid(ncid,'a',ivpres_a) 
  if ( ier == nf_noerr ) then
    ier = nf_get_var_real(ncid,ivpres_a,plev(1:nplev))
    call netcdferror(ier)
    ier = nf_inq_varid(ncid,'b',ivpres_b)  
    call netcdferror(ier)
    ier = nf_get_var_real(ncid,ivpres_b,plev_b(1:nplev))
    call netcdferror(ier)
    ier = nf_inq_varid(ncid,'p0',ivpres_0)
    call netcdferror(ier)
    ier = nf_get_var_real(ncid,ivpres_0,plev0)
    call netcdferror(ier)
    plev0 = plev0/100. ! convert to hPa
    do k = 1,nplev
      plev_b(k) = plev_b(k)*plev0
      plev(k) = plev(k)*1000.
    end do
  else
    ier = nf_inq_varid(ncid,'ap',ivpres_a)  
    call netcdferror(ier)
    ier = nf_get_var_real(ncid,ivpres_a,plev_b(1:nplev)) ! note ap is for plev_b
    call netcdferror(ier)
    ier = nf_inq_varid(ncid,'b',ivpres_b)
    call netcdferror(ier)
    ier = nf_get_var_real(ncid,ivpres_b,plev(1:nplev)) ! note b is for plev
    call netcdferror(ier)
    do k = 1,nplev
      plev_b(k) = plev_b(k)/100. ! convert to hPa
      plev(k) = plev(k)*1000.
    end do    
  end if    
  presunits="hPa"
else if ( .01<xplev .and. xplev<800.  ) then
  write(6,*)"^^^^^^^^^actualy sigma levels^^^^^^^ fix plevs"
  osig_in = .true.
  do k = 1,nplev
    plev(k) = plev(k)*1000.
  end do
  presunits="hPa"
else if ( xplev <= 0.01 ) then
  write(6,*)"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ fix plevs"
  write(6,*) 'xplev < 0.01 in cdfvidar'
  call finishbanner
  stop -1
end if
if ( osig_in ) then
  orev = plev(nplev)+plev_b(nplev) > plev(1)+plev_b(1)
else
  orev = plev(nplev)>plev(1)
end if
if ( .not.orev ) then
  do k = 1,nplev
    datan(k) = plev(k)
  end do
  do k = 1,nplev
    plev(k) = datan(nplev+1-k)
  end do
  do k = 1,nplev
    datan(k) = plev_b(k)
  end do
  do k = 1,nplev
    plev_b(k) = datan(nplev+1-k)
  end do
end if
if ( presunits(1:2)=="Pa" ) then
  plev = plev/100.
  plev_b = plev_b/100.
  presunits="hPa"
end if

if ( presunits(1:3)/="hPa" ) then
  write(6,*) "ERROR: Could not convert vertical levels to hPa"
  write(6,*) "Vertical level units was read as ",trim(presunits)
  call finishbanner
  stop -1
end if
  
return
end subroutine readpress

end module outcdf_m
