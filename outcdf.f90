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
    
      subroutine outcdf(ihr,idy,imon,iyr,iout,nt,time,mtimer,sig,cdffile,ddss,il,kl)

      use netcdf_m
      
      implicit none
      
      integer il,jl,kl,ifull
      integer ihr,idy,imon,iyr,iout
      integer nt
      real ddss
      real du,tanl,rnml,stl1,stl2

      !include 'newmpar.h'
!     include 'darcdf.h'   ! idnc,ncid,idifil  - stuff for netcdf
!     include 'dates.h' ! ktime,kdate,timer,timeg,xg,yg,mtimer
!     include 'filnames.h'  ! list of files, read in once only
!     include 'kuocom.h'
!     include 'liqwpar.h'  ! ifullw
      common/mapproj/du,tanl,rnml,stl1,stl2
!     include 'parm.h'
!     include 'parmdyn.h'  
!     include 'parmhor.h'  ! mhint, m_bs, nt_adv, ndept
!     include 'parmvert.h'

      character rundate*10

      integer, save :: idnc0, idnc1, idncm1
      integer nextout, itype, nihead, nrhead
      integer kdate, ktime, iarch
      integer idnc, ier, imode, ixp, iyp
      integer idlev, idnt, ktau, icy, icm, icd
      integer ich, icmi, ics, mtimer, idv
      parameter(nextout=0,itype=1)
      parameter(nihead=54)
      integer nahead(nihead)
      integer, dimension(1) :: dimids

      parameter(nrhead=14)
      real ahead(nrhead)

      real sig(kl)
      real time,dt,ds

      character cdffile*80

      common/cdfind/ixp,iyp,idlev,idnt

      integer dim(4)
      integer xdim,ydim,zdim,tdim
      integer oldmode
      character timorg*20
      character grdtim*33
      character*3 month(12)
      data month/'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'/

      data idnc1/0/, idnc0/0/, idncm1/0/
      data rundate/"ncepavnanl"/

      jl=6*il
      ifull=il*jl

      write(6,*)"outcdf ihr,idy,imon,iyr,iout=",ihr,idy,imon,iyr,iout
      write(6,*)"time=",time
!     write(6,*)"sig=",sig
!     write(6,*)"ddss=",ddss

      dt=0
      kdate=iyr*10000+imon*100+idy
      ktime=ihr*100
      write(6,*)"kdate,ktime=",kdate,ktime

! itype=1 outfile
      iarch=nt
      idnc=idnc1

      write(6,'("outcdf itype,idnc,iarch,cdffile=",3i5," ",a80)') itype,idnc,iarch,cdffile

!#######################################################################
! netcdf output
!#######################################################################

      if ( iarch.lt.1 ) stop "wrong iarch in outcdf"
      if ( iarch.eq.1 ) then
        write(6,*)'nccre of ',cdffile
#ifdef usenc3
      !idnc = nccre(cdffile, ncclob, ier)
      !ier = nf_create(cdffile,nf_clobber,idnc)
      ier = nf_create(cdffile,NF_64BIT_OFFSET,idnc)
#else
	!ier=nf_create(cdffile,NF_NOCLOBBER,idnc)
	!ier=nf_create(cdffile,NF_64BIT_OFFSET,idnc)
	!ier=nf_create(cdffile,NF_NETCDF4.or.NF90_CLASSIC_MODEL,idnc)
	ier=nf_create(cdffile,NF_NETCDF4,idnc)
#endif
        write(6,*)'idnc,ier=',idnc,ier
! Turn off the data filling
        ier = nf_set_fill(idnc,nf_nofill,oldmode)
        !write(6,*)'imode=',imode
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
        ier = nf_def_var(idnc,'time',nf_int,1,dimids,idnt)
        write(6,*)'idnt=',idnt
        ier = nf_put_att_text(idnc,idnt,'point_spacing',4,'even')

        !write(6,*)'kdate,ktime,ktau=',kdate,ktime,ktau
        write(6,*)'kdate,ktime=',kdate,ktime

        icy=kdate/10000
        icm=max(1,min(12,(kdate-icy*10000)/100))
        icd=max(1,min(31,(kdate-icy*10000-icm*100)))
        ! disabled by MJT
        !if(icy.lt.100)icy=icy+1900
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

        dim(1) = xdim
        dim(2) = ydim
        dim(3) = zdim
        dim(4) = tdim

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

        dimids = 0
        ier = nf_def_var(idnc,'ds',nf_float,0,dimids,idv)
        if(ier.ne.0)write(6,*)"ncvdef ds idnc,ier=",idnc,ier
        ier = nf_def_var(idnc,'du',nf_float,0,dimids,idv)
        if(ier.ne.0)write(6,*)"ncvdef du idnc,ier=",idnc,ier
        ier = nf_def_var(idnc,'rnml',nf_float,0,dimids,idv)
        if(ier.ne.0)write(6,*)"ncvdef rnml idnc,ier=",idnc,ier
        ier = nf_def_var(idnc,'tanl',nf_float,0,dimids,idv)
        if(ier.ne.0)write(6,*)"ncvdef tanl idnc,ier=",idnc,ier
        ier = nf_def_var(idnc,'stl1',nf_float,0,dimids,idv)        
        if(ier.ne.0)write(6,*)"ncvdef stl1 idnc,ier=",idnc,ier
        ier = nf_def_var(idnc,'stl2',nf_float,0,dimids,idv)
        if(ier.ne.0)write(6,*)"ncvdef stl2 idnc,ier=",idnc,ier
        ier = nf_def_var(idnc,'dt',nf_float,0,dimids,idv)        
        if(ier.ne.0)write(6,*)"ncvdef dt idnc,ier=",idnc,ier
      endif ! ( iarch=1 ) then

      write(6,*)'call openhist for itype= ',itype
      call openhist(idnc,iarch,itype,dim,sig,kdate,ktime,time,mtimer,il,kl)

      ier = nf_sync(idnc)
      if(ier.ne.0)write(6,*)"ncsnc idnc,ier=",idnc,ier

      if ( itype.eq.1 ) then
!       itype=1 outfile
        idnc1=idnc
      elseif ( itype.eq.0 ) then
!       itype=0 climcdf
        idnc0=idnc
      elseif ( itype.eq.-1 ) then
!       itype=-1 restfile
        idncm1=idnc
      endif ! ( itype.eq.1 ) then

      return ! outcdf
      end
!=======================================================================
      subroutine openhist(idnc,iarch,itype,dim,sig,kdate,ktime,time,mtimer,il,kl)

      use cll_m
      use netcdf_m
      use sigdata_m
      use xyzinfo_m, only : em,f

!     this routine creates attributes and writes output

      integer il,jl,kl,ifull

      !include 'newmpar.h'
!     include 'aalat.h'
!     include 'arrays.h'
!     include 'darcdf.h'   ! idnc,ncid,idifil  - stuff for netcdf
!     include 'dates.h' ! ktime,kdate,timer,timeg,xg,yg,mtimer
!     include 'extraout.h'
!     include 'filnames.h' ! list of files, read in once only
!     include 'kuocom.h'
!     include 'liqwpar.h'  ! ifullw
      !include 'map.h'
!     include 'mapproj.h'
!     include 'morepbl.h'
!     include 'nsibd.h' ! rsmin,ivegt,sigmf,tgg,tgf,ssdn,res,rmc,isoilm,ico2em
      !include 'parm.h'
!     include 'parmdyn.h'
!     include 'parmvert.h'
!     include 'pbl.h'
!     include 'prec.h'
!     include 'scamdim.h'
!     include 'screen.h'
!     include 'sigs.h'
      !common/cll/clon(ifull),clat(ifull)

      !include 'sigdata.h'
!n    common/sigdata/pmsl(ifull),sfct(ifull),zs(ifull),ps(ifull)
!n   &             ,us(ifull,kl)    ,vs(ifull,kl)    ,ts(ifull,kl)
!n   &             ,rs(ifull,kl)    ,hs(ifull,kl)    ,psg_m(ifull)
!n   &             ,zsi_m(ifull)
!o    common/sigdata/pmsl(ifull),sfct(ifull),zs(ifull),ps(ifull)
!o   &             ,u(ifull,kl)    ,v(ifull,kl)    ,t(ifull,kl)
!o   &             ,qg(ifull,kl)   ,hs(ifull,kl)

      character lname*50,expdesc*50
      integer dim(4)
      integer idim2(3)
      integer, dimension(1) :: ivals
      integer, dimension(1) :: start, ncount
      integer nrun
      real xpnt(il),ypnt(6*il)
      real sig(kl)

      common/cdfind/ixp,iyp,idlev,idnt
      real tst(6*il*il),tsb(6*il*il)
!       *** qscrn_ave not presently written     
      real aa(6*il*il),bb(6*il*il),cc(6*il*il)
      real cfrac(6*il*il,kl)
      real, dimension(1) :: rvals

      jl=6*il
      ifull=il*jl

!     character*3 mon(12)
!     data mon/'JAN','FEB','MAR','APR','MAY','JUN'
!    &        ,'JUL','AUG','SEP','OCT','NOV','DEC'/

      write(6,*)'openhist iarch,idnc,itype=',iarch,idnc,itype

!     if(itype.ne.-1)then  ! don't scale up for restart file as done already
!       insert stuff here if re-scaling clouds etc
!     endif  ! (itype.ne.-1)

!     if this is the first archive, set up some global attributes
      if(iarch.eq.1) then
        write(6,*)'dim=',dim
        idim2(1)=dim(1)
        idim2(2)=dim(2)
        idim2(3)=dim(4)
        write(6,*)'idim2=',idim2

!       Create global attributes
!       Model run number
        nrun = 1
        ivals = nrun
        ier = nf_put_att_int(idnc,nf_global,'nrun',nf_int,1,ivals)
        write(6,*)"nrun=",nrun," ier=",ier

!       Experiment description
        expdesc = 'CCAM model run'
        ier = nf_put_att_text(idnc,nf_global,'expdesc',50,expdesc)
        write(6,*)"expdesc=",expdesc," ier=",ier

!       Sigma levels
!       write(6,*)'sig=',sig
        ier = nf_put_att_real(idnc,nf_global,'sigma',nf_float,kl,sig)

        lname = 'year-month-day at start of run'
        ier = nf_def_var(idnc,'kdate',nf_int,1,dim(4:4),idkdate)
        ier = nf_put_att_text(idnc,idkdate,'long_name',lngstr(lname),lname)

        lname = 'hour-minute at start of run'
        ier = nf_def_var(idnc,'ktime',nf_int,1,dim(4:4),idktime)
        ier = nf_put_att_text(idnc,idktime,'long_name',lngstr(lname),lname)

        lname = 'timer (hrs)'
        ier = nf_def_var(idnc,'timer',nf_float,1,dim(4:4),idnter)
        ier = nf_put_att_text(idnc,idnter,'long_name',lngstr(lname),lname)

        lname = 'mtimer (mins)'
        ier = nf_def_var(idnc,'mtimer',nf_int,1,dim(4:4),idmtimer)
        ier = nf_put_att_text(idnc,idmtimer,'long_name',lngstr(lname),lname)

        lname = 'timeg (UTC)'
        ier = nf_def_var(idnc,'timeg',nf_float,1,dim(4:4),idnteg)
        ier = nf_put_att_text(idnc,idnteg,'long_name',lngstr(lname),lname)

        lname = 'number of time steps from start'
        ier = nf_def_var(idnc,'ktau',nf_int,1,dim(4:4),idktau)
        ier = nf_put_att_text(idnc,idktau,'long_name',lngstr(lname),lname)

        ier = nf_def_var(idnc,'sigma',nf_float,1,dim(3:3),idv)
        ier = nf_put_att_text(idnc,idv,'positive',4,'down')

        write(6,*)'define attributes of variables'

        lname ='Scaled Log Surface pressure'
        call attrib(idnc,idim2,3,'psf',lname,'none',-1.3,0.2)

        lname ='Mean sea level pressure'
        call attrib(idnc,idim2,3,'pmsl',lname,'hPa',800.,1200.)
        lname = 'Surface geopotential'
        call attrib(idnc,idim2,2,'zht',lname,'m2/s2',-2.e3,128.e3) ! MJT lsmask
!       call attrib(idnc,idim2,2,'zht',lname,'m2/s2',-100.,90.e3)  ! MJT lsmask
!       call attrib(idnc,idim2,2,'zht',lname,'m2/s2',-1.e6,90.e3) ! ocean too
        lname = 'Soil type'                                        ! MJT lsmask
        call attrib(idnc,idim2,2,'soilt',lname,'none',0.,65.e3)    ! MJT lsmask

!       For time invariant surface fields
        lname = 'Map factor'
        call attrib(idnc,idim2,2,'map',lname,'none',0.,20.)
        lname = 'Coriolis factor'
        call attrib(idnc,idim2,2,'cor',lname,'1/sec',-1.5e-4,1.5e-4)
!       lname = 'Initial wetness fraction layer 3'
!       call attrib(idnc,idim2,2,'wetfrac',lname,'none',-2.,2.)
        lname = 'clon'
        call attrib(idnc,idim2,2,'clon',lname,'none',-360.,360.)
        lname = 'clat'
        call attrib(idnc,idim2,2,'clat',lname,'none',-90.,90.)

!       For time varying surface fields
        lname = 'Surface temperature'
        call attrib(idnc,idim2,3,'tsu',lname,'K',0.,350.)
!       lname = 'Runoff'
!       call attrib(idnc,idim2,3,'runoff',lname,'mm/day',0.,1000.)
!       lname = 'Sea ice depth (Instantaneous)'
!       call attrib(idnc,idim2,3,'siced',lname,'cm',0.,500.)
!       lname = 'snow depth (liquid water)'
!       call attrib(idnc,idim2,3,'snd',lname,'cm',0.,1000.)

!       lname = 'Soil moisture as frac FC levels 1-2'
!       call attrib(idnc,idim2,3,'wbfshal',lname,'frac',0.,4.)
!       lname = 'Soil moisture as frac FC levels 3-4'
!       call attrib(idnc,idim2,3,'wbfroot',lname,'frac',0.,4.)
!       lname = 'Soil moisture as frac FC levels 1-6'
!       call attrib(idnc,idim2,3,'wbftot',lname,'frac',0.,4.)

!       lname = 'Soil temperature lev 1'
!       call attrib(idnc,idim2,3,'tgg1',lname,'K',100.,400.)
!       lname = 'Soil temperature lev 2'
!       call attrib(idnc,idim2,3,'tgg2',lname,'K',100.,400.)
!       lname = 'Soil temperature lev 3'
!       call attrib(idnc,idim2,3,'tgg3',lname,'K',100.,400.)
!       lname = 'Soil temperature lev 4'
!       call attrib(idnc,idim2,3,'tgg4',lname,'K',100.,400.)
!       lname = 'Soil temperature lev 5'
!       call attrib(idnc,idim2,3,'tgg5',lname,'K',100.,400.)
!       lname = 'Soil temperature lev 6'
!       call attrib(idnc,idim2,3,'tgg6',lname,'K',100.,400.)
        lname = 'Soil temperature top'
        call attrib(idnc,idim2,3,'tb3',lname,'K',100.,400.)
        lname = 'Soil temperature bottom'
        call attrib(idnc,idim2,3,'tb2',lname,'K',100.,400.)
        lname = 'Soil moisture top'
        call attrib(idnc,idim2,3,'wfg',lname,'none',0.,.4)
        lname = 'Soil moisture bottom'
        call attrib(idnc,idim2,3,'wfb',lname,'none',0.,.4)

        if (any(fracice.ge.0.)) then
          lname = 'Sea ice fraction'
          call attrib(idnc,idim2,3,'fracice',lname,'none',0.,6.5)
        end if

        write(6,*)'3d variables'
        call attrib(idnc,dim,4,'temp','Air temperature','K',100.,350.)
        call attrib(idnc,dim,4,'u','x-component wind','m/s',-150.,150.)
        call attrib(idnc,dim,4,'v','y-component wind','m/s',-150.,150.)
        lname= 'Water mixing ratio'
        call attrib(idnc,dim,4,'mixr',lname,'kg/kg',0.,.05)
        !if(ifullw.eq.ifull)then
          call attrib(idnc,dim,4,'qfg','Frozen water','kg/kg',0.,.02)
          call attrib(idnc,dim,4,'qlg','Liquid water','kg/kg',0.,.02)
          call attrib(idnc,dim,4,'cfrac','Cloud fraction','none',0.,1.)
        !endif

        write(6,*)'finished defining attributes'
!       Leave define mode
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

        ier = nf_inq_varid(idnc,'sigma',idv)
        start = 1
        ncount = kl
        ier = nf_put_vara_real(idnc,idv,start,ncount,sig)

        !ier = nf_inq_varid(idnc,'lev',idv)
        !start = 1
        !ncount = kl
        !ier = nf_put_vara_real(idnc,idv,start,ncount,sig)
!       do k = kl,1,-1
!         write(6,*)"k=",k," sig=",sig(k)
!       enddo

!        ier = nf_inq_varid(idnc,'ds',idv)
!        rvals = ds
!        ier = nf_put_var_real(idnc,idv,rvals)
!        ier = nf_inq_varid(idnc,'tanl',idv)
!        rvals = tanl
!        ier = nf_put_var_real(idnc,idv,rvals)
!        ier = nf_inq_varid(idnc,'rnml',idv)
!        rvals = rnml
!        ier = nf_put_var_real(idnc,idv,rvals)
!        ier = nf_inq_varid(idnc,'du',idv)
!        rvals = du
!        ier = nf_put_var_real(idnc,idv,rvals)
!        ier = nf_inq_varid(idnc,'stl1',idv)
!        rvals = stl1
!        ier = nf_put_var_real(idnc,idv,rvals)
!        ier = nf_inq_varid(idnc,'stl2',idv)
!        rvals = stl2
!        ier = nf_put_var_real(idnc,idv,rvals)
!        ier = nf_inq_varid(idnc,'dt',idv)
!        rvals = dt
!        ier = nf_put_var_real(idnc,idv,rvals)
      endif ! iarch.eq.1
     
!------------------------------------------------------------------      
      ktau=0
      !timer=0
      timer=time/60. ! MJT quick fix
      !timeg=0
      timeg=mod(time/60.,24.) ! MJT quick fix
      write(6,*)'outcdf processing kdate,ktime,ktau,time,mtimer: ',kdate,ktime,ktau,time,mtimer

!     set time to number of minutes since start
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
!       write time-invariant fields      
        call histwrt3(real(em),'map',idnc,iarch,il)
        call histwrt3(real(f),'cor',idnc,iarch,il)
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

      call histwrt3(sfct,'tsu',idnc,iarch,il)

      !call histwrt3(sfct,'tgg1',idnc,iarch,il)

      !call histwrt3(ts(1,2),'tgg2',idnc,iarch,il)
      !call histwrt3(ts(1,2),'tgg3',idnc,iarch,il)
      !call histwrt3(ts(1,2),'tgg4',idnc,iarch,il)
      !call histwrt3(ts(1,2),'tgg5',idnc,iarch,il)
      !call histwrt3(ts(1,2),'tgg6',idnc,iarch,il)

      call histwrt3(ts(:,2),'tb3',idnc,iarch,il) ! top
      call histwrt3(ts(:,2),'tb2',idnc,iarch,il) ! bottom

      aa(1:ifull)=0.14
      !call histwrt3(aa,'wbfshal',idnc,iarch,il)
      !call histwrt3(aa,'wbfroot',idnc,iarch,il)
      !call histwrt3(aa,'wbftot',idnc,iarch,il)
      call histwrt3(aa,'wfg',idnc,iarch,il)
      call histwrt3(aa,'wfb',idnc,iarch,il)

      if ( any( fracice >= 0. ) ) then
        call histwrt3(fracice,'fracice',idnc,iarch,il)
      end if

!     call histwrt3(sicedep,'siced',idnc,iarch)
!     call histwrt3(snowd,'snd',idnc,iarch)
      
      write(6,*)'netcdf save of 3d variables'
      call histwrt4(ts,'temp',idnc,iarch,il,kl)
      call histwrt4(us,'u',idnc,iarch,il,kl)
      call histwrt4(vs,'v',idnc,iarch,il,kl)
      call histwrt4(rs,'mixr',idnc,iarch,il,kl)
      !write(6,*)"ifullw,ifull=",ifullw,ifull

      !if(ifullw.eq.ifull)then
      !  call histwrt4(qfg,'qfg',idnc,iarch,il,kl)
      !  call histwrt4(qlg,'qlg',idnc,iarch,il,kl)
      !  call histwrt4(cfrax,'cfrac',idnc,iarch,il,kl)
      !endif

      return ! subroutine openhist(idnc,iarch,itype,dim,sig
      end
!=======================================================================
      subroutine attrib(cdfid,dim,ndim,name,lname,units,xmin,xmax)

      use netcdf_m

      integer*2 minv, maxv, missval   ! was integer*2
      parameter(minv = -32500, maxv = 32500, missval = -32501)
      integer ndim
      integer cdfid, idv, dim(ndim)
      character name*(*), lname*(*), units*(*)
      real xmin, xmax
      real, dimension(1) :: rvals
      integer lngstr
      integer(kind=2), dimension(1) :: i2vals

      ier = nf_def_var(cdfid, name, nf_int2, ndim, dim, idv)
      if ( ier.ne.0 ) then
        write(6,*)ier,' Error in variable declaration ', name
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
!     scalef=(xmax-xmin)/float(maxv - minv)
      scalef=(xmax-xmin)/(real(maxv)-real(minv)) ! jlm fix for precision problems
      addoff=xmin-scalef*minv
      rvals = addoff
      ier = nf_put_att_real(cdfid,idv,'add_offset',nf_float,1,rvals)
      rvals = scalef
      ier = nf_put_att_real(cdfid,idv,'scale_factor',nf_float,1,rvals)
      ier = nf_put_att_text(cdfid,idv,'FORTRAN_format',5,'G11.4')
      return
      end
!=======================================================================
      function lngstr( string )
      character*(*) string
      ilen = len(string)
!     print*,'string=',string
!     print*,'ilen=',ilen
      do 100 lngstr=ilen,1,-1
        if ( string(lngstr:lngstr) .ne. ' ' ) go to 99
  100 continue
      lngstr = 0
   99 continue
      return
      end
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
      integer(kind=2) ipack(il,6*il) ! was integer*2 
      character(len=*), intent(in) :: sname
!     character*8 sname
      integer(kind=2) minv, maxv, missval ! was integer*2 
      parameter(minv = -32500, maxv = 32500, missval = -32501)
      real addoff, scale_f
      real varn, varx, xmin, xmax, pvar
      integer, intent(in) :: idnc
      integer ier

      real var(il,6*il)

      jl=6*il
      ifull=il*jl

      write(6,*)"histwrt3 sname=",sname," iarch=",iarch," idnc=",idnc

      start(1) = 1
      start(2) = 1
      start(3) = iarch
      count(1) = il
      count(2) = jl
      count(3) = 1

! find variable index
      ier = nf_inq_varid(idnc,sname,mid)
      ier = nf_get_att_real(idnc,mid,'add_offset',addoff)
      ier = nf_get_att_real(idnc,mid,'scale_factor',scale_f)

      xmin=addoff+scale_f*minv
!     xmax=xmin+scale_f*float(maxv-minv)
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
      ier = nf_put_vara_int2(idnc, mid, start(1:ndims), count(1:ndims), ipack)
      if(ier.ne.0)stop "in histwrt3 ier not zero"

      write(6,'("histwrt3:",a7," nt=",i4," n=",f12.4," ij=",2i4," x=",f12.4," ij=",2i4)') sname,iarch,varn,imn,jmn,varx,imx,jmx

      return
      end ! histwrt3
!=======================================================================
      subroutine histwrt4(var,sname,idnc,iarch,il,kl)
! Write 3d+t fields from the savegrid array.

      use netcdf_m
      
      integer il,jl,kl,ifull

      !include 'newmpar.h'
      !include 'parm.h'

      integer mid, start(4), count(4)
      integer*2 ipack(il,6*il,kl) ! was integer*2 
      character(len=*), intent(in) :: sname
!     character*8 sname
      integer*2 minv, maxv, missval ! was integer*2 
      parameter(minv = -32500, maxv = 32500, missval = -32501)
      real addoff, scale_f

      real var(il,6*il,kl)

      jl=6*il
      ifull=il*jl

      write(6,*)"histwrt4 sname=",sname," iarch=",iarch," idnc=",idnc

      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = iarch
      count(1) = il
      count(2) = jl
      count(3) = kl
      count(4) = 1

! find variable index
      ier = nf_inq_varid(idnc,sname,mid)
      ier = nf_get_att_real(idnc,mid,'add_offset',addoff)
      ier = nf_get_att_real(idnc,mid,'scale_factor',scale_f)

      xmin=addoff+scale_f*minv
!     xmax=xmin+scale_f*float(maxv-minv)
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

      ier = nf_put_vara_int2(idnc, mid, start, count, ipack)

      write(6,'("histwrt4:",a7," nt=",i4," n=",f12.4," ijk=",3i4," x=",f12.4," ijk=",3i4)') &
          sname,iarch,varn,imn,jmn,kmn,varx,imx,jmx,kmx

      return
      end ! histwrt4
     
      subroutine mtimerget(mtimer,kdate1,ktime1,kdate2,ktime2) ! jlm
!     returns mtimer in minutes, corr. to (kdate2,ktime2) -  (kdate1,ktime1)    
      dimension ndoy(12)   ! days from beginning of year (1st Jan is 0)
      data ndoy/ 0,31,59,90,120,151,181,212,243,273,304,334/
      common/leap_yr/leap  ! 1 to allow leap years
 
      if(leap.ne.0)stop 'leap years not catered for in mtimerget'
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
      end
