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
    
      program cdfvidar

      use ccinterp
      use cll_m
      use comsig_m      
      use latlong_m
      use netcdf_m
      use outcdf_m
      use sigdata_m
      
      implicit none
      
      include 'version.h'
      
      integer iout,in
      integer kvin,kuin,kbin
      integer ism,nrh
      integer krin,ktin,khin
      integer mtimer,id,jd
      integer mxcyc,nvsig,ntimes
      integer insm,l
      integer ilx,jlx,i,j,iq
      integer nplev,ilonx,ilatx
      integer ier,ijd,ilt,jlk
      integer irecd,ngatts,nvars
      integer ndims,idv,lonid,latid
      integer k,ivpres
      integer ivpres_a, ivpres_b
      integer ivpres_0
      integer narch,idtim,idpres
      integer ilonn,iyr,imn,idy
      integer ivtim,ihr,imi
      integer icmonth_to_imn
      integer nt,iarch,khout,igout
      integer idvar,kn,kx
      integer nopnts,nlpnts,ijgd
      integer jk,imidpan2
      integer ilatn,if2,inlsavn
      integer inzsavn,ints,inzs
      integer io_out,igd,jgd
      
      real ptop,grdx,grdy
      real g,pi,rlonx,rlatx
      real rlonn,rlatn
      real rlong0,rlat0,schmidt
      real elon,elat,xplev,tohpa
      real xa,an
      real spval,sinlong,coslat
      real coslong,polenz,poleny
      real polenx,sinlat,sx,cx
      real zonx,zony,sn,cn
      real sinth,conth,den,zonz
      real uzon,costh,ull,vll
      real vmer,ds
      real rlon
      real llrng,minlon,maxlon,minlat,maxlat
      real plev0
      real(kind=8) time, newtime

      character(len=1024) inf
      character(len=1024) t_file, rh_file, u_file, v_file, z_file
      character(len=1024) lsm_file, zs_file, ps_file, psl_file, ts_file
      character(len=1024) sic_file, snod_file, soiltemp_file, soilmois_file
      
      integer t_ncid, rh_ncid, u_ncid, v_ncid, z_ncid
      integer lsm_ncid, zs_ncid, ps_ncid, psl_ncid, ts_ncid
      integer sic_ncid, snod_ncid, soiltemp_ncid, soilmois_ncid
      integer kdate, ktime

      include 'nplevs.h' ! maxplev
      real, dimension(maxplev) :: plev,cplev,bplev
      real, dimension(maxplev) :: plev_b

      include 'lmax.h'

      parameter ( pi=3.1415926536 )
      parameter ( g=9.80616 )
      
      real, dimension(lmax) :: dsg,sgml

! dimension ids
      integer  dimil,dimjl,dimkl,dimtim
! variable ids
      integer  idil,idjl,idkl,idnt,ix,iy
      integer  il,jl,kl,ifull
      integer iernc,liernc,lncid,varid,dimid
      integer idsoillvl, ivsoillvl, nsoillvl, ksearch, ktest
      integer, dimension(2) :: ccdim
      integer, dimension(2) :: start, ncount

      include 'vidar.h'
      logical sdiag

      real dst, xfrac
      real, dimension(2) :: lonlat
      real, dimension(:,:,:), allocatable :: rlld,xyz,axyz,bxyz      
      real, dimension(:,:), allocatable :: grid
      real, dimension(:), allocatable :: ax,ay,az,bx,by,bz,x,y,z
      real, dimension(:), allocatable :: datan

      real plevin(maxplev)
      real, dimension(:,:,:), allocatable :: hgt
      real, dimension(:,:,:), allocatable :: temp
      real, dimension(:,:,:), allocatable :: u
      real, dimension(:,:,:), allocatable :: v
      real, dimension(:,:,:), allocatable :: rh
      real, dimension(:,:,:), allocatable :: soildim
      real, dimension(:,:), allocatable :: twodim
      real, dimension(:), allocatable :: lsmg_m
      real, dimension(:), allocatable :: zsg
      real, dimension(:), allocatable :: validlevcc
      real, dimension(:), allocatable :: soildepth_in
      real, dimension(6), parameter :: soildepth_ccam = (/ 0.011, 0.051, 0.157, 0.4385, 1.1855, 3.164 /)
      real, dimension(6), parameter :: soilthick_ccam = (/ 0.022, 0.058, 0.154, 0.409,  1.085,  2.872 /)
      real, parameter :: rhowater = 1000.
      real fill_float

      character(len=1) in_type
      character(len=2) moist_var

      logical ofirst, ogbl, orev, olsm_gbl
      character(len=60) timorg
      character(len=60) cu
      character(len=3) cmonth
      character(len=80) zsavn,lsavn
      character(len=10) header,moistvar
      character(len=10) soilunits, geopotunits, presunits
      character(len=20) calendar
      character(len=80) presname
 
      namelist/gnml/inf,vfil,ds,inzs,zsfil                          &
                   ,ints,tsfil, ogbl,zsavn,inzsavn,lsavn,inlsavn    &
                   ,plevin,orev,io_out,igd,jgd,id,jd,mtimer,ntimes  &
                   ,spline,mxcyc,nvsig,nrh                          &
                   ,oesig,sgml,dsg,ptop,debug,notop,opre,have_gp    &
                   ,in,calout                                       &
                   ,iout,oform,sdiag                                &
                   ,insm,smfil                                      &
                   ,splineu,splinev,splinet,zerowinds               &
                   ,grdx,grdy,slon,slat                             &
                   ,moistvar,kl,llrng                               &
                   ,t_file,rh_file,u_file,v_file,z_file             &
                   ,lsm_file,zs_file,ps_file,psl_file,ts_file       &
                   ,sic_file,snod_file,soiltemp_file,soilmois_file

      data khin/0/,kuin/0/,kvin/0/,ktin/0/,krin/0/
      data igd/1/,jgd/1/,id/1/,jd/1/,mtimer/0/
      data ofirst/.true./,io_out/3/
      data ogbl/.true./,orev/.false./
      data inzs/10/, ints/0/, insm/0/, in/50/, iout/70/
      data ptop/0./, ntimes/1/, nvsig/4/
      data mxcyc/20/, oform/.true./
      data spline/.true./, notop/.false./, opre/.false./
      data debug/.false./,oesig/.true./, calout/.true./, nrh/0/
      data inzsavn/11/ , zsavn/'zsavn.ff'/
      data inlsavn/11/ , lsavn/'lsavn.ff'/
      data moistvar/''/
      data zsfil/'/tmp/csjjk/topog5'/
      data tsfil/'/tmp/csjjk/sfct'/
      data smfil/'/tmp/csjjk/smfil'/
      data sgml/lmax*0./, dsg/lmax*0./
      data plevin/maxplev*0./
      data splineu/.true./, splinev/.true./, splinet/.true./
      data sdiag/.false./
      data have_gp/.true./
      data zerowinds/.true./
      data grdx/1./
      data grdy/-1./
      data llrng/0./

      ! Start banner
      write(6,*) "=============================================================================="
      write(6,*) "CCAM: Starting cdfvidar"
      write(6,*) "============================================================================="
      write(6,*) version

#ifndef stacklimit
      ! For linux only - removes stacklimit on all processors
      call setstacklimit(-1)
#endif 
   
      slon=0.
      slat=90.
      
      inf = " "
      t_file = " "
      rh_file = " "
      u_file = " "
      v_file = " "
      z_file = " "
      lsm_file = " "
      zs_file = " "
      ps_file = " "
      psl_file = " "
      ts_file = " "
      sic_file = " "
      snod_file = " "
      soiltemp_file = " "
      soilmois_file = " "
      
      ds = 0.

!####################### read namelists ############################
      write(6,*)'read namelist'
      read (5, gnml)
      write(6,nml=gnml)
!####################### read namelist ############################

      if ( t_file == " " ) t_file=inf
      if ( rh_file == " " ) rh_file=inf
      if ( u_file == " " ) u_file=inf
      if ( v_file == " " ) v_file=inf
      if ( z_file == " " ) z_file=inf
      if ( lsm_file == " " ) lsm_file=inf
      if ( zs_file == " " ) zs_file=inf
      if ( ps_file == " " ) ps_file=inf
      if ( psl_file == " " ) psl_file=inf
      if ( ts_file == " " ) ts_file=inf
      if ( sic_file == " " ) sic_file=inf
      if ( snod_file == " " ) snod_file=inf
      if ( soiltemp_file == " " ) soiltemp_file=inf
      if ( soilmois_file == " " ) soilmois_file=inf
      
      ! open files
      write(6,*)'inf='
      write(6,*) trim(t_file)
      ier = nf_open(t_file,nf_nowrite,t_ncid)
      write(6,*)'ncid=',t_ncid
      if(ier.ne.0) then
        write(6,*)' cannot open netCDF file; error code ',ier
        call finishbanner
        stop -1
      end if
      ier = nf_open(rh_file,nf_nowrite,rh_ncid)
      ier = nf_open(u_file,nf_nowrite,u_ncid)
      ier = nf_open(v_file,nf_nowrite,v_ncid)
      ier = nf_open(z_file,nf_nowrite,z_ncid)
      ier = nf_open(lsm_file,nf_nowrite,lsm_ncid)
      ier = nf_open(zs_file,nf_nowrite,zs_ncid)
      ier = nf_open(ps_file,nf_nowrite,ps_ncid)
      ier = nf_open(psl_file,nf_nowrite,psl_ncid)
      ier = nf_open(ts_file,nf_nowrite,ts_ncid)
      ier = nf_open(sic_file,nf_nowrite,sic_ncid)
      ier = nf_open(snod_file,nf_nowrite,snod_ncid)
      ier = nf_open(soiltemp_file,nf_nowrite,soiltemp_ncid)
      
      if (kl>lmax) then
        write(6,*) "ERROR: kl is greater than lmax"
        write(6,*) "kl= ",kl," lmax=",lmax
        call finishbanner
        stop -1
      end if

      call comsigalloc(kl)
      dsgx = dsg(1:kl)
      sgmlx = sgml(1:kl)

      spline = splineu .or. splinev .or. splinet

! set up what sigma levels the outgoing data will have
! assumes top down, ie. sg(1)=0., dsg>0

      sgx(1)=0.
      if ( sgmlx(kl/2).gt.0. ) then
! dsg=0, sgml>0
         do l=2,kl
           sgx(l)=.5*(sgmlx(l-1)+sgmlx(l))
         end do ! l=2,kl
         do l=2,kl+1
           dsgx(l-1)=sgx(l)-sgx(l-1)
         end do ! l=2,kl+1
      elseif ( dsgx(kl/2).gt.0. ) then
! sgml=0, dsg>0
         do l=2,kl-1
           sgx(l)=sgx(l-1)+dsgx(l-1)
         end do ! l=2,kl-1
         do l=1,kl
           sgmlx(l)=.5*(sgx(l)+sgx(l+1))
         end do ! l=1,kl
      elseif ( kl.eq.35 ) then
         call calcsig(sgx,kl)
         do l=1,kl
           sgmlx(l)=.5*(sgx(l)+sgx(l+1))
           dsgx(l)=sgx(l+1)-sgx(l)
         end do ! l=1,kl
      elseif ( oesig ) then
         do l=1,kl
           dsgx(l)=1./float(kl)
           sgmlx(l)=(l-.5)/float(kl)
           sgx(l+1)=sgx(l)+dsgx(l)
         end do ! l=1,kl
      else
         write(6,*)"Wrong sigma specification: STOP"
         call finishbanner
         stop -1
      endif
      sgx(kl+1)=1.

!####################### read topography data ############################
      write(6,*)'open ',inzs,' zsfil=',zsfil

      write(6,*)"set up cc geometry"
        
      liernc=nf_open(zsfil,nf_nowrite,lncid)
      if (liernc==nf_noerr) then
        iernc=nf_inq_dimid(lncid,"longitude",dimid)
        iernc=nf_inq_dimlen(lncid,dimid,ilx)
        iernc=nf_inq_dimid(lncid,"latitude",dimid)
        iernc=nf_inq_dimlen(lncid,dimid,jlx)
        iernc=nf_get_att_real(lncid,nf_global,"lon0",rlong0)
        if (iernc/=0) then
          write(6,*) "ERROR reading lon0 ",iernc
          call finishbanner
          stop -1
        end if
        iernc=nf_get_att_real(lncid,nf_global,"lat0",rlat0)
        iernc=nf_get_att_real(lncid,nf_global,"schmidt",schmidt)
        ds=0.
        header=''
      else
        write(6,*) "ERROR: cdfvidar requires netcdf version of zsfil"
        call finishbanner
        stop -1
      end if
      write(6,*)"gbl mapproj=",ilx,jlx,rlong0,rlat0,schmidt,ds

      il=ilx
      jl=6*il
      ifull=il*jl

      call latlongalloc(il)
      call cllalloc(il)
      call sigdataalloc(il,kl)
      allocate(hgt(il,jl,maxplev),temp(il,jl,maxplev))
      allocate(u(il,jl,maxplev),v(il,jl,maxplev))
      allocate(rh(il,jl,maxplev))
      allocate(twodim(il,jl),soildim(il,jl,6))
      allocate(lsmg_m(ifull))
      allocate(zsg(ifull))
      allocate(x(ifull),y(ifull),z(ifull))
      allocate(ax(ifull),ay(ifull),az(ifull))
      allocate(bx(ifull),by(ifull),bz(ifull))
      allocate(validlevcc(ifull))
      validlevcc = 1.
        
      ! set-up CC grid
      ccdim(1)=il
      ccdim(2)=jl
      lonlat(1)=rlong0
      lonlat(2)=rlat0
      allocate(rlld(ccdim(1),ccdim(2),2),grid(ccdim(1),ccdim(2)))
      allocate(xyz(ccdim(1),ccdim(2),3),axyz(ccdim(1),ccdim(2),3))
      allocate(bxyz(ccdim(1),ccdim(2),3))        
      !call cgg2(rlld,grid,ccdim,lonlat,schmidt,dst,iin,iie,iis,iiw)
      Call getcc(rlld,grid,xyz,axyz,bxyz,ccdim,lonlat,schmidt,dst)
      do i=1,ccdim(1)
        do j=1,ccdim(2)
          iq=i+(j-1)*ccdim(1)
          rlong(iq)=rlld(i,j,1)
          rlat(iq)=rlld(i,j,2)
          x(iq)=xyz(i,j,1)
          y(iq)=xyz(i,j,2)
          z(iq)=xyz(i,j,3)
          ax(iq)=axyz(i,j,1)
          ay(iq)=axyz(i,j,2)
          az(iq)=axyz(i,j,3)
          bx(iq)=bxyz(i,j,1)
          by(iq)=bxyz(i,j,2)
          bz(iq)=bxyz(i,j,3)
        end do
      end do
      deallocate(xyz,axyz,bxyz)
      deallocate (rlld,grid)
      !call setxyz

      rlatx=-1.e29
      rlatn= 1.e29
      rlonx=-1.e29
      rlonn= 1.e29

      do iq=1,ifull

!     convert conformal cubic lats & longs to degrees (-90 to 90) & (0 to 360)
!     used in sint16; N.B. original rlong is -pi to pi
        if(rlong(iq).lt.0.)rlong(iq)=rlong(iq)+360.
        if(rlat(iq).gt.rlatx)then
          rlatx=rlat(iq)
          ilatx=iq
        endif
        if(rlong(iq).gt.rlonx)then
          rlonx=rlong(iq)
          ilonx=iq
        endif
        if(rlat(iq).lt.rlatn)then
          rlatn=rlat(iq)
          ilatn=iq
        endif
        if(rlong(iq).lt.rlonn)then
          rlonn=rlong(iq)
          ilonn=iq
        endif

      enddo  ! iq loop

      write(6,*)"rlong,rlat(1,1)=",rlong(1),rlat(1)
      write(6,*)"rlong:x,n=",rlonx,ilonx,rlonn,ilonn
      write(6,*)"rlatg:x,n=",rlatx,ilatx,rlatn,ilatn
        
      write(6,*)'read model grid zs'
      iernc=nf_inq_varid(lncid,"zs",varid)
      call netcdferror(iernc)
      iernc=nf_get_var_real(lncid,varid,zs)
      call netcdferror(iernc)
      write(6,*) "zs max,min ",maxval(zs),minval(zs)
      write(6,*)'read model grid land-sea mask (0=ocean, 1=land)'
      iernc=nf_inq_varid(lncid,"lsm",varid)
      call netcdferror(iernc)
      iernc=nf_get_var_real(lncid,varid,lsm_m)
      call netcdferror(iernc)
      write(6,*) "lsm_m max,min ",maxval(lsm_m),minval(lsm_m)
      iernc=nf_close(lncid)
              
      ijd=id+il*(jd-1)
      write(6,*)"ijd=",ijd," zs(m)=",zs(ijd)," lsm_m=",lsm_m(ijd)
!####################### read vertical data ############################

      call readpress(t_ncid,in_type,plev,plev_b,nplev,osig_in,orev)
      
      write(6,*)"input nplev=",nplev
      write(6,*)"plevs=",(plev(k),k=1,nplev)
      

!########## get number of times in input netcdf file ###########

      ier = nf_inq_dimid(t_ncid,'time',idtim)
      write(6,*)"ier=",ier," idtim=",idtim

      ier= nf_inq_dimlen(t_ncid,idtim,narch)
      write(6,*)"ier=",ier," narch=",narch
      narch=min(narch,ntimes)

      ier = nf_inq_varid(t_ncid,'time',ivtim)
      write(6,*)"ier=",ier," ivtim=",ivtim

      ier = nf_get_att_text(t_ncid,ivtim,'units',timorg)
      write(6,*)"ier=",ier," timorg=",trim(timorg)
      !! new code to handle date time info (like onthefly)
      !if (ier.ne.0) then
      !  ier = nf_get_att_text(t_ncid,ivtim,'time_origin',timorg)
      !  write(6,*)"ier=",ier," timorg=",timorg
      !  if (ier.ne.0) then
      !     write(6,*)"cannot find valid timorg"
      !     call finishbanner
      !     stop -1
      !  endif
      !end if
      calendar=""
      ier = nf_get_att_text(t_ncid,ivtim,'calendar',calendar)

      i=scan(timorg,' ')-1
      cu=''  ! clear string to ensure blank
      cu(1:i)=timorg(1:i)
      if ( cu(1:i) == "since" ) then
        cu="hours"
      endif

      call processdatestring(timorg,iyr,imn,idy,ihr,imi)

      write(6,'("iyr,imn,idy,ihr,imi=",5i4)')iyr,imn,idy,ihr,imi
      write(6,*)"cu=",cu
      
      write(6,*)"++++++++++++++++++++++++++++++++++++++++++++++++++++++"

      write(6,*)" nplev=",nplev
      write(6,*)" inlsavn=",inlsavn

      write(6,*)' reading variables '

!***********************************************************************
      do iarch=1,narch
!***********************************************************************

      sdiag=.false.

      ier = nf_inq_varid(t_ncid,'time',ivtim)
      start = iarch
      ier = nf_get_var1_double(t_ncid,ivtim,start,time)
      nt=1
      
      select case(cu) ! MJT quick fix
        case('days')
          time=time*1440._8 
        case('hours')
          time=time*60._8 
        case('minutes')
          ! no change	
        case DEFAULT
          write(6,*) "cannot convert unknown time unit ",trim(cu)
          call finishbanner
          stop -1
      end select

      write(6,*)"time=",time
      
      kdate = iyr*10000 + imn*100 + idy
      ktime = ihr*100 + imi
      newtime = time
      call datefix(kdate,ktime,newtime,calendar)
      
      write(6,*)" input levels are bottom-up"
      write(6,*)" model levels in vidar are top-down"
      write(6,*)" nplev=",nplev

      write(6,*)"==================================================hgt"
      
      geopotunits="gpm"
      call readvar(z_ncid,"hgt",kdate,ktime,iarch,sdiag,in_type,plev(1:nplev),hgt(:,:,1:nplev),ier,       &
                   units=geopotunits)
      if ( ier/=nf_noerr ) then
        call readvar(z_ncid,"z",kdate,ktime,iarch,sdiag,in_type,plev(1:nplev),hgt(:,:,1:nplev),ier,       &
                     units=geopotunits)        
      end if
      if ( ier/=nf_noerr ) then
        call readvar(z_ncid,"geop_ht",kdate,ktime,iarch,sdiag,in_type,plev(1:nplev),hgt(:,:,1:nplev),ier, &
                     units=geopotunits)
      end if
      if ( ier==nf_noerr ) then
        if ( geopotunits=="m**2 s**-2" ) then
          write(6,*) "Converting from m**2 s**-2 to gpm"
          hgt(:,:,1:nplev) = hgt(:,:,1:nplev)/g
        end if
      end if

      write(6,*)"==================================================u"

      call readvar(u_ncid,"u",kdate,ktime,iarch,sdiag,in_type,plev(1:nplev),u(:,:,1:nplev),ier)
      if ( ier/=nf_noerr ) then
        call readvar(u_ncid,"ua",kdate,ktime,iarch,sdiag,in_type,plev(1:nplev),u(:,:,1:nplev),ier)        
      end if
      if ( ier/=nf_noerr ) then
        call readvar(u_ncid,"zonal_wnd",kdate,ktime,iarch,sdiag,in_type,plev(1:nplev),u(:,:,1:nplev),ier)        
      end if
      if ( ier/=nf_noerr ) then
        write(6,*) "ERROR: Cannot read u"
        call finishbanner
        stop -1
      end if  

      write(6,*)"==================================================v"

      call readvar(v_ncid,"v",kdate,ktime,iarch,sdiag,in_type,plev(1:nplev),v(:,:,1:nplev),ier)
      if ( ier/=nf_noerr ) then
        call readvar(v_ncid,"va",kdate,ktime,iarch,sdiag,in_type,plev(1:nplev),v(:,:,1:nplev),ier)        
      end if
      if ( ier/=nf_noerr ) then
        call readvar(v_ncid,"merid_wnd",kdate,ktime,iarch,sdiag,in_type,plev(1:nplev),v(:,:,1:nplev),ier)        
      end if
      if ( ier/=nf_noerr ) then
        write(6,*) "ERROR: Cannot read v"
        call finishbanner
        stop -1
      end if

      write(6,*)"==================================================temp"

      call readvar(t_ncid,"temp",kdate,ktime,iarch,sdiag,in_type,plev(1:nplev),temp(:,:,1:nplev),ier)
      if ( ier/=nf_noerr ) then
        call readvar(t_ncid,"ta",kdate,ktime,iarch,sdiag,in_type,plev(1:nplev),temp(:,:,1:nplev),ier)        
      end if
      if ( ier/=nf_noerr ) then
        call readvar(t_ncid,"air_temp",kdate,ktime,iarch,sdiag,in_type,plev(1:nplev),temp(:,:,1:nplev),ier)        
      end if
      if ( ier/=nf_noerr ) then
        write(6,*) "ERROR: Cannot read temp"
        call finishbanner
        stop -1
      end if
      
      write(6,*)"================================================rh/q"

      call readvar(rh_ncid,"rh",kdate,ktime,iarch,sdiag,in_type,plev(1:nplev),rh(:,:,1:nplev),ier)
      moist_var="rh"
      if ( ier/=nf_noerr ) then
         call readvar(rh_ncid,"relhum",kdate,ktime,iarch,sdiag,in_type,plev(1:nplev),rh(:,:,1:nplev),ier)
         moist_var="rh"
      end if
      if ( ier/=nf_noerr ) then
         call readvar(rh_ncid,"mix_rto",kdate,ktime,iarch,sdiag,in_type,plev(1:nplev),rh(:,:,1:nplev),ier)
         moist_var="mr"
      endif
      if ( ier/=nf_noerr ) then
         call readvar(rh_ncid,"hus",kdate,ktime,iarch,sdiag,in_type,plev(1:nplev),rh(:,:,1:nplev),ier)
         moist_var="mr"
         write(6,*) "Convert hus to mr"
         rh(:,:,1:nplev) = rh(:,:,1:nplev)/(1.-rh(:,:,1:nplev)) ! convert from specific humidity to mixing ratio
      endif
      if ( ier/=nf_noerr ) then
        write(6,*) "ERROR: Cannot read rh/q"
        call finishbanner
        stop -1
      end if

      write(6,*)"##################################moist_var=",moist_var

      xa = maxval(rh(:,:,1:nplev))
      if ( xa<.1 ) then
        moist_var="mr"
        write(6,*)"################################moist_var=",moist_var
      endif ! ( xa .lt. 1.1 ) then

      rh(:,:,1:nplev)=max(0.,min(100.,rh(:,:,1:nplev)))
      
      if ( moist_var .eq. "rh" .and. xa<1.1 ) then
        write(6,*)"######################convert rh from 0-1 to 0-100"
        rh(:,:,1:nplev)=max(0.,min(100.,rh(:,:,1:nplev)*100.))
      endif ! ( moist_var .eq. "rh" .and. xa .lt. 1.1 ) then

      call findxn(hgt(:,:, 1),ifull,-1.e29,xa,kx,an,kn)
      call findxn(hgt(:,:,nplev),ifull,-1.e29,xa,kx,an,kn)
      write(6,*)"nplev=",nplev     
      
!############################################################################
! sfc data
!############################################################################

      write(6,*)"================================================mslp"

      call readvar(psl_ncid,'mslp',kdate,ktime,iarch,sdiag,twodim,ier)
      if ( ier/=nf_noerr ) then
        call readvar(psl_ncid,'psl',kdate,ktime,iarch,sdiag,twodim,ier)  
      end if
      if ( ier/=nf_noerr ) then
        call readvar(psl_ncid,'pmsl',kdate,ktime,iarch,sdiag,twodim,ier)  
      end if
      if ( ier==nf_noerr ) then
        pmsl = reshape( twodim, (/ ifull /) )
        an = minval(pmsl)
        if ( an > 2000. ) then
           write(6,*)"#########################convert pmsl to hPa"
           pmsl(1:ifull)=pmsl(1:ifull)/100. ! to convert to hPa
        endif ! ( an .gt. 2000. ) then
      else
        write(6,*)"No pmsl data found, setting to 0"
        pmsl(1:ifull) = 0.
      endif ! ier

      call prt_pan(pmsl,il,jl,2,'pmsl')

      write(6,*)"================================================zs"

      geopotunits="gpm"
      call readvar(zs_ncid,'zs',sdiag,twodim,ier,units=geopotunits)
      if ( ier/=nf_noerr .and. inf==' ' ) then
        call readvar(zs_ncid,'z',sdiag,twodim,ier,units=geopotunits)  
      end if
      if ( ier/=nf_noerr ) then
        call readvar(zs_ncid,'orog',sdiag,twodim,ier,units=geopotunits)  
      end if
      if ( ier/=nf_noerr ) then
        call readvar(zs_ncid,'zsfc',sdiag,twodim,ier,units=geopotunits)  
      end if
      if ( ier/=nf_noerr ) then
        call readvar(zs_ncid,'topo',sdiag,twodim,ier,units=geopotunits)  
      end if
      if ( ier/=nf_noerr ) then
        call readvar(zs_ncid,'topog',sdiag,twodim,ier,units=geopotunits)  
      end if
      if ( ier==nf_noerr ) then
        zsi_m = reshape( twodim, (/ ifull /) )  
        if ( geopotunits=="m**2 s**-2" ) then
          write(6,*) "Converting from m**2 s**-2 to gpm"
          zsi_m(:) = zsi_m(:)/g
        end if
      else
        write(6,*)"No zs data found, setting to -999."
        zsi_m(1:ifull) = -999.
      endif ! ier

      call prt_pan(zs,il,jl,2,'zs(m)')
      call prt_pan(zsi_m,il,jl,2,'zsi_m(m)')
      
      write(6,*)"================================================land"
      write(6,*)"===================================== 1=land 0=ocean"

      call readvar(lsm_ncid,'land',sdiag,twodim,ier,island=.true.)
      if ( ier/=nf_noerr ) then
        call readvar(lsm_ncid,'lsm',sdiag,twodim,ier,island=.true.) 
      endif
      if ( ier/=nf_noerr ) then
        call readvar(lsm_ncid,'sftlf',sdiag,twodim,ier,island=.true.) 
      endif
      if ( ier/=nf_noerr ) then
        call readvar(lsm_ncid,'sfc_lsm',sdiag,twodim,ier,island=.true.) 
      endif
      if ( ier/=nf_noerr ) then
        call readvar(lsm_ncid,'land_mask',sdiag,twodim,ier,island=.true.) 
      endif
      if ( ier/=nf_noerr ) then
        call readvar(lsm_ncid,'lnd_mask',sdiag,twodim,ier,island=.true.) 
      endif
      if ( ier==nf_noerr ) then
        olsm_gbl = .true.
        lsmg_m = reshape( twodim, (/ ifull /) )
      else
        write(6,*)"No landmask data found, setting to -999."
        lsmg_m(1:ifull)=-999.
      endif ! ier

      call prt_pan(lsmg_m,il,jl,2,'lsmg_m')

      write(6,*)"================================================ps"

      call readvar(ps_ncid,'ps',kdate,ktime,iarch,sdiag,twodim,ier)
      if ( ier/=nf_noerr ) then
        call readvar(ps_ncid,'sfc_pres',kdate,ktime,iarch,sdiag,twodim,ier)  
      endif
      if ( ier==nf_noerr ) then
         psg_m = reshape( twodim, (/ ifull /) )
         if (any(psg_m>2000)) then
           psg_m=psg_m/100.
         end if
         ! remove levels below surface
         if ( .not.osig_in ) then
           do iq = 1,size(psg_m)
             do k = 1,nplev
               if ( psg_m(iq)>plev(k) ) then
                 validlevcc(iq) = real(k)
                 exit
               end if
             end do
           end do 
         end if  
         psg_m(1:ifull) = psg_m(1:ifull)*100. ! convert to Pa
      else
         write(6,*)"No sfcp data found, setting to -999."
         psg_m(1:ifull)=-999.
      end if ! ier

      call prt_pan(psg_m,il,jl,2,'psg_m')

      write(6,*)"================================================tss"

      call readsst(ts_ncid,lsm_ncid,'tss',kdate,ktime,iarch,sdiag,lsm_m,twodim,0,ier)
      if ( ier/=nf_noerr ) then
        call readsst(ts_ncid,lsm_ncid,'tos',kdate,ktime,iarch,sdiag,lsm_m,twodim,0,ier)  
      end if
      if ( ier/=nf_noerr ) then
        call readsst(ts_ncid,lsm_ncid,'sfc_temp',kdate,ktime,iarch,sdiag,lsm_m,twodim,0,ier)  
      end if
      if ( ier/=nf_noerr ) then
        call readsst(ts_ncid,lsm_ncid,'sst',kdate,ktime,iarch,sdiag,lsm_m,twodim,0,ier)  
      end if  
      if ( ier/=nf_noerr ) then
        call readsst(ts_ncid,lsm_ncid,'skt',kdate,ktime,iarch,sdiag,lsm_m,twodim,0,ier)
      end if 
      if ( ier==nf_noerr ) then
        sfct = reshape( twodim, (/ ifull /) )
        sfct = min( max( sfct, 100. ), 425. )
      end if  

      call prt_pan(sfct,il,jl,2,'tss')
      
      write(6,*)"============================================fracice"

      fracice=-1.
      call readsst(sic_ncid,lsm_ncid,'fracice',kdate,ktime,iarch,sdiag,lsm_m,twodim,1,ier)
      if ( ier/=nf_noerr ) then
         call readsst(sic_ncid,lsm_ncid,'sic',kdate,ktime,iarch,sdiag,lsm_m,twodim,1,ier) 
      endif
      if ( ier/=nf_noerr ) then
         call readsst(sic_ncid,lsm_ncid,'seaice',kdate,ktime,iarch,sdiag,lsm_m,twodim,1,ier) 
      endif
      if ( ier/=nf_noerr ) then
         call readsst(sic_ncid,lsm_ncid,'siconc',kdate,ktime,iarch,sdiag,lsm_m,twodim,1,ier)
      endif
      if ( ier==nf_noerr ) then
        fracice = reshape( twodim, (/ ifull /) )
      end if 
      if ( any(fracice>1.) ) then
        fracice = fracice/100.
      end if

      if ( ier==nf_noerr ) then ! we have fracice data
        where (lsm_m>=0.5)
          fracice=0.
        end where
      else
        write(6,*)"###############no fracice data in input dataset"
        write(6,*)"###############setting fracice data to -1.!!!!!!!"
        fracice(1:ifull) = -1.
      endif ! ier eq 0 , sfct

      call prt_pan(fracice,il,jl,2,'fracice')
      
      write(6,*)"================================================snod"

      snod = 0.
      call readvar(snod_ncid,'snod',kdate,ktime,iarch,sdiag,twodim,ier)
      if ( ier/=nf_noerr ) then
         call readvar(snod_ncid,'snow_amt_lnd',kdate,ktime,iarch,sdiag,twodim,ier) 
      endif
      if ( ier==nf_noerr ) then
         snod = reshape( twodim, (/ ifull /) )
         ! convert to equiv water
         snod(:) = snod(:)*100./1000. 
         where ( lsm_m(:)<0.5 )
           snod(:)=0.
         end where
      else
         write(6,*)"No snod data found, setting to 0."
         snod(:)=0.
      endif ! ier
      
      call prt_pan(snod,il,jl,2,'snod')
      
!############################################################################
! end sfc data
!############################################################################

!############################################################################
! soil data
!############################################################################
      
      write(6,*)"============================================soil_temp"

      soiltemp(:,:)=-1.
      call readsoil(soiltemp_ncid,"soil_temp",kdate,ktime,iarch,sdiag,soildepth_ccam,soildim,ier)
      if ( ier==nf_noerr ) then
        soiltemp = reshape( soildim, (/ ifull, 6 /) )  
      else
        write(6,*) "soil_temp not found"  
      end if
      
      write(6,*)"============================================soil_mois"

      soilmoist(:,:)=-1.
      call readsoil(soilmois_ncid,"soil_mois",kdate,ktime,iarch,sdiag,soildepth_ccam,soildim,ier) 
      if ( ier==nf_noerr ) then
        soilmoist = reshape( soildim, (/ ifull, 6 /) )
      else
        write(6,*) "soil_mois not found"  
      end if
      
!############################################################################
! end soil data
!############################################################################
      
      write(6,*)"check of temp data to ensure all is going okay"
      write(6,*)" findxn model temp(1)"
      call findxn(temp(:,:,1),ifull,-1.e29,xa,kx,an,kn)
      write(6,*)" findxn model temp(nplev)"
      call findxn(temp(:,:,nplev),ifull,-1.e29,xa,kx,an,kn)

      write(6,*)"nplev=",nplev
! constrain rh to 0-100
        do k=1,nplev
         do j=1,jl
          do i=1,il
           rh(i,j,k)=min(100.,max(0.,rh(i,j,k)))
          enddo ! i
         enddo ! j
        enddo ! k

!############### fix winds if CC grid ###############################

!     here use unstaggered lats and lons for u and v
!     For calculating zonal and meridional wind components, use the
!     following information, where theta is the angle between the
!     (ax,ay,az) vector [along the xg axis] and the zonal-component-vector:
!     veczon = k x r, i.e. (-y,x,0)/sqrt(x**2 + y**2)
!     vecmer = r x veczon, i.e. (-xz,-yz,x**2 + y**2)/sqrt(x**2 + y**2)
!     costh is (veczon . a) = (-y*ax + x*ay)/sqrt(x**2 + y**2)
!     sinth is (vecmer . a) = [-xz*ax - yz*ay + (x**2 + y**2)*az]/sqrt
!      using (r . a)=0, sinth collapses to az/sqrt(x**2 + y**2)
!     For rotated coordinated version, see JMcG's notes

      coslong=cos(rlong0*pi/180.)
      sinlong=sin(rlong0*pi/180.)
      coslat=cos(rlat0*pi/180.)
      sinlat=sin(rlat0*pi/180.)
      polenx=-coslat
      poleny=0.
      polenz=sinlat

      write(6,*)'polenx,poleny,polenz ',polenx,poleny,polenz
      write(6,*)'x(1),y(1),z(1)',x(1),y(1),z(1)
      write(6,*)'ax(1),ay(1),az(1)',ax(1),ay(1),az(1)
      write(6,*)'bx(1),by(1),bz(1)',bx(1),by(1),bz(1)

      write(6,*)'before zon/meridional'

      call maxmin(u,' u',0,1.,il,nplev)
      call maxmin(v,' v',0,1.,il,nplev)

      write(6,*)"convert winds to CCAM grid convention"

      cx=-1.e29
      sx=-1.e29
      cn= 1.e29
      sn= 1.e29


      do j=1,jl
        do i=1,il
          iq=i+(j-1)*il
!         set up unit zonal vector components
          zonx=            -polenz*y(iq)
          zony=polenz*x(iq)-polenx*z(iq)
          zonz=polenx*y(iq)
          den=sqrt( max(zonx**2 + zony**2 + zonz**2,1.e-7) )  ! allow for poles
          costh= (zonx*ax(iq)+zony*ay(iq)+zonz*az(iq))/den
          sinth=-(zonx*bx(iq)+zony*by(iq)+zonz*bz(iq))/den
          cx=max(cx,costh)
          sx=max(sx,sinth)
          cn=min(cn,costh)
          sn=min(sn,sinth)
          do k=1,nplev
             uzon = u(i,j,k)
             vmer = v(i,j,k)
             u(i,j,k)= costh*uzon+sinth*vmer
             v(i,j,k)=-sinth*uzon+costh*vmer
          end do  ! k loop
        end do
      end do

      write(6,*)'cx,cn,sx,sn=',cx,cn,sx,sn
      write(6,*)'after zon/meridional'

      call maxmin(u(:,:,1:nplev),' u',0,1.,il,nplev)
      call maxmin(v(:,:,1:nplev),' v',0,1.,il,nplev)
      call maxmin(hgt(:,:,1:nplev),' hgt',0,.001,il,nplev)
      call maxmin(temp(:,:,1:nplev),' temp',0,1.,il,nplev)
      call maxmin(rh(:,:,1:nplev),' rh',0,1.,il,nplev)


      i=il/2
      j=jl/2
      write(6,'(6a10," at i,j=",2i3)') "p","z","t","u","v","rh",i,j
      write(6,*)"invert order of plev  nplev=",nplev
      write(6,'(6a10)')"plev","hgt","temp","u","v","rh/mr"
      xplev = maxval( plev(1:nplev) )
      tohpa = 1.
      if ( xplev > 1500. ) then
        tohpa = 0.01
      end if
      do k=1,nplev
          cplev(k)=tohpa*plev(k)
          bplev(k)=tohpa*plev_b(k)
          write(6,'(5f10.2,f10.5)') cplev(k),hgt(i,j,k),temp(i,j,k),u(i,j,k),v(i,j,k),rh(i,j,k)
      enddo ! k=1,nplev

      iq=il/2+(jl/2-1)*il
      write(6,'("pmsl=",f12.2," sfct=",f12.2)') pmsl(iq),sfct(iq)

      write(6,*) "validlevcc max,min ",maxval(validlevcc),minval(validlevcc)

      write(6,*)"calling vidar now!! ntimes,iarch=",ntimes,iarch

      if2=0

!#######################################################################
      call vidar(nplev,hgt,temp,u,v,rh,validlevcc,iyr,imn,idy,ihr,iarch,time,mtimer,cplev,bplev, &
                 io_out,il,kl,minlon,maxlon,minlat,maxlat,llrng,moist_var,calendar,rlong0,rlat0, &
                 schmidt)
!#######################################################################

      enddo ! narch

      ! close files
      ier = nf_close(t_ncid)
      ier = nf_close(rh_ncid)
      ier = nf_close(u_ncid)
      ier = nf_close(v_ncid)
      ier = nf_close(z_ncid)
      ier = nf_close(lsm_ncid)
      ier = nf_close(zs_ncid)
      ier = nf_close(ps_ncid)
      ier = nf_close(psl_ncid)
      ier = nf_close(ts_ncid)
      ier = nf_close(sic_ncid)
      ier = nf_close(snod_ncid)
      ier = nf_close(soiltemp_ncid)
      ier = nf_close(soilmois_ncid)
      
      !deallocate(datan,zs_gbl,lsm_gbl)
      !deallocate(glon,glat)
      !deallocate(validlevhost)

      write(6,*)'*********** Finished cdfvidar ************************'
      
      deallocate(x,y,z,ax,ay,az,bx,by,bz)
      deallocate(hgt,temp,u,v,rh)
      deallocate(twodim,soildim)
      deallocate(lsmg_m)
      deallocate(zsg)
      call sigdatadealloc    
      call latlongdealloc      
      call clldealloc
      call comsigdealloc

      ! Complete
      write(6,*) "CCAM: cdfvidar completed successfully"
      
      ! End banner
      call finishbanner
      
      stop
      end ! cdfvidar
!***************************************************************************

!***********************************************************************
      function icmonth_to_imn(cmonth)

      integer icmonth_to_imn
      character*(*) cmonth

      write(6,*)"icmonth_to_imn cmonth=",cmonth

      icmonth_to_imn=0
      if ( cmonth.eq.'jan' ) icmonth_to_imn=1
      if ( cmonth.eq.'feb' ) icmonth_to_imn=2
      if ( cmonth.eq.'mar' ) icmonth_to_imn=3
      if ( cmonth.eq.'apr' ) icmonth_to_imn=4
      if ( cmonth.eq.'may' ) icmonth_to_imn=5
      if ( cmonth.eq.'jun' ) icmonth_to_imn=6
      if ( cmonth.eq.'jul' ) icmonth_to_imn=7
      if ( cmonth.eq.'aug' ) icmonth_to_imn=8
      if ( cmonth.eq.'sep' ) icmonth_to_imn=9
      if ( cmonth.eq.'oct' ) icmonth_to_imn=10
      if ( cmonth.eq.'nov' ) icmonth_to_imn=11
      if ( cmonth.eq.'dec' ) icmonth_to_imn=12

      write(6,*)"icmonth_to_imn=",icmonth_to_imn

      return 
      end
!***********************************************************************
      
      subroutine filldat(datan,ix,iy,plev)
      
      implicit none
      
      integer, intent(in) :: ix,iy,plev
      integer i,j,k,is,ie,iqk,iqn
      integer ncount
      real, dimension(ix*iy*plev), intent(inout) :: datan
      real, dimension(ix*iy) :: datatemp
      real datasum

      if ( any(datan(:)>9.e9) ) then
        write(6,*) "Using filldat to remove missing values"
      end if
      
      do k = 1,plev
          
        is = ix*iy*(k-1) + 1
        ie = ix*iy*k
        if (all(datan(is:ie)>9.E9)) then
          write(6,*) "ERROR: No valid data on level"
          call finishbanner
          stop -1
        end if
        
        datatemp(:) = datan(is:ie)  
        
        do while ( any(datatemp(:)>9.E9) )
          
          do j = 1,iy
            do i = 1,ix
                
              iqk = i + ix*(j-1) ! no vertical dimension for datatemp
              if (datatemp(iqk)>9.e9) then
                ncount = 0
                datasum = 0.
                if ( j<iy ) then
                  iqn = i + ix*j + ix*iy*(k-1)
                  if ( datan(iqn)<=9.e9 ) then
                    datasum = datasum + datan(iqn)
                    ncount = ncount + 1
                  end if
                end if
                if ( j>1 ) then
                  iqn = i + ix*(j-2) + ix*iy*(k-1)
                  if ( datan(iqn)<=9.e9 ) then
                    datasum = datasum + datan(iqn)
                    ncount = ncount + 1
                  end if
                end if
                if ( i<ix ) then
                  iqn = i + 1 + ix*(j-1) + ix*iy*(k-1)
                else
                  iqn = 1 + ix*(j-1) + ix*iy*(k-1)	
                end if
                if ( datan(iqn)<=9.e9 ) then
                  datasum = datasum + datan(iqn)
                  ncount = ncount + 1
                end if
                if ( i>1 ) then
                  iqn = i - 1 + ix*(j-1) + ix*iy*(k-1)
                else
                  iqn = ix + ix*(j-1) + ix*iy*(k-1)
                end if
                if ( datan(iqn)<=9.e9 ) then
                  datasum = datasum + datan(iqn)
                  ncount = ncount + 1
                end if
                if ( ncount>0 ) then
                  datatemp(iqk) = datasum/real(ncount)
                end if
              end if
            end do
          end do
          
          datan(is:ie) = datatemp(:)
          
        end do
      end do
      
      return
      end subroutine filldat

subroutine processdatestring(datestring,yyyy,mm,dd,hh,mt)

implicit none

integer, intent(out) :: yyyy, mm, dd, hh, mt
integer iposa, iposb, ierx
character(len=*), intent(in) :: datestring

!if ( datestring(1:7)/='minutes' ) then
!  write(6,*) "ERROR: Time units expected to be minutes"
!  write(6,*) "Found ",trim(datestring)
!  stop -1
!end if

! process year
iposa = index(trim(datestring),'since')
iposa = iposa + 5 ! skip 'since'
iposb = index(trim(datestring(iposa:)),'-')
iposb = iposa + iposb - 2 ! remove '-'
read(datestring(iposa:iposb),FMT=*,iostat=ierx) yyyy
if ( ierx/=0 ) then
  write(6,*) "ERROR reading time units.  Expecting year but found ",datestring(iposa:iposb)
  call finishbanner
  stop -1
end if

! process month
iposa = iposb + 2 ! skip '-'
iposb = index(trim(datestring(iposa:)),'-')
iposb = iposa + iposb - 2 ! remove '-'
read(datestring(iposa:iposb),FMT=*,iostat=ierx) mm
if ( ierx/=0 ) then
  write(6,*) "ERROR reading time units.  Expecting month but found ",datestring(iposa:iposb)
  call finishbanner
  stop -1
end if

! process day
iposa = iposb + 2 ! skip '-'
iposb = index(trim(datestring(iposa:)),' ')
iposb = iposa + iposb - 2 ! remove ' '
read(datestring(iposa:iposb),FMT=*,iostat=ierx) dd
if ( ierx/=0 ) then
  write(6,*) "ERROR reading time units.  Expecting day but found ",datestring(iposa:iposb)
  call finishbanner
  stop -1
end if

! process hour
iposa = iposb + 2 ! skip ' '
iposb = index(trim(datestring(iposa:)),':')
iposb = iposa + iposb - 2 ! remove ':'
read(datestring(iposa:iposb),FMT=*,iostat=ierx) hh
if ( ierx/=0 ) then
  write(6,*) "ERROR reading time units.  Expecting hour but found ",datestring(iposa:iposb)
  call finishbanner
  stop -1
end if

! process mins
iposa = iposb + 2 ! skip ':'
iposb = index(trim(datestring(iposa:)),':')
iposb = iposa + iposb - 2 ! remove ':'
read(datestring(iposa:iposb),FMT=*,iostat=ierx) mt
if ( ierx/=0 ) then
  write(6,*) "ERROR reading time units.  Expecting minutes but found ",datestring(iposa:iposb)
  call finishbanner
  stop -1
end if

return
end subroutine processdatestring

subroutine finishbanner

implicit none

! End banner
write(6,*) "=============================================================================="
write(6,*) "CCAM: Finished cdfvidar"
write(6,*) "=============================================================================="

return
end

