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
    
      program cdfvidar

      use ccinterp
      use cll_m
      use comsig_m      
      use latlong_m
      use netcdf_m
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
      real time,xa,an
      real spval,sinlong,coslat
      real coslong,polenz,poleny
      real polenx,sinlat,sx,cx
      real zonx,zony,sn,cn
      real sinth,conth,den,zonz
      real uzon,costh,ull,vll
      real vmer,ther,ds,du,tanl
      real rnml,stl1,stl2
      real rlon
      real llrng,minlon,maxlon,minlat,maxlat
      real plev0

      character(len=1024) inf

      common/mapproj/du,tanl,rnml,stl1,stl2
      include 'nplevs.h' ! maxplev
      real, dimension(maxplev) :: plev,cplev,bplev
      real, dimension(maxplev) :: plev_b
      common/levpre/nplev,plev,plev_b ! does this need to be a common block?

      include 'lmax.h'

      parameter ( pi=3.1415926536 )
      parameter ( g=9.80616 )
      
      real, dimension(lmax) :: dsg,sgml

      common/ncdfids/dimil,dimjl,dimkl,dimtim,idil,idjl,idkl,idnt
! dimension ids
      integer  dimil,dimjl,dimkl,dimtim
! variable ids
      integer  idil,idjl,idkl,idnt,ix,iy
      integer  il,jl,kl,ifull
      integer ncid
      integer iernc,liernc,lncid,varid,dimid
      integer idsoillvl, ivsoillvl, nsoillvl, ksearch, ktest
      integer procformat_nproc
      integer, dimension(2) :: ccdim
      integer, dimension(2) :: start, ncount

      common/lconther/ther
      include 'vidar.h'
      logical sdiag

      real dst, xfrac
      real, dimension(2) :: lonlat
      real, dimension(:,:,:), allocatable :: rlld,xyz,axyz,bxyz      
      real, dimension(:,:), allocatable :: grid
      real, dimension(:), allocatable :: ax,ay,az,bx,by,bz,x,y,z
      real, dimension(:), allocatable :: datan, datatemp
      real, dimension(:), allocatable :: zs_gbl, lsm_gbl
      real, dimension(:), allocatable :: glon, glat

      real plevin(maxplev)
      real, dimension(:,:,:), allocatable :: hgt
      real, dimension(:,:,:), allocatable :: temp
      real, dimension(:,:,:), allocatable :: u
      real, dimension(:,:,:), allocatable :: v
      real, dimension(:,:,:), allocatable :: rh
      real, dimension(:), allocatable :: sfcto_m
      real, dimension(:), allocatable :: lsmg_m
      real, dimension(:), allocatable :: zsg
      real, dimension(:), allocatable :: validlevhost
      real, dimension(:), allocatable :: validlevcc
      real, dimension(:), allocatable :: soildepth_in
      real, dimension(6), parameter :: soildepth_ccam = (/ 0.011, 0.051, 0.157, 0.4385, 1.1855, 3.164 /)
      real, dimension(6), parameter :: soilthick_ccam = (/ 0.022, 0.058, 0.154, 0.409,  1.085,  2.872 /)
      real, parameter :: rhowater = 1000.
      real fill_float

      common/datatype/moist_var,in_type
      character*1 in_type
      character*2 moist_var

      logical ofirst, ogbl, orev, olsm_gbl
      character*60 timorg
      character*60 cu
      character*3 cmonth
      character*80 zsavn,lsavn
      character*10 header,moistvar
      character*10 soilunits, geopotunits, presunits
      character*80 presname

      namelist/gnml/inf,vfil,ds,du,tanl,rnml,stl1,stl2,inzs,zsfil   &
                   ,ints,tsfil, ogbl,zsavn,inzsavn,lsavn,inlsavn    &
                   ,plevin,orev,io_out,igd,jgd,id,jd,mtimer,ntimes  &
                   ,spline,mxcyc,nvsig,nrh                          &
                   ,oesig,sgml,dsg,ptop,debug,notop,opre,have_gp    &
                   ,in,calout                                       &
                   ,iout,oform,sdiag                                &
                   ,insm,smfil                                      &
                   ,splineu,splinev,splinet,zerowinds               &
                   ,grdx,grdy,slon,slat                             &
                   ,moistvar,kl,llrng,procformat_nproc

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
      data llrng/0./,procformat_nproc/0/

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

!####################### read namelists ############################
      write(6,*)'read namelist'
      read (5, gnml)
      write(6,nml=gnml)
!####################### read namelist ############################

      if (kl>lmax) then
        write(6,*) "ERROR: kl is greater than lmax"
        write(6,*) "kl= ",kl," lmax=",lmax
        call finishbanner
        stop -1
      end if

      call comsigalloc(kl)
      dsgx=dsg(1:kl)
      sgmlx=sgml(1:kl)

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

      if ( ogbl ) then
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
          open(unit=inzs,file=zsfil,status='old',form='formatted')
          read(inzs,*)ilx,jlx,rlong0,rlat0,schmidt,ds,header
        end if
        du=rlong0
        tanl=rlat0
        rnml=schmidt
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
        allocate(sfcto_m(ifull),lsmg_m(ifull))
        allocate(zsg(ifull))
        allocate(x(ifull),y(ifull),z(ifull))
        allocate(ax(ifull),ay(ifull),az(ifull))
        allocate(bx(ifull),by(ifull),bz(ifull))
        allocate(validlevcc(ifull))
        
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

!       convert conformal cubic lats & longs to degrees (-90 to 90) & (0 to 360)
!       used in sint16; N.B. original rlong is -pi to pi
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
      else ! ( not ogbl ) then
        liernc=-1
        open(inzs,file=zsfil,form='formatted',recl=il*7,status='old')
        write(6,*)'read zsfil header'
        read(inzs,*)ilt,jlk,ds,du,tanl,rnml,stl1,stl2
        if(ilt.eq.0.or.jlk.eq.0)then
           write(6,*)'no header in newtopo file'
        else
           write(6,*)'Header information for topofile'
           write(6,*)'ilt,jlk,ds,du,tanl,rnml,stl1,stl2',ilt,jlk,ds,du,tanl,rnml,stl1,stl2
           if(ilt/=il.or.jlk/=jl) then
             write(6,*) 'wrong topofile supplied'
             call finishbanner
             stop -1
           end if
        endif     ! (ilt.eq.0.or.jlk.eq.0)
        write(6,*)"set up model grid params by calling lconset ds=",ds
        call lconset(ds)
      endif ! ( ogbl ) then

      if (liernc==nf_noerr) then
        write(6,*)'read model grid zs'
        iernc=nf_inq_varid(lncid,"zs",varid)
        iernc=nf_get_var_real(lncid,varid,zs)
        write(6,*)'read model grid land-sea mask (0=ocean, 1=land)'
        iernc=nf_inq_varid(lncid,"lsm",varid)
        iernc=nf_get_var_real(lncid,varid,lsm_m)
        iernc=nf_close(lncid)
      else
        write(6,*)'read model grid zsg = g*zs'
        read(inzs,*)zsg
        write(6,*)'read model grid land-sea mask (0=ocean, 1=land)'
        read(inzs,*)lsm_m
        close(inzs)
        write(6,*)'convert g*zs to zs(m)'
        zs(1:ifull)=zsg(1:ifull)/g ! convert ascii read in zs*g to zs(m)
      end if
      
      ijd=id+il*(jd-1)
      write(6,*)"ijd=",ijd," zs(m)=",zs(ijd)," lsm_m=",lsm_m(ijd)
!####################### read topography data ############################

!####################### open input netcdf file ############################
      write(6,*)'inf='
      write(6,*) trim(inf)
      ier = nf_open(inf,nf_nowrite,ncid)
      write(6,*)'ncid=',ncid
      if(ier.ne.0) then
        write(6,*)' cannot open netCDF file; error code ',ier
        call finishbanner
        stop -1
      end if

!####################### get attributes of input netcdf file ############################
      ier = nf_inq(ncid,ndims,nvars,ngatts,irecd)
      write(6,'("ndims,nvars,ngatts,irecd,ier")')
      write(6,'(5i6)') ndims,nvars,ngatts,irecd,ier

! Get dimensions
      write(6,*) "get dim1 ncid=",ncid
      ier = nf_inq_dimid(ncid,'lon',lonid)
      write(6,*)"lon ncid,lonid,ier=",ncid,lonid,ier
      if ( ier.eq.0 ) then
        write(6,*)"ncid,lonid=",ncid,lonid
        ier= nf_inq_dimlen(ncid,lonid,ix)
        write(6,*)"input ix,ier=",ix,ier
        ier = nf_inq_dimid(ncid,'lat',latid)
        ier= nf_inq_dimlen(ncid,latid,iy)
        write(6,*)"input iy,ier=",iy,ier
        ier = nf_inq_varid(ncid,'lon',idv)

        allocate(glon(ix),glat(iy))

! get glon from input dataset
        ier = nf_get_var_real(ncid,idv,glon)
        ier = nf_inq_varid(ncid,'lat',idv)
        ier = nf_get_var_real(ncid,idv,glat)
      else
        write(6,*)"now try longitude"
        ier = nf_inq_dimid(ncid,'longitude',lonid)
        write(6,*)"lonid=",lonid," ier=",ier
        ier= nf_inq_dimlen(ncid,lonid,ix)
        write(6,*)"input ix=",ix," ier=",ier

        ier = nf_inq_dimid(ncid,'latitude',latid)
        ier= nf_inq_dimlen(ncid,latid,iy)
        write(6,*)"input iy=",iy

        allocate(glon(ix),glat(iy))

        ier = nf_inq_varid(ncid,'longitude',idv)
        ier = nf_get_var_real(ncid,idv,glon)
        write(6,*)"glon=",(glon(i),i=1,ix)
        ier = nf_inq_varid(ncid,'latitude',idv)
        ier = nf_get_var_real(ncid,idv,glat)
        write(6,*)"glat=",(glat(i),i=1,iy)

      endif ! ( ier .eq. 0 ) then

! find min/max  for input data set
      slon = glon(1)
      elon = glon(ix)
      slat = glat(1)
      elat = glat(iy)
      if ( slon<elon ) then
        minlon = glon(2)
        maxlon = glon(ix-2)
      else
        minlon = glon(ix-2)
        maxlon = glon(2)
      end if    
      if ( slat<elat ) then
        minlat = glat(2)
        maxlat = glat(iy-2)
      else
        minlat = glat(iy-2)
        maxlat = glat(2)
      end if
      write(6,*)"==================> slon=",slon," elon=",elon," ix=",ix
      write(6,*)"==================> slat=",slat," elat=",elat," iy=",iy

      ier = nf_inq_varid(ncid,'pres',ivpres)
      ier = nf_inq_dimid(ncid,'pres',idpres)
      in_type="p"
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'plev',ivpres) 
         ier = nf_inq_dimid(ncid,'plev',idpres)
         in_type="p"
      endif
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'lev',ivpres) 
         ier = nf_inq_dimid(ncid,'lev',idpres)
         in_type="p"
      endif
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'lvl',ivpres) 
         ier = nf_inq_dimid(ncid,'lvl',idpres)
         in_type="s"
      endif
      ier = nf_get_att_text(ncid,ivpres,'long_name',presname)
      if ( presname(1:32) == "hybrid sigma pressure coordinate" ) then
        in_type="h"  
      end if
      write(6,*)"ier=",ier," idpres=",idpres," in_type=",in_type
      
      ier= nf_inq_dimlen(ncid,idpres,nplev)
      write(6,*)"ier=",ier," nplev=",nplev

      write(6,*)"ier=",ier," ivpres=",ivpres

      ier = nf_get_var_real(ncid,ivpres,plev)
      write(6,*)"ier=",ier," ivpres=",ivpres
      write(6,*)"input nplev=",nplev
      write(6,*)"plevs=",(plev(k),k=1,nplev)

      ier = nf_get_att_text(ncid,ivpres,'units',presunits)
      write(6,*)"ier=",ier," presunits=",trim(presunits)     
      
      orev = plev(nplev).gt.plev(1)
      write(6,*)"#################################### orev=",orev

      ier = nf_inq_dimid(ncid,'soil_lvl',idsoillvl)
      write(6,*)"ier=",ier," idsoillvl=",idsoillvl
      if ( ier==0 ) then
        ier = nf_inq_dimlen(ncid,idsoillvl,nsoillvl)
        write(6,*)"ier=",ier," nsoillvl=",nsoillvl
      else
        nsoillvl = 0  
      end if
      if ( nsoillvl>0 ) then    
        allocate( soildepth_in(nsoillvl) )
        ier = nf_inq_varid(ncid,'soil_lvl',ivsoillvl)
        write(6,*)"ier=",ier," ivsoillvl=",ivsoillvl
        ier = nf_get_var_real(ncid,ivsoillvl,soildepth_in)
        write(6,*)"soildepth_in=",soildepth_in(1:nsoillvl)
      end if
      
      write(6,*)"####################################"
      
      write(6,*)"allocate arrays ix,iy,nplevs=",ix,iy,nplev
      if ( nsoillvl>nplev ) then
        allocate(datan(ix*iy*nsoillvl))
      else
        allocate(datan(ix*iy*nplev))
      end if
      allocate(datatemp(ix*iy))
      allocate(zs_gbl(ix*iy))
      allocate(lsm_gbl(ix*iy))
      allocate(validlevhost(ix*iy))
      validlevhost(:) = 1.
 
      if ( orev ) then
        do k = 1,nplev
          datan(k) = plev(k)
        end do
        do k = 1,nplev
          plev(k) = datan(nplev+1-k)
        end do
      end if

      xplev = maxval( plev(1:nplev) )
      write(6,*)"xplev=",xplev

      plev_b = 0.
      
      osig_in = .false.
      if ( in_type == "h" ) then
        write(6,*)"^^^^^^^^^hybrid sigma levels^^^^^^^^"
        osig_in = .true.
        ier = nf_inq_varid(ncid,'a',ivpres_a)
        ier = nf_get_var_real(ncid,ivpres_a,plev)
        ier = nf_inq_varid(ncid,'b',ivpres_b)
        ier = nf_get_var_real(ncid,ivpres_b,plev_b)
        ier = nf_inq_varid(ncid,'p0',ivpres_0)
        ier = nf_get_var_real(ncid,ivpres_0,plev0)
        plev0 = plev0/100. ! convert to hPa
        write(6,*) "a=",plev(1:nplev)
        write(6,*) "b=",plev_b(1:nplev)
        write(6,*) "p0=",plev0
        do k = 1,nplev
          plev_b(k) = plev_b(k)*plev0  
          plev(k) = plev(k)*1000. ! updated later  
        end do
        presunits="hPa"  
      else if ( .01<xplev .and. xplev<800.  ) then
        write(6,*)"^^^^^^^^^actualy sigma levels^^^^^^^ fix plevs"
        osig_in = .true.
        plev_b = 0.
        do k=1,nplev
          plev(k)=plev(k)*1000. !
        enddo
        presunits="hPa"
        write(6,*)"plevs=",(plev(k),k=1,nplev)
      else if ( xplev .le. .01  ) then
        write(6,*)"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ fix plevs"
        do k=1,nplev
          plev(k)=plevin(k)
        end do
        write(6,*)"plevs=",(plev(k),k=1,nplev)
        write(6,*) 'xplev < 0.01 in cdfvidar'
        call finishbanner
        stop -1
      end if

      if ( presunits=="Pa" ) then
        write(6,*) "Converting pressure levels from Pa to hPa"  
        plev = plev/100.
        presunits="hPa"
      end if
      
      if ( presunits/="hPa" ) then
        write(6,*) "ERROR: Could not convert vertical levels to hPa"
        write(6,*) "Vertical level units was read as ",trim(presunits)
        call finishbanner
        stop -1
      end if
      
      ier = nf_inq_dimid(ncid,'time',idtim)
      write(6,*)"ier=",ier," idtim=",idtim

!########## get number of times in input netcdf file ###########

      ier= nf_inq_dimlen(ncid,idtim,narch)
      write(6,*)"ier=",ier," narch=",narch
      narch=min(narch,ntimes)

      ier = nf_inq_varid(ncid,'time',ivtim)
      write(6,*)"ier=",ier," ivtim=",ivtim

      !call ncagtc(ncid,ivtim,"time_origin",timorg,20,ier)
      ier = nf_get_att_text(ncid,ivtim,'units',timorg) ! MJT quick fix
      write(6,*)"ier=",ier," timorg=",trim(timorg)


! new code to handle date time info (like onthefly)
      if (ier.ne.0) then
        ier = nf_get_att_text(ncid,ivtim,'time_origin',timorg)
        write(6,*)"ier=",ier," timorg=",timorg
        if (ier.ne.0) then
           write(6,*)"cannot find valid timorg"
           call finishbanner
           stop -1
        endif
      end if

      i=scan(timorg,' ')-1
      cu=''  ! clear string to ensure blank
      cu(1:i)=timorg(1:i)
      if ( cu(1:i) == "since" ) then
        cu="hours"
      endif

      call processdatestring(timorg,iyr,imn,idy,ihr,imi)

      write(6,'("iyr,imn,idy,ihr,imi=",5i4)')iyr,imn,idy,ihr,imi
      write(6,*)"cu=",cu

      do j=1,iy
       do i=1,ix
         datan(i+(j-1)*ix     )=glon(i)
         datan(i+(j-1)*ix+ix*iy)=glat(j)
       enddo ! i
      enddo ! j


! printout of glon
      do j=1,iy,iy-1
        do i=1,ix,ix-1
          write(6,*)i,j,datan(i+(j-1)*ix)
        enddo
      enddo

       call prt_pan(rlong,il,jl,2,'rlong')
       call prt_pan(rlat ,il,jl,2,'rlat')

      write(6,*)"============= sintp16 clon++++++++++++++++++++++++++++"
      write(6,*)" nplev=",nplev

      call sintp16(datan(1:ix*iy),ix,iy,clon,glon,glat,sdiag,il)

      call prt_pan(clon,il,jl,2,'clon')

! printout of glat
      do j=1,iy,iy-1
        do i=1,ix,ix-1
          write(6,*)i,j,datan(i+(j-1)*ix+ix*iy)
        enddo
      enddo

      write(6,*)"============= sintp16 clat++++++++++++++++++++++++++++"
      write(6,*)" nplev=",nplev

      call sintp16(datan(1+ix*iy:2*ix*iy),ix,iy,clat,glon,glat,sdiag,il)

      call prt_pan(clat,il,jl,2,'clat pan2')

      write(6,'("ix,iy,nplev,narch=",4i5)')ix,iy,nplev,narch

      write(6,*)"++++++++++++++++++++++++++++++++++++++++++++++++++++++"

      write(6,*)" nplev=",nplev
      write(6,*)" inlsavn=",inlsavn

      write(6,*)'read land-sea mask (0=ocean, 1=land) inlsavn=',inlsavn
      if ( inlsavn>0 ) then
        olsm_gbl=.true.
        write(6,*)"AVN land sea mask"
        write(6,*)"open sfc. orog. for avn data"
        open(inlsavn,file=lsavn,form='formatted',status='old')
!       rewind inlsavn
        write(6,*)"note that it runs north to south"
        read(inlsavn,*)((lsm_gbl(i+(j-1)*ix),i=1,ix),j=1,iy)
        close(inlsavn)
        call amap ( lsm_gbl, ix, iy, 'gbl lsmsk', 0., 0. )
      else
        olsm_gbl=.false.
        write(6,*)"######################## WARNING!!!!!!!!!!!!!!!!"
        write(6,*)"######################## since inlsavn le 0 ####"
        write(6,*)"######################## setting input lsm == 1!"
        write(6,*)" nplev=",nplev,ix,iy
        do j=1,iy
         do i=1,ix
          lsm_gbl(i+(j-1)*ix)=1.
         enddo ! i
        enddo ! j
      endif
      write(6,*)" nplev=",nplev

      write(6,*)"Now deal with sfc. zs inzsavn=",inzsavn
      if ( inzsavn>0 ) then
        write(6,*)"open sfc. orog. for avn data"
        open(inzsavn,file=zsavn,form='formatted',status='old')
! read in free formatted avn sfc. orography (28-mar-2000)
! note that it runs north to south
!       rewind inzsavn
        read(inzsavn,*)((zs_gbl(i+(j-1)*ix),i=1,ix),j=1,iy)
        call amap ( zs_gbl, ix, iy, 'gbl sfczs(m)', 0., 0. )
        write(6,*) 'close unit inzsavn=',inzsavn
        close(inzsavn)
      endif

      write(6,*)' reading variables '

!***********************************************************************
      do iarch=1,narch
!***********************************************************************

      sdiag=.false.

      ier = nf_inq_varid(ncid,'time',ivtim)
      start = iarch
      ier = nf_get_var1_real(ncid,ivtim,start,time)
      nt=1
      
      select case(cu) ! MJT quick fix
        case('days')
          time=time*1440. 
        case('hours')
          time=time*60. 
        case('minutes')
          ! no change	
        case DEFAULT
          write(6,*) "cannot convert unknown time unit ",trim(cu)
          call finishbanner
          stop -1
      end select

      write(6,*)"time=",time
      
      if ( time>1.E9 ) then
        write(6,*) "ERROR: Time too large for real variable"
        write(6,*) "Consider adjusting base date"
        call finishbanner
        stop -1
      end if

      write(6,*)" input levels are bottom-up"
      write(6,*)" model levels in vidar are top-down"
      write(6,*)" nplev=",nplev

      write(6,*)"==================================================hgt"

      ier = nf_inq_varid(ncid,'hgt',idvar)
      write(6,*)"ier=",ier," idvar=",idvar
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'geop_ht',idvar)
         write(6,*)"ier=",ier," idvar=",idvar
      endif

      if ( ier .eq. 0 ) then

        write(6,*) ncid,iarch,idvar,ix,iy,nplev
        call ncread_3d(ncid,iarch,idvar,ix,iy,nplev,datan)

        ! MJT quick fix
        ier = nf_get_att_real(ncid,idvar,'_FillValue',fill_float)
        if ( ier .ne. 0 ) then
          ier = nf_get_att_real(ncid,idvar,'missing_value',fill_float)  
        end if
        if ( ier .eq. 0 ) then
          where ( datan==fill_float )
            datan = 1.e10
          end where
        end if
        call getvalidlev(validlevhost,datan,ix,iy,nplev)
        call filldat(datan,ix,iy,nplev)

        call amap (datan(1:ix*iy),ix,iy,'input hgt',0.,0.)
        call amap (datan(1+ix*iy*(nplev-1):ix*iy*nplev),ix,iy,'input hgt',0.,0.)

        do k=1,nplev
         khin=k
         khout=nplev+1-k
         if(orev)khout=k
         igout=ix/2+ix*(iy/2-1)+ix*iy*(khin-1)
         write(6,*)"************************************************k=",k
         write(6,*)"===> khin,datan(igout)=",khin,datan(igout)
         call sintp16(datan(1+ix*iy*(khin-1):ix*iy*khin),ix,iy,hgt(:,:,khout),glon,glat,sdiag,il)
         write(6,*)"khout,hgt(il/2,jl.2,khout)=",khout,hgt(il/2,jl/2,khout)
         write(6,*)'<=== model hgt(m) khin,khout=',khin,khout
         call findxn(hgt(:,:,khout),ifull,-1.e29,xa,kx,an,kn)
        enddo ! k

      endif ! ier = 0

      write(6,*)"==================================================u"

      ier = nf_inq_varid(ncid,'u',idvar)
      write(6,*)"ier=",ier," idvar=",idvar
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'zonal_wnd',idvar)
         write(6,*)"ier=",ier," idvar=",idvar
      endif

      call ncread_3d(ncid,iarch,idvar,ix,iy,nplev,datan)

      ! MJT quick fix
      ier = nf_get_att_real(ncid,idvar,'_FillValue',fill_float)
      if ( ier .ne. 0 ) then
        ier = nf_get_att_real(ncid,idvar,'missing_value',fill_float)  
      end if
      if ( ier .eq. 0 ) then
        where ( datan==fill_float )
          datan = 1.e10
        end where
      end if
      call getvalidlev(validlevhost,datan,ix,iy,nplev)
      call filldat(datan,ix,iy,nplev)

      call amap ( datan(1:ix*iy), ix, iy, 'input u', 0., 0. )

      do k=1,nplev
        write(6,*)"************************************************k=",k
        khin=k
        khout=nplev+1-k
        if(orev)khout=k
        call sintp16(datan(1+ix*iy*(khin-1):ix*iy*khin),ix,iy,u(:,:,khout),glon,glat,sdiag,il)
        call findxn(u(:,:,khout),ifull,-1.e29,xa,kx,an,kn)
      enddo

      write(6,*)"==================================================v"

      ier = nf_inq_varid(ncid,'v',idvar)
      write(6,*)"ier=",ier," idvar=",idvar
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'merid_wnd',idvar)
         write(6,*)"ier=",ier," idvar=",idvar
      endif

      call ncread_3d(ncid,iarch,idvar,ix,iy,nplev,datan)

      ! MJT quick fix
      ier = nf_get_att_real(ncid,idvar,'_FillValue',fill_float)
      if ( ier .ne. 0 ) then
        ier = nf_get_att_real(ncid,idvar,'missing_value',fill_float)  
      end if
      if ( ier .eq. 0 ) then
        where ( datan==fill_float )
          datan = 1.e10
        end where
      end if
      call getvalidlev(validlevhost,datan,ix,iy,nplev)
      call filldat(datan,ix,iy,nplev)

      call amap ( datan(1:ix*iy), ix, iy, 'input v', 0., 0. )

      do k=1,nplev
        write(6,*)"************************************************k=",k
        khin=k
        khout=nplev+1-k
        if(orev)khout=k
        call sintp16(datan(1+ix*iy*(khin-1):ix*iy*khin),ix,iy,v(:,:,khout),glon,glat,sdiag,il)
        call findxn(v(:,:,khout),ifull,-1.e29,xa,kx,an,kn)
      enddo

      write(6,*)"==================================================temp"

      ier = nf_inq_varid(ncid,'temp',idvar)
      write(6,*)"ier=",ier," idvar=",idvar
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'air_temp',idvar)
         write(6,*)"ier=",ier," idvar=",idvar
      endif

      call ncread_3d(ncid,iarch,idvar,ix,iy,nplev,datan)

      ! MJT quick fix 
      ier = nf_get_att_real(ncid,idvar,'_FillValue',fill_float)
      if ( ier .ne. 0 ) then
        ier = nf_get_att_real(ncid,idvar,'missing_value',fill_float)  
      end if
      if ( ier .eq. 0 ) then
        where ( datan==fill_float )
          datan = 1.e10
        end where
      end if
      call getvalidlev(validlevhost,datan,ix,iy,nplev)
      call filldat(datan,ix,iy,nplev)

      call amap ( datan(1:ix*iy), ix, iy, 'input temp', 0., 0. )

      do k=1,nplev
        write(6,*)"************************************************k=",k
        khin=k
        khout=nplev+1-k
        if(orev)khout=k
        call sintp16(datan(1+ix*iy*(khin-1):ix*iy*khin),ix,iy,temp(:,:,khout),glon,glat,sdiag,il)
        call findxn(temp(:,:,khout),ifull,-1.e29,xa,kx,an,kn)
      enddo

      write(6,*)"================================================rh/q"

      ier = nf_inq_varid(ncid,'rh',idvar)
      write(6,*)"ier=",ier," idvar=",idvar
      if ( ier .eq. 0 ) then
         moist_var="rh"
      endif
      
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'relhum',idvar)
         write(6,*)"ier=",ier," idvar=",idvar
         moist_var="rh"
      end if

      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'mix_rto',idvar)
         write(6,*)"ier=",ier," idvar=",idvar
         moist_var="mr"
      endif

      if ( ier .ne. 0 .and. moistvar .ne. '' ) then
         ier = nf_inq_varid(ncid,moistvar,idvar)
         write(6,*)"ier=",ier," idvar=",idvar
         moist_var="mr"
      endif

      write(6,*)"##################################moist_var=",moist_var

      call ncread_3d(ncid,iarch,idvar,ix,iy,nplev,datan)

      ! MJT quick fix 
      ier = nf_get_att_real(ncid,idvar,'_FillValue',fill_float)
      if ( ier .ne. 0 ) then
        ier = nf_get_att_real(ncid,idvar,'missing_value',fill_float)  
      end if
      if ( ier .eq. 0 ) then
        where ( datan==fill_float )
          datan = 1.e10
        end where
      end if
      call getvalidlev(validlevhost,datan,ix,iy,nplev)
      call filldat(datan,ix,iy,nplev)

      call findxn(datan,ix*iy*nplev,-1.e29,xa,kx,an,kn)

      if ( xa .lt. .1 ) then
        moist_var="mr"
        write(6,*)"################################moist_var=",moist_var
      endif ! ( xa .lt. 1.1 ) then

      call amap ( datan(1:ix*iy), ix, iy, 'input '//moist_var, 0., 0. )

      do k=1,nplev
        write(6,*)"************************************************k=",k
        khin=k
        khout=nplev+1-k
        if(orev)khout=k
        call sintp16(datan(1+ix*iy*(khin-1):ix*iy*khin),ix,iy,rh(:,:,khout),glon,glat,sdiag,il)

        write(6,*)"make sure data is always between 0 and 100!"
        write(6,*)"for both mixr and rh"
        rh(:,:,khout)=max(0.,min(100.,rh(:,:,khout)))
        call findxn(rh(:,:,khout),ifull,-1.e29,xa,kx,an,kn)

        if ( moist_var .eq. "rh" .and. xa .lt. 1.1 ) then

          write(6,*)"######################convert rh from 0-1 to 0-100"
          rh(:,:,khout)=max(0.,min(100.,rh(:,:,khout)*100.))
          call findxn(rh(:,:,khout),ifull,-1.e29,xa,kx,an,kn)

        endif ! ( moist_var .eq. "rh" .and. xa .lt. 1.1 ) then

      enddo

      call findxn(hgt(:,:, 1),ifull,-1.e29,xa,kx,an,kn)
      call findxn(hgt(:,:,nplev),ifull,-1.e29,xa,kx,an,kn)
      write(6,*)"nplev=",nplev

!############################################################################
! sfc data
!############################################################################

      write(6,*)"================================================mslp"

      ier = nf_inq_varid(ncid,'mslp',idvar)
      write(6,*)"ier=",ier," idvar=",idvar
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'pmsl',idvar)
         write(6,*)"ier=",ier," idvar=",idvar
      endif

      if ( ier .eq. 0 ) then
        call ncread_2d(ncid,iarch,idvar,ix,iy,datan(1:ix*iy))
        call amap ( datan(1:ix*iy), ix, iy, 'gbl mslp', 0., 0. )
        call sintp16(datan(1:ix*iy),ix,iy,pmsl,glon,glat,sdiag,il)
        write(6,*)" findxn model mslp(Pa)"
        call findxn(pmsl,ifull,-1.e29,xa,kx,an,kn)
        if ( an .gt. 2000. ) then
           write(6,*)"#########################convert pmsl to hPa"
           pmsl(1:ifull)=pmsl(1:ifull)/100. ! to convert to hPa
        endif ! ( an .gt. 2000. ) then
      else
        write(6,*)"No pmsl data found, setting to 0"
        pmsl(1:ifull) = 0.
      endif ! ier

      call prt_pan(pmsl,il,jl,2,'pmsl')

      write(6,*)"================================================zs"

      ier = nf_inq_varid(ncid,'zs',idvar) ! from input netcdf file
      write(6,*)"ier=",ier," idvar=",idvar
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'topo',idvar)
         write(6,*)"ier=",ier," idvar=",idvar
      endif
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'topog',idvar)
         write(6,*)"ier=",ier," idvar=",idvar
      endif
      
      if ( ier .eq. 0 ) then
         geopotunits="gpm"
         ier = nf_get_att_text(ncid,idvar,'units',geopotunits)
         call ncread_2d(ncid,iarch,idvar,ix,iy,datan(1:ix*iy))  ! zsi(m)
         call findxn(zsi_m,ifull,-1.e29,xa,kx,an,kn)
         if ( geopotunits=="m" ) then
            write(6,*) "Converting from meters to gpm"
            datan(1:ix*iy) = datan(1:ix*iy)*g
         end if
         call amap ( datan(1:ix*iy), ix, iy, 'gbl zs', 0., 0. )
         call sintp16(datan(1:ix*iy),ix,iy,zsi_m,glon,glat,sdiag,il)  ! (m)
         write(6,*)" findxn zsi_m(m)"
!        if ( an .gt. 2000. ) then
!           write(6,*)"#########################convert m2/s2 to m"
!           do i=1,ifull
!             zsi_m(i)=zsi_m(i)/9.80616
!           enddo ! i=1,ifull
!        endif ! ( an .gt. 2000. ) then
      else
           write(6,*)"No zs data found, setting to -999."
           zsi_m(1:ifull)=-999.
      endif ! ier

      call prt_pan(zs,il,jl,2,'zs(m)')
      call prt_pan(zsi_m,il,jl,2,'zsi_m(m)')

      write(6,*)"================================================land"
      write(6,*)"===================================== 1=land 0=ocean"

      ier = nf_inq_varid(ncid,'land',idvar)
      write(6,*)"ier=",ier," idvar=",idvar
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'sfc_lsm',idvar)
         write(6,*)"ier=",ier," idvar=",idvar
      endif
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'land_mask',idvar)
         write(6,*)"ier=",ier," idvar=",idvar
      endif
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'lnd_mask',idvar)
         write(6,*)"ier=",ier," idvar=",idvar
      endif

      if ( ier .eq. 0 ) then
         olsm_gbl = .true.
         call ncread_2d(ncid,1,idvar,ix,iy,datan(1:ix*iy))	! MJT quick fix 
         datan(1:ix*iy)=abs(datan(1:ix*iy))                 ! MJT quick fix
         lsm_gbl(1:ix*iy)=datan(1:ix*iy)
         ! MJT quick fix 
         where (abs(lsm_gbl(1:ix*iy))>=1.e10)
           lsm_gbl(1:ix*iy)=0.
         end where
         call amap ( lsm_gbl, ix, iy, 'lsm_gbl', 0., 0. )
         call sintp16(lsm_gbl,ix,iy,lsmg_m,glon,glat,sdiag,il)
         write(6,*)" findxn model lsmg_m"
         call findxn(lsmg_m,ifull,-1.e29,xa,kx,an,kn)
!        if ( xa .gt. 1.5 ) then
!           write(6,*)"#################convert  so that 0=ocean/1=land"
!           do i=1,ifull
!             lsmg_m(i)=lsmg_m(i)
!           enddo ! i=1,ifull
!        endif ! ( an .gt. 2000. ) then
      else
           write(6,*)"No landmask data found, setting to -999."
           lsmg_m(1:ifull)=-999.
      endif ! ier

      call prt_pan(lsmg_m,il,jl,2,'lsmg_m')

      write(6,*)"================================================ps"

      ier = nf_inq_varid(ncid,'ps',idvar)
      write(6,*)"ier=",ier," idvar=",idvar
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'sfc_pres',idvar)
         write(6,*)"ier=",ier," idvar=",idvar
      endif

      if ( ier .eq. 0 ) then
         call ncread_2d(ncid,iarch,idvar,ix,iy,datan(1:ix*iy))
         if (any(datan(1:ix*iy)>2000)) then
           datan(1:ix*iy)=datan(1:ix*iy)/100.
         end if
         call amap ( datan(1:ix*iy), ix, iy, 'gbl sfcp', 0., 0. )
         call sintp16(datan(1:ix*iy),ix,iy,psg_m,glon,glat,sdiag,il)
         ! remove levels below surface
         do j = 1,iy
           do i = 1,ix
             iq = i + ix*(j-1)
             do k = nint(validlevhost(iq)),nplev
               if ( datan(iq)>plev(k) ) then
                 validlevhost(iq) = real(k)
                 exit
               end if
             end do
           end do
         end do
         psg_m(1:ifull) = psg_m(1:ifull)*100. ! convert to Pa
      else
         write(6,*)"No sfcp data found, setting to -999."
         psg_m(1:ifull)=-999.
      end if ! ier

      call prt_pan(psg_m,il,jl,2,'psg_m')

      write(6,*)"================================================tss"

      ier = nf_inq_varid(ncid,'tss',idvar)
      write(6,*)"ier=",ier," idvar=",idvar
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'sfc_temp',idvar)
         write(6,*)"ier=",ier," idvar=",idvar
      endif

      if ( ier .eq. 0 ) then ! we have sfc temp data
        write(6,*)"input data has sfc temp data, now read in"
        call ncread_2d(ncid,iarch,idvar,ix,iy,datan(1:ix*iy))
        spval=-1.e10
        write(6,*)"spval=",spval
        if (any(datan(1:ix*iy)<100..or.datan(1:ix*iy)>400.)) then
          write(6,*) "Missing data found in sfc temp"
          where (datan(1:ix*iy)<100..or.datan(1:ix*iy)>400.)
            datan(1:ix*iy)=spval
          end where
          call fill(datan(1:ix*iy),ix,iy,.1*spval,datan(1+2*ix*iy:3*ix*iy))
        end if
        call amap ( datan(1:ix*iy), ix, iy, 'gbl sfct', 0., 0. )

        write(6,*)"###################### do we have olsm_gbl=",olsm_gbl
        ijgd=igd+ix*(jgd-1)
        write(6,*)"igd,jgd,ijgd=",igd,jgd,ijgd
        ijd=id+il*(jd-1)
        write(6,*)"id,jd,ijd=",id,jd,ijd

        write(6,*)"prepare to interp. tss for sea and land separately"
        write(6,*)"igd,jgd,gtss=",igd,jgd,datan(ijgd)
        write(6,*)"putting only land values into datan"
        write(6,*)"putting only ocean values into datan(+ix*iy)"
!not done since tss already at sea level     write(6,*)"First: reduce tss to sea level"

        nlpnts=0
        nopnts=0
        do j=1,iy
         do i=1,ix
          iq = i+(j-1)*ix
          datan(iq+ix*iy)=datan(iq)                          ! for ocean pts
          if(olsm_gbl)then
            if ( lsm_gbl(iq) .lt. .5 ) then
              datan(iq)=spval                               ! land, fill in ocean pts
              nlpnts=nlpnts+1
            else !!!  ( lsm_gbl(iq) .lt. .5 ) then
              datan(iq+ix*iy)=spval                          ! ocean, fill in land pts
              nopnts=nopnts+1
            endif ! ( lsm_gbl(iq) .gt. .5 ) then
          endif!(olsm_gbl)then
         enddo ! ix
        enddo ! iy

        write(6,*)"two global tss arrays with spval=", spval
        write(6,*)"igd,jgd,lgtss=",igd,jgd,datan(ijgd)
        write(6,*)"igd,jgd,ogtss=",igd,jgd,datan(ijgd+ix*iy)

        write(6,*)"fill in missing values nlpnts,nopnts=",nlpnts,nopnts

        write(6,*)"=======> for land array, fill in tss ocean values"
        call fill(datan(1:ix*iy),ix,iy,.1*spval,datan(1+2*ix*iy:3*ix*iy))

        if(olsm_gbl .and. nopnts.gt.0)then
           write(6,*)"=======> for ocean array, fill in tss land values"
           call fill(datan(1+ix*iy:2*ix*iy),ix,iy,.1*spval,datan(1+2*ix*iy:3*ix*iy))
        endif!(olsm_gbl)then

        write(6,*)"igd,jgd,lgtss=",igd,jgd,datan(ijgd)
        write(6,*)"igd,jgd,ogtss=",igd,jgd,datan(ijgd+ix*iy)


!not done since tss already at sea level     write(6,*)"tss at sea level here"

        write(6,*)"=========================> now interp. land data"
        call sintp16(datan(1:ix*iy),ix,iy,sfct,glon,glat,sdiag,il)                 ! land

        if(olsm_gbl .and. nopnts.gt.0)then
          write(6,*)"=========================> now interp. ocean data"
           call sintp16(datan(1+ix*iy:2*ix*iy),ix,iy,sfcto_m,glon,glat,sdiag,il)   ! ocean
        endif!(olsm_gbl)then

        call prt_pan(sfct   ,il,jl,2,'tss')
        call prt_pan(sfcto_m,il,jl,2,'tsso')

        write(6,*)"id,jd,ltss=",id,jd,sfct(ijd)
        write(6,*)"id,jd,otss=",id,jd,sfcto_m(ijd)

        if(olsm_gbl)then

          write(6,*)"now recombine two (land/ocean) fields"
!not done since tss already at sea level     write(6,*)"Also need to recompute tss at zs"

          do j=1,jl
           do i=1,il
            iq=i+(j-1)*il
!not done since tss already at sea level  sfct(i,j)=sfct(i,j)-zs(iq)*.0065
!not done since tss already at sea level  sfctl_m(i,j)=sfctl_m(i,j)-zs(iq)*.0065

! remeber, land < .5 is an ocean point
            if ( lsm_m(iq) .lt. .5 ) sfct(iq)=sfcto_m(iq)  ! set to ocean interp pnt

           enddo ! i
          enddo ! j

        endif!(olsm_gbl)then

        write(6,*)"id,jd,sfct=",id,jd,sfct(ijd)

        write(6,*)" findxn model sfct"
        call findxn(sfct,ifull,-1.e29,xa,kx,an,kn)

      else
        write(6,*)"###################no sfc temp data in input dataset"
        write(6,*)"###################setting sfc temp data to 0!!!!!!!"
        sfct(1:il*jl) = 0.
      endif ! ier eq 0 , sfct

      call prt_pan(sfct,il,jl,2,'sfct')

      write(6,*)"============================================fracice"

      fracice=-1.
      ier = nf_inq_varid(ncid,'fracice',idvar)
      write(6,*)"ier=",ier," idvar=",idvar
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'seaice',idvar)
         write(6,*)"ier=",ier," idvar=",idvar
      endif

      if ( ier .eq. 0 ) then ! we have fracice data
        write(6,*)"input data has fracice data, now read in"
        call ncread_2d(ncid,iarch,idvar,ix,iy,datan(1:ix*iy))
        where (datan(1:ix*iy).gt.1.01)
          datan(1:ix*iy)=0.
        end where

        call amap ( datan(1:ix*iy), ix, iy, 'gbl fice', 0., 0. )

        spval=-1.e10
        write(6,*)"spval=",spval
        write(6,*)"###################### do we have olsm_gbl=",olsm_gbl
        ijgd=igd+ix*(jgd-1)
        write(6,*)"igd,jgd,ijgd=",igd,jgd,ijgd
        ijd=id+il*(jd-1)
        write(6,*)"id,jd,ijd=",id,jd,ijd

        write(6,*)"prepare to interp. fracice"
        write(6,*)"igd,jgd,gtss=",igd,jgd,datan(ijgd)

        nlpnts=0
        nopnts=0
        do j=1,iy
         do i=1,ix
          iq = i+(j-1)*ix
          datan(iq+ix*iy)=datan(iq)                          ! for ocean pts
          if(olsm_gbl)then
            if ( lsm_gbl(iq) .ge. .5 ) then
              datan(iq+ix*iy)=spval                          ! ocean, fill in land pts
              nopnts=nopnts+1
            else
              nlpnts=nlpnts+1
            endif ! ( lsm_gbl(iq) .ge. .5 ) then
          endif!(olsm_gbl)then
         enddo ! ix
        enddo ! iy

        write(6,*)"global fracice array with spval=", spval
        write(6,*)"igd,jgd,ogtss=",igd,jgd,datan(ijgd+ix*iy)

        write(6,*)"fill in missing values nlpnts,nopnts=",nlpnts,nopnts

        write(6,*)"===> for ocean array, fill in fracice land values"
        call fill(datan(1+ix*iy:2*ix*iy),ix,iy,.1*spval,datan(1+2*ix*iy:3*ix*iy))

        write(6,*)"igd,jgd,ogtss=",igd,jgd,datan(ijgd+ix*iy)

        write(6,*)"=========================> now interp. ocean data"
        call sintp16(datan(1+ix*iy:2*ix*iy),ix,iy,fracice,glon,glat,sdiag,il)   ! ocean

        call prt_pan(fracice,il,jl,2,'fracice')

        where (lsm_m.ge.0.5)
          fracice=0.
        end where

        write(6,*)"id,jd,fracice=",id,jd,fracice(ijd)

        write(6,*)" findxn model fracice"
        call findxn(fracice,ifull,-1.e29,xa,kx,an,kn)

      else
        write(6,*)"###############no fracice data in input dataset"
        write(6,*)"###############setting fracice data to -1.!!!!!!!"
        fracice(1:ifull) = -1.
      endif ! ier eq 0 , sfct

      call prt_pan(fracice,il,jl,2,'fracice')

      write(6,*)"================================================snod"

      ier = nf_inq_varid(ncid,'snod',idvar)
      write(6,*)"ier=",ier," idvar=",idvar
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'snow_amt_lnd',idvar)
         write(6,*)"ier=",ier," idvar=",idvar
      endif

      snod(:) = 0.
      if ( ier .eq. 0 ) then
         call ncread_2d(ncid,iarch,idvar,ix,iy,datan(1:ix*iy))
         ! fix for missing values
         where ( datan(1:ix*iy) > 1000. )
           datan(1:ix*iy)=0.
         end where
         call amap ( datan(1:ix*iy), ix, iy, 'snod', 0., 0. )
         call sintp16(datan(1:ix*iy),ix,iy,snod,glon,glat,sdiag,il)
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
      ier = nf_inq_varid(ncid,'soil_temp',idvar)
      write(6,*)"ier=",ier," idvar=",idvar," nsoillvl=",nsoillvl

      if ( ier==0 ) then ! we have soil_temp data
        call ncread_3d(ncid,iarch,idvar,ix,iy,nsoillvl,datan)
        call filldat(datan,ix,iy,nsoillvl)
        call amap ( datan, ix, iy, 'input soil_temp', 0., 0. )        
        ! interpolate to CCAM soil levels
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
          call sintp16(datatemp(:),ix,iy,soiltemp(:,k),glon,glat,sdiag,il)
        end do
      end if
      
      write(6,*)"============================================soil_mois"

      soilmoist(:,:)=-1.
      ier = nf_inq_varid(ncid,'soil_mois',idvar)
      write(6,*)"ier=",ier," idvar=",idvar," nsoillvl=",nsoillvl

      if ( ier==0 ) then ! we have soil_moist data
        ! check units  
        ier = nf_get_att_text(ncid,idvar,'units',soilunits)
        write(6,*)"ier=",ier," soilunits=",trim(soilunits)        
          
        call ncread_3d(ncid,iarch,idvar,ix,iy,nsoillvl,datan)  
        call filldat(datan,ix,iy,nsoillvl)
        call amap ( datan, ix, iy, 'input soil_mois', 0., 0. )        
        ! interpolate to CCAM soil levels
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
          if ( trim(soilunits)=="kg m-2" ) then
            ! fix units from kg/m2 to m3/m3
            write(6,*) "Converting units from kg/m2 to m3/m3"
            datatemp(:) = datatemp(:)/soilthick_ccam(k)/rhowater  
          end if
          call sintp16(datatemp(:),ix,iy,soilmoist(:,k),glon,glat,sdiag,il)
        end do
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
       if ( ogbl ) then

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


!############### fix winds if DARLAM grid ###############################
       else ! not ogbl

! convert e-w/n-s lat/lon winds to model winds
! loop over all model grid points
        do k=1,nplev
         write(6,*)k,temp(1,1,k),u(1,1,k),v(1,1,k)
         do j=1,jl
          do i=1,il
! get lat lon of model grid ( just to get ther )
           call lconll(rlon,rlat,float(i),float(j))
! calculate ucmp l.c.winds
! ulc=v@u*s(th)+u@u*c(th)
           ull = u(i,j,k)
           vll = v(i,j,k)
           u(i,j,k)=vll*sin(ther)+ull*cos(ther)
! calculate vcmp l.c.winds
! vlc=v@v*c(th)-u@v*s(th)
           v(i,j,k)=vll*cos(ther)-ull*sin(ther)
          enddo ! i
         enddo ! j
         write(6,*)k,temp(1,1,k),u(1,1,k),v(1,1,k)
        enddo ! k

      endif ! not ogbl

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
          cplev(k)=tohpa*plev(nplev+1-k)
          bplev(k)=tohpa*plev_b(nplev+1-k)
          write(6,'(5f10.2,f10.5)') cplev(k),hgt(i,j,k),temp(i,j,k),u(i,j,k),v(i,j,k),rh(i,j,k)
      enddo ! k=1,nplev

      iq=il/2+(jl/2-1)*il
      write(6,'("pmsl=",f12.2," sfct=",f12.2)') pmsl(iq),sfct(iq)

      ! intepolate valid level
      call amap ( validlevhost, ix, iy, 'gbl levl', 0., 0. )
      call sintp16(validlevhost,ix,iy,validlevcc,glon,glat,sdiag,il)

      write(6,*)"calling vidar now!! ntimes,iarch=",ntimes,iarch

      if2=0

!#######################################################################
      call vidar(nplev,hgt,temp,u,v,rh,validlevcc,iyr,imn,idy,ihr,iarch,time,mtimer,cplev,bplev, &
                 io_out,il,kl,minlon,maxlon,minlat,maxlat,llrng,procformat_nproc)
!#######################################################################

      enddo ! narch

      deallocate(datan,zs_gbl,lsm_gbl)
      deallocate(glon,glat)
      deallocate(datatemp,validlevhost)

      write(6,*)'*********** Finished cdfvidar ************************'
      
      deallocate(x,y,z,ax,ay,az,bx,by,bz)
      deallocate(hgt,temp,u,v,rh)
      deallocate(sfcto_m,lsmg_m)
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
      subroutine ncread_2d(idhist,iarch,idvar,il,jl,var)

      use netcdf_m

      implicit none
      
!     include 'gblparm.h'

      integer idhist,il,jl,iarch,idvar
      integer i,j,ij,ier,itype
      integer start(3),count(3)

      real var(il*jl), addoff, sf
      real dx, dn

      integer*2, dimension(:), allocatable :: ivar
      double precision, dimension(:), allocatable :: dvar
!     integer*2 ivar(nnx*nny)

      character*30 name

      write(6,*)"=============================> ncread_2d idhist=",idhist
      write(6,*)"iarch=",iarch," idvar=",idvar
      write(6,*)"il=",il," jl=",jl

! read name
      ier = nf_inq_varname(idhist,idvar,name)
      write(6,*)"ier=",ier," name=",name

      if(ier.eq.0)then

!       if(il*jl.gt.nnx*nny)stop "ncread_2d il*jl.gt.nnx*nny"

        start(1) = 1
        start(2) = 1
        start(3) = iarch
        count(1) = il
        count(2) = jl
        count(3) = 1

        write(6,'("start=",4i4)') start
        write(6,'("count=",4i4)') count

      ier = nf_inq_vartype(idhist,idvar,itype)
      write(6,*)"itype=",itype," ier=",ier

      if ( itype .eq. nf_short ) then
         allocate(ivar(il*jl))
         write(6,*)"variable is short"
         ier = nf_get_vara_int2(idhist,idvar,start,count,ivar)
         write(6,*)"ivar(1)=",ivar(1)," ier=",ier
         write(6,*)"ivar(il*jl)=",ivar(il*jl)
      else if ( itype .eq. nf_float ) then
         write(6,*)"variable is float"
         ier = nf_get_vara_real(idhist,idvar,start,count,var)
         write(6,*)"var(1)=",var(1)," ier=",ier
         write(6,*)"var(il*jl)=",var(il*jl)
      else if ( itype .eq. nf_double ) then
         allocate(dvar(il*jl))
         write(6,*)"variable is double"
         ier = nf_get_vara_double(idhist,idvar,start,count,dvar)
         write(6,*)"dvar(1)=",dvar(1)," ier=",ier
         write(6,*)"dvar(il*jl)=",dvar(il*jl) 
      else
         write(6,*)"variable is unknown"
         call finishbanner
         stop -1
      endif

! obtain scaling factors and offsets from attributes
        ier = nf_get_att_real(idhist,idvar,'add_offset',addoff)
        if ( ier.ne.0 ) addoff=0.
        write(6,*)"ier=",ier," addoff=",addoff

        ier = nf_get_att_real(idhist,idvar,'scale_factor',sf)
        if ( ier.ne.0 ) sf=1.
        write(6,*)"ier=",ier," addoff=",addoff

      else!(ier.eq.0)then
! no data found
        do i=1,il*jl
         var(i)=0
        enddo
        sf=0.
        addoff=0.
      endif!(ier.eq.0)then

! unpack data
      dx=-1.e29
      dn= 1.e29
      do j=1,jl
        do i=1,il
          ij=i+(j-1)*il
      	  if ( itype .eq. nf_short ) then
           if(i.eq.1.and.j.eq.1) write(6,*)"ivar,sf,addoff=",ivar(ij),sf,addoff
            var(ij) = ivar(ij)*sf + addoff
          else if ( itype .eq. nf_float ) then
           if(i.eq.1.and.j.eq.1) write(6,*)"var,sf,addoff=",var(ij),sf,addoff
            var(ij) = var(ij)*sf + addoff
          else
           if(i.eq.1.and.j.eq.1) write(6,*)"var,sf,addoff=",dvar(ij),sf,addoff
            var(ij) = dvar(ij)*sf + addoff 
          endif
          dx=max(dx,var(ij))
          dn=min(dn,var(ij))
        end do
      end do

      write(6,*)"ncread_2d idvar=",idvar," iarch=",iarch
      write(6,*)"ncread_2d dx=",dx," dn=",dn

      if ( itype .eq. nf_short ) then
         deallocate(ivar)
      endif
      if ( itype .eq. nf_double ) then
         deallocate(dvar)
      end if

      return ! ncread_2d
      end
!***************************************************************************
      subroutine ncread_3d(idhist,iarch,idvar,il,jl,kl,var)
!             call ncread_3d(ncid,iarch,idvar,ix,iy,nplev,datan)

      use netcdf_m
      
      implicit none
      
!     include 'gblparm.h'

      integer, dimension(4) :: start,count
      integer, intent(in) :: idhist, iarch, idvar
      integer, intent(in) :: il, jl, kl
      integer ier, itype, i, j, k, ijk
      
      real addoff, sf
      real dx, dn

      integer*2, dimension(:), allocatable :: ivar
      double precision, dimension(:), allocatable :: dvar
      real, dimension(il*jl*kl), intent(out) :: var
      character*30 name

      write(6,*)"ncread_2d idhist=",idhist
      write(6,*)"iarch=",iarch," idvar=",idvar
      write(6,*)"il=",il," jl=",jl

      ier = nf_inq_varname(idhist,idvar,name)
      write(6,*)"ier=",ier," name=",name

      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = iarch

      count(1) = il
      count(2) = jl
      count(3) = kl
      count(4) = 1

      write(6,'("start=",4i4)') start
      write(6,'("count=",4i4)') count

! read data
      write(6,*)"idhist=",idhist," idvar=",idvar
      ier = nf_inq_vartype(idhist,idvar,itype)
      write(6,*)"ier=",ier," itype=",itype

      addoff=0.
      sf=1.
      
      if ( itype .eq. nf_short ) then
         write(6,*)"variable is short"
         allocate( ivar(il*jl*kl) )
         ier = nf_get_vara_int2(idhist,idvar,start,count,ivar)
      else if ( itype .eq. nf_float ) then
         write(6,*)"variable is float"
         ier = nf_get_vara_real(idhist,idvar,start,count,var)
      else if ( itype .eq. nf_double ) then
         write(6,*)"variable is double"
         allocate( dvar(il*jl*kl) )
         ier = nf_get_vara_double(idhist,idvar,start,count,dvar)
         var=real(dvar)
      else
         write(6,*)"variable is unknown"
         call finishbanner
         stop -1
      endif


      if ( itype .eq. nf_short ) then
      write(6,*)"obtain scaling factors and offsets from attributes"
      ier = nf_get_att_real(idhist,idvar,'add_offset',addoff)
      if ( ier.ne.0 ) addoff=0.
      write(6,*)"ier=",ier," addoff=",addoff

      ier = nf_get_att_real(idhist,idvar,'scale_factor',sf)
      if ( ier.ne.0 ) sf=1.
      write(6,*)"ier=",ier," sf=",sf
      endif

! unpack data
      dx=-1.e29
      dn= 1.e29
      do k=1,kl
       do j=1,jl
        do i=1,il
          ijk=i+(j-1)*il+(k-1)*il*jl
          if(i.eq.1.and.j.eq.1.and.k.eq.1) write(6,*)"i,j,k,ijk=",i,j,k,ijk
      	  if ( itype .eq. nf_short ) then
            if(i.eq.1.and.j.eq.1.and.k.eq.1) write(6,*)"ivar,sf,addoff=",ivar(ijk),sf,addoff
            var(ijk) = real(ivar(ijk))*sf + addoff
          else
            if(i.eq.1.and.j.eq.1.and.k.eq.1) write(6,*)"var,sf,addoff=",var(ijk),sf,addoff
            var(ijk) = var(ijk)*sf + addoff
          endif
          if(i.eq.1.and.j.eq.1.and.k.eq.1) write(6,*)"var=",var(ijk)
          dx=max(dx,var(ijk))
          dn=min(dn,var(ijk))
        end do
       end do
      end do
      
      if ( itype .eq. nf_short ) then
        deallocate( ivar )    
      end if
      if ( itype .eq. nf_double ) then
        deallocate( dvar )  
      end if

      write(6,*)"ncread_3d idvar=",idvar," iarch=",iarch
      write(6,*)"ncread_3d dx=",dx," dn=",dn

      return ! ncread_3d
      end
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

      subroutine getvalidlev(validlev,datan,ix,iy,plev)
      
      implicit none
      
      integer, intent(in) :: ix,iy,plev
      integer i,j,k,iq,iqk
      real, dimension(ix*iy), intent(inout) :: validlev
      real, dimension(ix*iy*plev), intent(in) :: datan
      
      do j = 1,iy
        do i = 1,ix
          iq = i + ix*(j-1)
          do k = nint(validlev(iq)),plev	
            iqk = i + ix*(j-1) + ix*iy*(k-1)
            if ( datan(iqk)<1.E10 ) then
              validlev(iq) = real(k)
              exit
            end if
          end do
        end do
      end do
      
      return
      end subroutine getvalidlev
      
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

