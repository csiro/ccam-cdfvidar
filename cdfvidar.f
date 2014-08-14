      program cdfvidar

      use ccinterp
      use cll_m
      use comsig_m      
      use latlong_m
      use sigdata_m
      
      implicit none
      
      integer iout,in
      integer kvin,kuin,kbin
      integer ism,nrh
      integer krin,ktin,khin
      integer mtimer,id,jd
      integer mxcyc,nvsig,ntimes
      integer insm,klmax,l
      integer ilx,jlx,i,j,iq
      integer nplev,ilonx,ilatx
      integer ier,ijd,ilt,jlk
      integer irecd,ngatts,nvars
      integer ndims,idv,lonid,latid
      integer k,ivpres
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
      real elon,elat,xplev
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

      character*80 inf

      common/mapproj/du,tanl,rnml,stl1,stl2
      include 'nplevs.h' ! maxplev
      real, dimension(maxplev) :: plev,cplev
      common/levpre/nplev,plev

      include 'netcdf.inc'

      parameter ( pi=3.1415926536 )
      parameter ( g=9.80616 )
      parameter ( klmax = 100 )
      
      real, dimension(klmax) :: dsg,sgml

      common/ncdfids/dimil,dimjl,dimkl,dimtim
     &              ,idil,idjl,idkl,idnt
* dimension ids
      integer  dimil,dimjl,dimkl,dimtim
* variable ids
      integer  idil,idjl,idkl,idnt,ix,iy
      integer  il,jl,kl,ifull
      integer ncid
      integer iernc,liernc,lncid,varid,dimid
      integer, dimension(2) :: ccdim

      common/lconther/ther
      include 'vidar.h'
      logical sdiag

      real dst
      real, dimension(2) :: lonlat
      real, dimension(:,:,:), allocatable :: rlld,xyz,axyz,bxyz      
      real, dimension(:,:), allocatable :: grid
      real, dimension(:), allocatable :: ax,ay,az,bx,by,bz,x,y,z
      real, dimension(:), allocatable :: datan
      real, dimension(:), allocatable :: zs_gbl, lsm_gbl
      real, dimension(:), allocatable :: glon, glat

      real plevin    (maxplev)
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

      common/datatype/moist_var,in_type
      character*1 in_type
      character*2 moist_var

      logical ofirst, ogbl, orev, olsm_gbl
      character*60 timorg
      character*60 cu
      character*3 cmonth
      character*80 zsavn,lsavn
      character*10 header,moistvar

      namelist/gnml/inf,vfil,ds,du,tanl,rnml,stl1,stl2,inzs,zsfil
     &             ,ints,tsfil, ogbl,zsavn,inzsavn,lsavn,inlsavn
     &             ,plevin,orev,io_out,igd,jgd,id,jd,mtimer,ntimes
     &             ,spline,mxcyc,nvsig,nrh
     &             ,oesig,sgml,dsg,ptop,debug,notop,opre,have_gp
     &             ,in,calout
     &             ,iout,oform,sdiag
     &             ,insm,smfil
     &             ,splineu,splinev,splinet,zerowinds
     &             ,grdx,grdy,slon,slat
     &             ,moistvar,kl

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
      data sgml/klmax*0./, dsg/klmax*0./
      data plevin/maxplev*0./
      data splineu/.true./, splinev/.true./, splinet/.true./
      data sdiag/.false./
      data have_gp/.true./
      data zerowinds/.true./
      data grdx/1./
      data grdy/-1./

      save

   
      slon=0.
      slat=90.

!####################### read namelists ############################
      write(6,*)'read namelist'
      read (5, gnml)
      write(6,nml=gnml)
! read and write namelist input for vidar
!     open  ( 98, file='vidar.nml',status='unknown' )
!     read  ( 98, nml=vi )
!     write ( unit=6, nml=vi)
!####################### read namelist ############################

      if (kl.gt.klmax) then
        write(6,*) "ERROR: kl is greater than klmax"
	stop
      end if

      call comsigalloc(kl)
      dsgx=dsg(1:kl)
      sgmlx=sgml(1:kl)

      spline = splineu .or. splinev .or. splinet

! set up what sigma levels the outgoing data will have
! assumes top down, ie. sg(1)=0., dsg>0

           sgx(1)=0.
           if ( sgmlx(kl/2).gt.0. ) then
c dsg=0, sgml>0
              do l=2,kl
                sgx(l)=.5*(sgmlx(l-1)+sgmlx(l))
              end do ! l=2,kl
              do l=2,kl+1
                dsgx(l-1)=sgx(l)-sgx(l-1)
              end do ! l=2,kl+1
           elseif ( dsgx(kl/2).gt.0. ) then
c sgml=0, dsg>0
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
              stop
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
            stop
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

c       convert conformal cubic lats & longs to degrees (-90 to 90) & (0 to 360)
c       used in sint16; N.B. original rlong is -pi to pi
          !rlat(iq)=rlat(iq)*180./pi
          !rlong(iq)=rlong(iq)*180./pi
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
        read(inzs,*,err=25)ilt,jlk,ds,du,tanl,rnml,stl1,stl2
 25     if(ilt.eq.0.or.jlk.eq.0)then
           write(6,*)'no header in newtopo file'
        else
           write(6,*)'Header information for topofile'
           write(6,*)'ilt,jlk,ds,du,tanl,rnml,stl1,stl2'
     &           ,ilt,jlk,ds,du,tanl,rnml,stl1,stl2
           if(ilt.ne.il.or.jlk.ne.jl)stop 'wrong topofile supplied'
        endif     ! (ilt.eq.0.or.jlk.eq.0)
        write(6,*)"set up model grid params by calling lconset ds=",ds
        call lconset(ds)
      endif ! ( ogbl ) then

      if (liernc==nf_noerr) then
        write(6,*)'read model grid zsg = g*zs'
        iernc=nf_inq_varid(lncid,"zs",varid)
        iernc=nf_get_var_real(lncid,varid,zsg)
        write(6,*)'read model grid land-sea mask (0=ocean, 1=land)'
        iernc=nf_inq_varid(lncid,"lsm",varid)
        iernc=nf_get_var_real(lncid,varid,lsm_m)
        iernc=nf_close(lncid)
      else
        write(6,*)'read model grid zsg = g*zs'
        read(inzs,*)zsg

        write(6,*)'convert g*zs to zs(m)'
        do iq=1,ifull
          zs(iq)=zsg(iq)/g ! convert ascii read in zs*g to zs(m)
        enddo !iq=1,ifull

        write(6,*)'read model grid land-sea mask (0=ocean, 1=land)'
        read(inzs,*)lsm_m
        close(inzs)
      end if

      ijd=id+il*(jd-1)
      write(6,*)"ijd=",ijd," zs(m)=",zs(ijd)," lsm_m=",lsm_m(ijd)
!####################### read topography data ############################

!####################### open input netcdf file ############################
      write(6,*)'inf='
      write(6,*)inf
      ncid = ncopn(inf,ncnowrit,ier)
      write(6,*)'ncid=',ncid
      if(ier.ne.0) then
        write(6,*)' cannot open netCDF file; error code ',ier
        stop
      end if

!####################### get attributes of input netcdf file ############################
      call ncinq(ncid,ndims,nvars,ngatts,irecd,ier)
      write(6,'("ndims,nvars,ngatts,irecd,ier")')
      write(6,'(5i6)') ndims,nvars,ngatts,irecd,ier

c Get dimensions
      write(6,*) "get dim1 ncid=",ncid
c turn OFF fatal netcdf errors
      call ncpopt(0)
      lonid = ncdid(ncid,'lon',ier)
      write(6,*)"lon ncid,lonid,ier=",ncid,lonid,ier
c turn on fatal netcdf errors
c     write(6,*)"NCVERBOS,NCFATAL=",NCVERBOS,NCFATAL
c     call ncpopt(NCVERBOS+NCFATAL)
      if ( ier.eq.0 ) then
        write(6,*)"ncid,lonid=",ncid,lonid
        ier= nf_inq_dimlen(ncid,lonid,ix)
        write(6,*)"input ix,ier=",ix,ier
        latid= ncdid(ncid,'lat',ier)
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
        lonid = ncdid(ncid,'longitude',ier)
        write(6,*)"lonid=",lonid," ier=",ier
        ier= nf_inq_dimlen(ncid,lonid,ix)
        write(6,*)"input ix=",ix," ier=",ier

        latid= ncdid(ncid,'latitude',ier)
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
      write(6,*)"==================> slon=",slon," elon=",elon," ix=",ix
      write(6,*)"==================> slat=",slat," elat=",elat," iy=",iy

      ier = nf_inq_dimid(ncid,'pres',idpres)
      in_type="p"
      if ( ier .ne. 0 ) then
         ier = nf_inq_dimid(ncid,'plev',idpres)
         in_type="p"
      endif
      if ( ier .ne. 0 ) then
         ier = nf_inq_dimid(ncid,'lvl',idpres)
         in_type="s"
      endif
      write(6,*)"ier=",ier," idpres=",idpres," in_type=",in_type
      
      ier= nf_inq_dimlen(ncid,idpres,nplev)
      write(6,*)"ier=",ier," nplev=",nplev

      ier = nf_inq_varid(ncid,'pres',ivpres)
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'plev',ivpres)
      endif
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'lvl',ivpres)
      endif
      write(6,*)"ier=",ier," ivpres=",ivpres

      ier = nf_get_var_real(ncid,ivpres,plev)
      write(6,*)"ier=",ier," ivpres=",ivpres
      write(6,*)"input nplev=",nplev
      write(6,*)"plevs=",(plev(k),k=1,nplev)

      orev = plev(nplev).gt.plev(1)
      write(6,*)"#################################### orev=",orev

      write(6,*)"allocate arrays ix,iy,nplevs=",ix,iy,nplev
      allocate(datan(ix*iy*nplev))
      allocate(zs_gbl(ix*iy))
      allocate(lsm_gbl(ix*iy))
      allocate(validlevhost(ix*iy))

      if(orev) then
        do k=1,nplev
          datan(k)=plev(k)
        enddo
        do k=1,nplev
          plev(k)=datan(nplev+1-k)
        enddo
      endif

      xplev = -1.
      do k=1,nplev
         xplev=max(xplev,plev(k))
      enddo
      write(6,*)"xplev=",xplev

      osig_in = .false.
      if ( .01 .lt. xplev .and. xplev .lt. 800.  ) then
        write(6,*)"^^^^^^^^^actualy sigma levels^^^^^^^ fix plevs"
        osig_in = .true.
        do k=1,nplev
          plev(k)=plev(k)*1000. !
        enddo
        write(6,*)"plevs=",(plev(k),k=1,nplev)
      else if ( xplev .le. .01  ) then
        write(6,*)"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ fix plevs"
        do k=1,nplev
          plev(k)=plevin(k)
        enddo
        write(6,*)"plevs=",(plev(k),k=1,nplev)
        stop 'xplev < 800 in cdfvidar'
      endif

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

      if (ier.eq.0) then
        i=index(timorg,'since')
      else
        timorg='hours'
        i=0
      end if

      if (i.ne.0) then
        i=scan(timorg,' ')-1
        cu=''
        cu(1:i)=timorg(1:i)
        timorg(1:19)=timorg(i+8:i+26)
        read(timorg(1:4),*) iyr
        read(timorg(6:7),*) imn
        read(timorg(9:10),*) idy
        read(timorg(12:13),*) ihr
        read(timorg(15:16),*) imi 
      else
        cu=timorg
        ier = nf_get_att_text(ncid,ivtim,'time_origin',timorg)
        write(6,*)"ier=",ier," timorg=",timorg
        if (ier.ne.0) stop "timorg"
        read(timorg,'(i2)') idy
        read(timorg,'(3x,a3)') cmonth
        write(6,*)"cmonth=",cmonth
        imn = icmonth_to_imn(cmonth)
        write(6,*)"imn=",imn
        read(timorg,'(9x,i2)') iyr
        read(timorg,'(12x,i2)') ihr
        read(timorg,'(15x,i2)') imi
      end if

      ! disabled by MJT
      !if ( iyr .lt. 10 ) iyr = iyr+2000
      !if ( iyr .lt. 100 ) iyr = iyr+1900

      write(6,'("iyr,imn,idy,ihr,imi=",5i4)')iyr,imn,idy,ihr,imi

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

      call sintp16(datan(1+ix*iy:2*ix*iy),ix,iy,clat,glon,glat,
     &             sdiag,il)

      call prt_pan(clat,il,jl,2,'clat pan2')

      write(6,'("ix,iy,nplev,narch=",4i5)')ix,iy,nplev,narch

      write(6,*)"++++++++++++++++++++++++++++++++++++++++++++++++++++++"

      write(6,*)" nplev=",nplev
      write(6,*)" inlsavn=",inlsavn

      write(6,*)'read land-sea mask (0=ocean, 1=land) inlsavn=',inlsavn
      if(inlsavn.gt.0)then
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
      if(inzsavn.gt.0)then
        write(6,*)"open sfc. orog. for avn data"
        open(inzsavn,file=zsavn,form='formatted',status='old')
! read in free formatted avn sfc. orography (28-mar-2000)
! note that it runs north to south
!       rewind inzsavn
        read(inzsavn,*)((zs_gbl(i+(j-1)*ix),i=1,ix),j=1,iy)
        call amap ( zs_gbl, ix, iy, 'gbl sfczs(m)', 0., 0. )
        !write(6,*)"interp. zsavn to output grid"
        !call sintp16(zs_gbl,ix,iy,zs,glon,glat,sdiag)
        !write(6,*) 'findxn model sfc.height (m)'
        !call findxn(zs,ifull,-1.e29,xa,kx,an,kn)
        write(6,*) 'close unit inzsavn=',inzsavn
        close(inzsavn)
      endif

      write(6,*)' reading variables '

c***********************************************************************
      do iarch=1,narch
c***********************************************************************

      sdiag=.false.

      ier = nf_inq_varid(ncid,'time',ivtim)
      ier = nf_get_var1_real(ncid,ivtim,iarch,time)
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
	    stop
      end select

      write(6,*)"time=",time
      
      if ( time>1.E9 ) then
        write(6,*) "ERROR: Time too large for real variable"
        write(6,*) "Consider adjusting base date"
        stop
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

      write(6,*)ncid,iarch,idvar,ix,iy,nplev
      call ncread_3d(ncid,iarch,idvar,ix,iy,nplev,datan)

      call amap (datan(1:ix*iy),ix,iy,'input hgt',0.,0.)
      call amap (datan(1+ix*iy*(nplev-1):ix*iy*nplev),ix,iy,
     &           'input hgt',0.,0.)

      do k=1,nplev
       khin=k
       khout=nplev+1-k
       if(orev)khout=k
       igout=ix/2+ix*(iy/2-1)+ix*iy*(khin-1)
       write(6,*)"************************************************k=",k
       write(6,*)"===> khin,datan(igout)=",khin,datan(igout)
       call sintp16(datan(1+ix*iy*(khin-1):ix*iy*khin),ix,iy,
     &              hgt(:,:,khout),glon,glat,sdiag,il)
       write(6,*)"khout,hgt(il/2,jl.2,khout)="
     &           ,khout,hgt(il/2,jl/2,khout)
       write(6,*)'<=== model hgt(m) khin,khout=',khin,khout
       call findxn(hgt(:,:,khout),ifull,-1.e29,xa,kx,an,kn)
      enddo ! k

      !call prt_pan(hgt(1,1, 1),il,jl,1,'hgt: 1')
      !call prt_pan(hgt(1,1, 1),il,jl,2,'hgt: 1')
      !call prt_pan(hgt(1,1, 1),il,jl,3,'hgt: 1')
      !call prt_pan(hgt(1,1, 1),il,jl,4,'hgt: 1')
      !call prt_pan(hgt(1,1, 1),il,jl,5,'hgt: 1')
      !call prt_pan(hgt(1,1, 1),il,jl,6,'hgt: 1')
      !call prt_pan(hgt(1,1,nplev),il,jl,1,'hgt: nplev')
      !call prt_pan(hgt(1,1,nplev),il,jl,2,'hgt: nplev')
      !call prt_pan(hgt(1,1,nplev),il,jl,3,'hgt: nplev')
      !call prt_pan(hgt(1,1,nplev),il,jl,4,'hgt: nplev')
      !call prt_pan(hgt(1,1,nplev),il,jl,5,'hgt: nplev')
      !call prt_pan(hgt(1,1,nplev),il,jl,6,'hgt: nplev')

      !if ( k.gt.0 ) stop

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
      validlevhost=1.
      call getvalidlev(validlevhost,datan,ix,iy,nplev)
      call filldat(datan,ix,iy,nplev)

      call amap ( datan(1:ix*iy), ix, iy, 'input u', 0., 0. )

      do k=1,nplev
        write(6,*)"************************************************k=",k
        khin=k
        khout=nplev+1-k
        if(orev)khout=k
c       igout=ix/2+ix*(iy/2-1)+ix*iy*(khin-1)
c       write(6,*)khin,datan(igout)

        call sintp16(datan(1+ix*iy*(khin-1):ix*iy*khin),ix,iy,
     &      u(:,:,khout),glon,glat,sdiag,il)

c       write(6,*)khout,u(il/2,jl/2,khout)
c       write(6,*)'model u(m) khin,khout=',khin,khout
        call findxn(u(:,:,khout),ifull,-1.e29,xa,kx,an,kn)
      enddo

      !call prt_pan(u(1,1, 1),il,jl,2,'u : 1')
      !call prt_pan(u(1,1, 1),il,jl,1,'u : 1')
      !call prt_pan(u(1,1,nplev),il,jl,2,'u : nplev')
      !call prt_pan(u(1,1,nplev),il,jl,1,'u : nplev')

      write(6,*)"==================================================v"

      ier = nf_inq_varid(ncid,'v',idvar)
      write(6,*)"ier=",ier," idvar=",idvar
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'merid_wnd',idvar)
         write(6,*)"ier=",ier," idvar=",idvar
      endif

      call ncread_3d(ncid,iarch,idvar,ix,iy,nplev,datan)

      ! MJT quick fix 
      call getvalidlev(validlevhost,datan,ix,iy,nplev)
      call filldat(datan,ix,iy,nplev)

      call amap ( datan(1:ix*iy), ix, iy, 'input v', 0., 0. )

      do k=1,nplev
        write(6,*)"************************************************k=",k
        khin=k
        khout=nplev+1-k
        if(orev)khout=k
c       igout=ix/2+ix*(iy/2-1)+ix*iy*(khin-1)
c       write(6,*)khin,datan(igout)

        call sintp16(datan(1+ix*iy*(khin-1):ix*iy*khin),ix,iy,
     &      v(:,:,khout),glon,glat,sdiag,il)

c       write(6,*)khout,v(il/2,jl/2,khout)
c       write(6,*)'model v(m) khin,khout=',khin,khout
        call findxn(v(:,:,khout),ifull,-1.e29,xa,kx,an,kn)
      enddo

      !call prt_pan(v(1,1, 1),il,jl,2,'v : 1')
      !call prt_pan(v(1,1, 1),il,jl,1,'v : 1')
      !call prt_pan(v(1,1,nplev),il,jl,2,'v : nplev')
      !call prt_pan(v(1,1,nplev),il,jl,1,'v : nplev')

      write(6,*)"==================================================temp"

      ier = nf_inq_varid(ncid,'temp',idvar)
      write(6,*)"ier=",ier," idvar=",idvar
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'air_temp',idvar)
         write(6,*)"ier=",ier," idvar=",idvar
      endif

      call ncread_3d(ncid,iarch,idvar,ix,iy,nplev,datan)

      ! MJT quick fix 
      call getvalidlev(validlevhost,datan,ix,iy,nplev)
      call filldat(datan,ix,iy,nplev)


      call amap ( datan(1:ix*iy), ix, iy, 'input temp', 0., 0. )

      do k=1,nplev
        write(6,*)"************************************************k=",k
        khin=k
        khout=nplev+1-k
        if(orev)khout=k
c       igout=ix/2+ix*(iy/2-1)+ix*iy*(khin-1)
c       write(6,*)khin,datan(igout)
       call sintp16(datan(1+ix*iy*(khin-1):ix*iy*khin),ix,iy,
     &      temp(:,:,khout),glon,glat,sdiag,il)
c       write(6,*)khout,temp(il/2,jl/2,khout)
c       write(6,*)'model temp(m) khin,khout=',khin,khout
        call findxn(temp(:,:,khout),ifull,-1.e29,xa,kx,an,kn)
      enddo

      !call prt_pan(temp(1,1, 1),il,jl,1,'temp:  1')
      !call prt_pan(temp(1,1, 1),il,jl,2,'temp:  1')
      !call prt_pan(temp(1,1, 1),il,jl,3,'temp:  1')
      !call prt_pan(temp(1,1, 1),il,jl,4,'temp:  1')
      !call prt_pan(temp(1,1, 1),il,jl,5,'temp:  1')
      !call prt_pan(temp(1,1, 1),il,jl,6,'temp:  1')
      !call prt_pan(temp(1,1,nplev),il,jl,1,'temp: nplev')
      !call prt_pan(temp(1,1,nplev),il,jl,2,'temp: nplev')
      !call prt_pan(temp(1,1,nplev),il,jl,3,'temp: nplev')
      !call prt_pan(temp(1,1,nplev),il,jl,4,'temp: nplev')
      !call prt_pan(temp(1,1,nplev),il,jl,5,'temp: nplev')
      !call prt_pan(temp(1,1,nplev),il,jl,6,'temp: nplev')

      write(6,*)"================================================rh/q"

      ier = nf_inq_varid(ncid,'rh',idvar)
      write(6,*)"ier=",ier," idvar=",idvar
      if ( ier .eq. 0 ) then
         moist_var="rh"
      endif

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

c       igout=ix/2+ix*(iy/2-1)+ix*iy*(khin-1)
c       write(6,*)khin,datan(igout)

        call sintp16(datan(1+ix*iy*(khin-1):ix*iy*khin),ix,iy,
     &      rh(:,:,khout),glon,glat,sdiag,il)

        write(6,*)"make sure data is always between 0 and 100!"
        write(6,*)"for both mixr and rh"
        !do i=1,ifull
           rh(:,:,khout)=max(0.,min(100.,rh(:,:,khout)))
        !enddo !i=1,ifull

c       write(6,*)khout,rh(il/2,jl/2,khout)
c       write(6,*)'model rh(m) khin,khout=',khin,khout

        call findxn(rh(:,:,khout),ifull,-1.e29,xa,kx,an,kn)

        if ( moist_var .eq. "rh" .and. xa .lt. 1.1 ) then

          write(6,*)"######################convert rh from 0-1 to 0-100"
          rh(:,:,khout)=max(0.,min(100.,rh(:,:,khout)*100.))
          call findxn(rh(:,:,khout),ifull,-1.e29,xa,kx,an,kn)

        endif ! ( moist_var .eq. "rh" .and. xa .lt. 1.1 ) then

      enddo


!############################################################################
! sfc data
!############################################################################
      call findxn(hgt(:,:, 1),ifull,-1.e29,xa,kx,an,kn)
      call findxn(hgt(:,:,nplev),ifull,-1.e29,xa,kx,an,kn)
      write(6,*)"nplev=",nplev

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
           do i=1,ifull
             pmsl(i)=pmsl(i)/100. ! to convert to hPa
           enddo ! i=1,ifull
        endif ! ( an .gt. 2000. ) then

      else
           write(6,*)"No pmsl data found, setting to 0"
           do i=1,ifull
             pmsl(i)=0.
           enddo ! i=1,ifull
      endif ! ier

      call prt_pan(pmsl,il,jl,2,'pmsl')
      !call prt_pan(pmsl,il,jl,1,'pmsl')

      write(6,*)"================================================zs"

      ier = nf_inq_varid(ncid,'zs',idvar) ! from input netcdf file
      write(6,*)"ier=",ier," idvar=",idvar
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'topo',idvar)
         write(6,*)"ier=",ier," idvar=",idvar
      endif

      if ( ier .eq. 0 ) then

         call ncread_2d(ncid,iarch,idvar,ix,iy,datan(1:ix*iy))  ! zsi(m)

         call amap ( datan(1:ix*iy), ix, iy, 'gbl zs', 0., 0. )

         call sintp16(datan(1:ix*iy),ix,iy,zsi_m,glon,glat,sdiag,il)  ! (m)

         write(6,*)" findxn zsi_m(m)"
         call findxn(zsi_m,ifull,-1.e29,xa,kx,an,kn)

!        if ( an .gt. 2000. ) then
!           write(6,*)"#########################convert m2/s2 to m"
!           do i=1,ifull
!             zsi_m(i)=zsi_m(i)/9.80616
!           enddo ! i=1,ifull
!        endif ! ( an .gt. 2000. ) then

      else
           write(6,*)"No zs data found, setting to -999."
           do i=1,ifull
             zsi_m(i)=-999.
           enddo ! i=1,ifull
      endif ! ier

      call prt_pan(zs,il,jl,2,'zs(m)')
      call prt_pan(zsi_m,il,jl,2,'zsi_m(m)')
      !call prt_pan(zsi_m,il,jl,1,'zsi_m')

      write(6,*)"================================================land"
      write(6,*)"===================================== 1=land 0=ocean"

      ier = nf_inq_varid(ncid,'land',idvar)
      write(6,*)"ier=",ier," idvar=",idvar
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'sfc_lsm',idvar)
         write(6,*)"ier=",ier," idvar=",idvar
      endif

      if ( ier .eq. 0 ) then

         olsm_gbl = .true.

!         call ncread_2d(ncid,iarch,idvar,ix,iy,lsm_gbl)
         call ncread_2d(ncid,1,idvar,ix,iy,datan(1:ix*iy))	! MJT quick fix 
	 datan(1:ix*iy)=abs(datan(1:ix*iy))                     ! MJT quick fix
         lsm_gbl=datan(1:ix*iy)

         ! MJT quick fix 
         where (abs(lsm_gbl(:)).ge.1.e10)
           lsm_gbl(:)=0.
         end where

         call amap ( lsm_gbl, ix, iy, 'lsm_gbl', 0., 0. )

         call sintp16(lsm_gbl,ix,iy,lsmg_m,glon,glat,sdiag,il)
         
         !lsmg_m=real(nint(lsmg_m))

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
           do i=1,ifull
             lsmg_m(i)=-999.
           enddo ! i=1,ifull
      endif ! ier

      call prt_pan(lsmg_m,il,jl,2,'lsmg_m')
      !call prt_pan(lsmg_m,il,jl,1,'lsmg_m')

      write(6,*)"================================================ps"

      ier = nf_inq_varid(ncid,'ps',idvar)
      write(6,*)"ier=",ier," idvar=",idvar
      if ( ier .ne. 0 ) then
         ier = nf_inq_varid(ncid,'sfc_pres',idvar)
         write(6,*)"ier=",ier," idvar=",idvar
      endif

      if ( ier .eq. 0 ) then

         call ncread_2d(ncid,iarch,idvar,ix,iy,datan(1:ix*iy))

         if (any(datan(1:ix*iy).gt.2000)) then
           datan(1:ix*iy)=datan(1:ix*iy)/100.
         end if

         call amap ( datan(1:ix*iy), ix, iy, 'gbl sfcp', 0., 0. )

         call sintp16(datan(1:ix*iy),ix,iy,psg_m,glon,glat,sdiag,il)

         ! remove levels below surface
         do j=1,iy
           do i=1,ix
             iq=i+ix*(j-1)
             do k=nint(validlevhost(iq)),nplev
               if (datan(iq).gt.plev(k)) then
                 validlevhost(iq)=real(k)
                 exit
               end if
             end do
           end do
         end do

      else
           write(6,*)"No sfcp data found, setting to -999."
           do i=1,ifull
             psg_m(i)=-999.
           enddo ! i=1,ifull
      endif ! ier

      call prt_pan(psg_m,il,jl,2,'psg_m')
      !call prt_pan(psg_m,il,jl,1,'psg_m')

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
	
	if (any(datan(1:ix*iy).gt.400.)) then
	  write(6,*) "Missing data found in sfc temp"
	  where (datan(1:ix*iy).gt.400.)
	    datan(1:ix*iy)=spval
	  end where
        call fill(datan(1:ix*iy),ix,iy,.1*spval,
     &            datan(1+2*ix*iy:3*ix*iy))
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

!         write(6,*)i,j,iq,datan(iq) ,lsm_gbl(iq)
!not done since tss already at sea level       datan(iq)=datan(iq)+zs_gbl(iq)*.0065

          datan(iq+ix*iy)=datan(iq)                          ! for ocean pts

          if(olsm_gbl)then

            !if ( lsm_gbl(iq) .gt. .5 ) then
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
        call fill(datan(1:ix*iy),ix,iy,.1*spval,
     &            datan(1+2*ix*iy:3*ix*iy))

        if(olsm_gbl .and. nopnts.gt.0)then
           write(6,*)"=======> for ocean array, fill in tss land values"
           call fill(datan(1+ix*iy:2*ix*iy),ix,iy,.1*spval,
     &               datan(1+2*ix*iy:3*ix*iy))
        endif!(olsm_gbl)then

        write(6,*)"igd,jgd,lgtss=",igd,jgd,datan(ijgd)
        write(6,*)"igd,jgd,ogtss=",igd,jgd,datan(ijgd+ix*iy)


!not done since tss already at sea level     write(6,*)"tss at sea level here"

        write(6,*)"=========================> now interp. land data"
        call sintp16(datan(1:ix*iy),ix,iy,sfct,glon,glat,sdiag,il)                 ! land

        if(olsm_gbl .and. nopnts.gt.0)then
          write(6,*)"=========================> now interp. ocean data"
           call sintp16(datan(1+ix*iy:2*ix*iy),ix,iy,sfcto_m,glon,glat,
     &                  sdiag,il)   ! ocean
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
        do j=1,jl
         do i=1,il
          iq=i+(j-1)*il
          sfct(iq)=0.
         enddo ! i
        enddo ! j

      endif ! ier eq 0 , sfct

      call prt_pan(sfct,il,jl,2,'sfct')
      !call prt_pan(sfct,il,jl,1,'sfct')

      write(6,*)"============================================fracice"

      fracice=-1.
      ier = nf_inq_varid(ncid,'fracice',idvar)
      write(6,*)"ier=",ier," idvar=",idvar

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
        call fill(datan(1+ix*iy:2*ix*iy),ix,iy,.1*spval,
     &            datan(1+2*ix*iy:3*ix*iy))

        write(6,*)"igd,jgd,ogtss=",igd,jgd,datan(ijgd+ix*iy)

        write(6,*)"=========================> now interp. ocean data"
        call sintp16(datan(1+ix*iy:2*ix*iy),ix,iy,fracice,glon,glat,
     &                  sdiag,il)   ! ocean

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
        fracice=-1.

      endif ! ier eq 0 , sfct

      call prt_pan(fracice,il,jl,2,'fracice')

!############################################################################
! end sfc data
!############################################################################

      write(6,*)"check of temp data to ensure all is going okay"
      write(6,*)" findxn model temp(1)"
      call findxn(temp(:,:,1),ifull,-1.e29,xa,kx,an,kn)
      write(6,*)" findxn model temp(nplev)"
      call findxn(temp(:,:,nplev),ifull,-1.e29,xa,kx,an,kn)

      write(6,*)"nplev=",nplev
c constrain rh to 0-100
        do k=1,nplev
         do j=1,jl
          do i=1,il
           rh(i,j,k)=min(100.,max(0.,rh(i,j,k)))
          enddo ! i
         enddo ! j
        enddo ! k

!############### fix winds if CC grid ###############################
       if ( ogbl ) then

c     here use unstaggered lats and lons for u and v
c     For calculating zonal and meridional wind components, use the
c     following information, where theta is the angle between the
c     (ax,ay,az) vector [along the xg axis] and the zonal-component-vector:
c     veczon = k x r, i.e. (-y,x,0)/sqrt(x**2 + y**2)
c     vecmer = r x veczon, i.e. (-xz,-yz,x**2 + y**2)/sqrt(x**2 + y**2)
c     costh is (veczon . a) = (-y*ax + x*ay)/sqrt(x**2 + y**2)
c     sinth is (vecmer . a) = [-xz*ax - yz*ay + (x**2 + y**2)*az]/sqrt
c      using (r . a)=0, sinth collapses to az/sqrt(x**2 + y**2)
c     For rotated coordinated version, see JMcG's notes

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

      call maxmin(u,' u',0,1.,il,kl)
      call maxmin(v,' v',0,1.,il,kl)
      !call maxmin(hgt,' hgt',0,.001,il,kl)
      !call maxmin(temp,' temp',0,1.,il,kl)
      !call maxmin(rh,' rh',0,1000.,il,kl)

      !call prt_pan(u(1,1, 1),il,jl,2,'u : 1')
      !call prt_pan(v(1,1, 1),il,jl,2,'v : 1')

      imidpan2 = il/2+(jk+jl/2-1)*il
      do k=1,nplev
        write(6,*)'k,u/v(imidpan2,1,k)',u(imidpan2,1,k),v(imidpan2,1,k)
      enddo

      write(6,*)"convert winds to CCAM grid convention"

      cx=-1.e29
      sx=-1.e29
      cn= 1.e29
      sn= 1.e29


      !do iq=1,ifull
      do j=1,jl
        do i=1,il
        iq=i+(j-1)*il
c       set up unit zonal vector components
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
           uzon = u(iq,1,k)
           vmer = v(iq,1,k)
           u(i,j,k)= costh*uzon+sinth*vmer
           v(i,j,k)=-sinth*uzon+costh*vmer
           if(iq.eq.imidpan2)then
             write(6,'("before zon/mer; k,u,v: ",i3,2f10.2)')k,uzon,vmer
             write(6,'("zonx,zony,zonz,den,costh,sinth",
     &                6f8.4)')zonx,zony,zonz,den,costh,sinth
             write(6,'("after zon/mer; k,u,v: ",i3,2f10.2)')
     &                        k,u(i,j,k),v(i,j,k)
           endif
        enddo  ! k loop
        end do
      end do
      !enddo      ! iq loop

      write(6,*)'cx,cn,sx,sn=',cx,cn,sx,sn
      write(6,*)'after zon/meridional'

      call maxmin(u(:,:,1:kl),' u',0,1.,il,kl)
      call maxmin(v(:,:,1:kl),' v',0,1.,il,kl)
      call maxmin(hgt(:,:,1:kl),' hgt',0,.001,il,kl)
      call maxmin(temp(:,:,1:kl),' temp',0,1.,il,kl)
      call maxmin(rh(:,:,1:kl),' rh',0,1.,il,kl)

      !call prt_pan(u(1,1, 1),il,jl,2,'u : 1')
      !call prt_pan(v(1,1, 1),il,jl,2,'v : 1')

!############### fix winds if DARLAM grid ###############################
       else ! not ogbl

c convert e-w/n-s lat/lon winds to model winds
c loop over all model grid points
        do k=1,nplev
         write(6,*)k,temp(1,1,k),u(1,1,k),v(1,1,k)
         do j=1,jl
          do i=1,il
c get lat lon of model grid ( just to get ther )
           call lconll(rlon,rlat,float(i),float(j))
c calculate ucmp l.c.winds
c ulc=v@u*s(th)+u@u*c(th)
           ull = u(i,j,k)
           vll = v(i,j,k)
           u(i,j,k)=vll*sin(ther)+ull*cos(ther)
c calculate vcmp l.c.winds
c vlc=v@v*c(th)-u@v*s(th)
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
      do k=1,nplev
          cplev(k)=plev(nplev+1-k)
          write(6,'(5f10.2,f10.5)') cplev(k),hgt(i,j,k),temp(i,j,k)
     &                       ,u(i,j,k),v(i,j,k),rh(i,j,k)
      enddo ! k=1,nplev

      iq=il/2+(jl/2-1)*il
      write(6,'("pmsl=",f12.2," sfct=",f12.2)') pmsl(iq),sfct(iq)

      ! intepolate valid level
      call amap ( validlevhost, ix, iy, 'gbl levl', 0., 0. )
      call sintp16(validlevhost,ix,iy,validlevcc,glon,glat,sdiag,il)

      write(6,*)"calling vidar now!! ntimes,iarch=",ntimes,iarch

      if2=0

!#######################################################################
      call vidar(nplev,hgt,temp,u,v,rh,validlevcc
     &     ,iyr,imn,idy,ihr,iarch,time,mtimer,cplev,io_out,il,kl)
!#######################################################################

      enddo ! narch

      deallocate(datan,zs_gbl,lsm_gbl)
      deallocate(glon,glat)

      write(6,*)'*********** Finished cdfvidar ************************'
      
      deallocate(x,y,z,ax,ay,az,bx,by,bz)
      deallocate(hgt,temp,u,v,rh)
      deallocate(sfcto_m,lsmg_m)
      deallocate(zsg)
      call sigdatadealloc    
      call latlongdealloc      
      call clldealloc
      call comsigdealloc

      stop
      end ! cdfvidar
c***************************************************************************
      subroutine ncread_2d(idhist,iarch,idvar,il,jl,var)

!     include 'gblparm.h'
      include 'netcdf.inc'

      integer start(3),count(3)

      real var(il*jl), addoff, sf

      integer*2, dimension(:), allocatable :: ivar
!     integer*2 ivar(nnx*nny)

      character*30 name

      write(6,*)"=============================> ncread_2d idhist=",idhist
      write(6,*)"iarch=",iarch," idvar=",idvar
      write(6,*)"il=",il," jl=",jl

c read name
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
         call ncvgt(idhist,idvar,start,count,ivar,ier)
         write(6,*)"ivar(1)=",ivar(1)," ier=",ier
         write(6,*)"ivar(il*jl)=",ivar(il*jl)
      else if ( itype .eq. nf_float ) then
         write(6,*)"variable is float"
         call ncvgt(idhist,idvar,start,count,var,ier)
         write(6,*)"var(1)=",var(1)," ier=",ier
         write(6,*)"var(il*jl)=",var(il*jl)
      else
         write(6,*)"variable is unknown"
         stop
      endif

c obtain scaling factors and offsets from attributes
        call ncagt(idhist,idvar,'add_offset',addoff,ier)
        if ( ier.ne.0 ) addoff=0.
        write(6,*)"ier=",ier," addoff=",addoff

        call ncagt(idhist,idvar,'scale_factor',sf,ier)
        if ( ier.ne.0 ) sf=1.
        write(6,*)"ier=",ier," addoff=",addoff

      else!(ier.eq.0)then
c no data found
        do i=1,il*jl
         var(i)=0
        enddo
        sf=0.
        addoff=0.
      endif!(ier.eq.0)then

c unpack data
      dx=-1.e29
      dn= 1.e29
      do j=1,jl
        do i=1,il
          ij=i+(j-1)*il
      	  if ( itype .eq. nf_short ) then
           if(i.eq.1.and.j.eq.1)
     &      write(6,*)"ivar,sf,addoff=",ivar(ij),sf,addoff
            var(ij) = ivar(ij)*sf + addoff
          else
           if(i.eq.1.and.j.eq.1)
     &      write(6,*)"var,sf,addoff=",var(ij),sf,addoff
            var(ij) = var(ij)*sf + addoff
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

      return ! ncread_2d
      end
c***************************************************************************
      subroutine ncread_3d(idhist,iarch,idvar,il,jl,kl,var)
c             call ncread_3d(ncid,iarch,idvar,ix,iy,nplev,datan)

!     include 'gblparm.h'
      include 'netcdf.inc'

      integer start(4),count(4)

!     integer*2 ivar(nmax*35)
      integer*2, dimension(:), allocatable :: ivar
      double precision, dimension(:), allocatable :: dvar

      real var(il*jl*kl)
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

c read data
      write(6,*)"idhist=",idhist," idvar=",idvar
      ier = nf_inq_vartype(idhist,idvar,itype)
      write(6,*)"ier=",ier," itype=",itype

      if ( itype .eq. nf_short ) then
         write(6,*)"variable is short"
	 if (.not.allocated(ivar)) allocate(ivar(il*jl*kl))
	 if (size(ivar).ne.il*jl*kl) then
	   deallocate(ivar)
	   allocate(ivar(il*jl*kl))
	 end if
         call ncvgt(idhist,idvar,start,count,ivar,ier)
      else if ( itype .eq. nf_float ) then
         write(6,*)"variable is float"
         call ncvgt(idhist,idvar,start,count,var,ier)
      else if ( itype .eq. nf_double ) then
         write(6,*)"variable is double"
	 if (.not.allocated(dvar)) allocate(dvar(il*jl*kl))
	 if (size(dvar).ne.il*jl*kl) then
	   deallocate(dvar)
	   allocate(dvar(il*jl*kl))
	 end if
         call ncvgt(idhist,idvar,start,count,dvar,ier)
	 var=real(dvar)
      else
         write(6,*)"variable is unknown"
         stop
      endif

      addoff=0.
      sf=1.
      if ( itype .eq. nf_short ) then
      write(6,*)"obtain scaling factors and offsets from attributes"
      call ncagt(idhist,idvar,'add_offset',addoff,ier)
      if ( ier.ne.0 ) addoff=0.
      write(6,*)"ier=",ier," addoff=",addoff

      call ncagt(idhist,idvar,'scale_factor',sf,ier)
      if ( ier.ne.0 ) sf=1.
      write(6,*)"ier=",ier," sf=",sf
      endif

c unpack data
      dx=-1.e29
      dn= 1.e29
      do k=1,kl
       do j=1,jl
        do i=1,il
          ijk=i+(j-1)*il+(k-1)*il*jl
          if(i.eq.1.and.j.eq.1.and.k.eq.1)
     &       write(6,*)"i,j,k,ijk=",i,j,k,ijk
      	  if ( itype .eq. nf_short ) then
           if(i.eq.1.and.j.eq.1.and.k.eq.1)
     &      write(6,*)"ivar,sf,addoff=",ivar(ijk),sf,addoff
            var(ijk) = ivar(ijk)*sf + addoff
          else
           if(i.eq.1.and.j.eq.1.and.k.eq.1)
     &      write(6,*)"var,sf,addoff=",var(ijk),sf,addoff
            var(ijk) = var(ijk)*sf + addoff
          endif
          if(i.eq.1.and.j.eq.1.and.k.eq.1)
     &      write(6,*)"var=",var(ijk)
          dx=max(dx,var(ijk))
          dn=min(dn,var(ijk))
        end do
       end do
      end do

      write(6,*)"ncread_3d idvar=",idvar," iarch=",iarch
      write(6,*)"ncread_3d dx=",dx," dn=",dn

      return ! ncread_3d
      end
c***********************************************************************
      subroutine filt_nc(var,il,jl,kl)

      real var(il,jl,kl)

      write(6,*) "filt_nc"

!     do k=1,kl
!      do j=1,jl
!       do i=1,il
!         var(i,j,k) = var(i,j,k)
!       end do
!      end do
!     end do

      return ! filt_nc
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
      
      do j=1,iy
        do i=1,ix
          iq=i+ix*(j-1)
          do k=nint(validlev(iq)),plev	
	    iqk=i+ix*(j-1)+ix*iy*(k-1)
	    if (datan(iqk).lt.1.E10) then
	      validlev(iq)=real(k)
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
      real, dimension(ix*iy*plev), intent(inout) :: datan
      real, dimension(ix*iy) :: datatemp

      do k=1,plev
        is=ix*iy*(k-1)+1
        ie=ix*iy*k
        if (all(datan(is:ie).gt.1.E10)) then
          write(6,*) "ERROR: No valid data on level"
          stop
        end if
        datatemp(:)=datan(is:ie)
        do while (any(datan(is:ie).gt.1.E10))
          do j=1,iy
            do i=1,ix
              iqk=i+ix*(j-1)+ix*iy*(k-1)
              if (datan(iqk).lt.1.E10) then
                if (j.lt.iy) then
                  iqn=i+ix*j
                  if (datatemp(iqn).gt.1.E10) then
                    datatemp(iqn)=datan(iqk)
                  end if
                end if
                if (i.lt.ix) then
                  iqn=i+1+ix*(j-1)
		else
                  iqn=1+ix*(j-1)	
		end if
                if (datatemp(iqn).gt.1.E10) then
                    datatemp(iqn)=datan(iqk)
                end if
	        if (j.gt.1) then
                  iqn=i+ix*(j-2)
                  if (datatemp(iqn).gt.1.E10) then
                    datatemp(iqn)=datan(iqk)
                  end if
		end if
	        if (i.gt.1) then
                  iqn=i-1+ix*(j-1)
		else
		  iqn=ix+ix*(j-1)
		end if
                if (datatemp(iqn).gt.1.E10) then
                  datatemp(iqn)=datan(iqk)
                end if
	      end if
	    end do
	  end do
          datan(is:ie)=datatemp(:)
        end do
      end do
      
      return
      end subroutine filldat
