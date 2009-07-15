      program one

      common/mapproj/du,tanl,rnml,stl1,stl2

      !parameter (nnx=360,nny=181,nmax=nnx*nny)  ! global, 1.0 degree
      include 'gblparm.h'

      real datavn(nmax)
      real datlnd(nmax)
      real avnlnd(nmax)
      real avnzs(nmax)

      include 'newmpar.h'
      include 'netcdf.inc'

      include 'parm.h'
      include 'xyzinfo.h'  ! x,y,z,rlat,rlong,wts
      include 'latlong.h'  ! x,y,z,rlat,rlong,wts
      include 'vecsuv.h'
      parameter ( pi=3.1415926536 )
      parameter ( g=9.80616 )

      real datl(il*jl)
      real lsmod(il*jl)
      real dat(ifull),zsmod(ifull)
      real tmp(ifull)

      logical opos,occam,ols
      character*80 zsavn,lsavn
      character*80 gbldat,zsfil,moddat

      namelist/gnml/igbldat,gbldat,imoddat,moddat
     &     ,inzsavn,zsavn,inlsavn,lsavn,inzs,zsfil
     &     ,ols,opos,occam
     &     ,igd,jgd,igd1,igd2,jgd1,jgd2
     &     ,id,jd,id1,id2,jd1,jd2

      data igd/1/,jgd/1/,id/1/,jd/1/
      data id1/30/,id2/40/,jd1/85/,jd2/90/
      data igd1/30/,igd2/40/,jgd1/85/,jgd2/90/
      data ix/nnx/,iy/nny/
      data ols/.false./
      data occam/.true./,opos/.false./
      data inzsavn/11/ , zsavn/'zsavn.ff'/
      data inlsavn/12/ , lsavn/'lsavn.ff'/
      data igbldat/13/ , gbldat/'gbldat.ff'/
      data inzs   /14/ , zsfil/'topog5'/

      save

      grdx=1.
      slon=0.
      grdy=-1.
      slat=90.
      npts=il*jl

      write(6,*)'read namelist'
      read (5, gnml)
      write(6,nml=gnml)

      write(6,*)'open ',inzs,' zsfil=',zsfil

      if ( occam ) then
        write(6,*)"set up cc geometry"
        open(inzs,file=zsfil,status='old')
        read(inzs,'(i3,i4,2f6.1,f5.2,f9.0,a47)')
     &          ilx,jlx,rlong0,rlat0,schmidt,ds,header
        du=rlong0
        tanl=rlat0
        rnml=schmidt
        write(6,*)"gbl mapproj=",ilx,jlx,rlong0,rlat0,schmidt,ds,header
        if(ilx.ne.il.or.jlx.ne.jl)
     &     stop 'wrong topo file supplied (il,jl) for cc model grid'
        call setxyz
        rlatx=-1.e29
        rlatn= 1.e29
        rlonx=-1.e29
        rlonn= 1.e29
        do iq=1,ifull
c       convert conformal cubic lats & longs to degrees (-90 to 90) & (0 to 360)
c       used in sint16; N.B. original rlong is -pi to pi
          rlat(iq)=rlat(iq)*180./pi
          rlong(iq)=rlong(iq)*180./pi
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
      else ! ( occam ) then
        open(inzs,file=zsfil,form='formatted',recl=il*7,status='old')
        write(6,*)'read zsfil header'
        read(inzs,*,err=25)ik,jk,ds,du,tanl,rnml,stl1,stl2
 25     if(ik.eq.0.or.jk.eq.0)then
           write(6,*)'no header in newtopo file'
        else
           write(6,*)'Header information for topofile'
           write(6,*)'ik,jk,ds,du,tanl,rnml,stl1,stl2'
     &           ,ik,jk,ds,du,tanl,rnml,stl1,stl2
           if(ik.ne.il.or.jk.ne.jl)stop 'wrong topofile supplied'
        endif     ! (ik.eq.0.or.jk.eq.0)
        write(6,*)"set up model grid params by calling lconset ds=",ds
        call lconset(ds)
      endif ! ( occam ) then

      write(6,*)'read model grid zs (g*zs)'
      read(inzs,*)zsmod
      write(6,*)'convert g*zs to zs(m)'
      do iq=1,ifull
        zsmod(iq)=zsmod(iq)/g
      enddo !iq=1,ifull
      write(6,*)'read model grid land-sea mask (0=ocean, 1=land)'
      read(inzs,*)lsmod ! this is zs
      !if(inzsavn.le.0)read(inzs,*)lsmod ! extra read is zs not read above

      ijd=id+il*(jd-1)
      write(6,*)"ijd=",ijd," zsmod=",zsmod(ijd)," lsmod=",lsmod(ijd)

      write(6,*)"zsmod on model grid"
      do j=jd1,jd2
       write(6,'(i4,20f10.2)')j,(zsmod(i+(j-1)*ix),i=id1,id2)
      enddo
      write(6,'(i4,20f10.2)')j,(real(i),i=id1,id2)

      write(6,*)"lsmod on model grid"
      do j=jd1,jd2
       write(6,'(i4,20f10.2)')j,(lsmod(i+(j-1)*ix),i=id1,id2)
      enddo
      write(6,'(i4,20f10.2)')j,(real(i),i=id1,id2)

      write(6,*)"Now deal with sfc. zs and land/sea mask"

      if(inzsavn.gt.0)then
        write(6,*)"open sfc. orog. for avn data"
        write(6,*)"zsavn=",inzsavn,zsavn
        open(inzsavn,file=zsavn,form='formatted',status='old')
! read in free formatted avn sfc. orography (28-mar-2000)
        rewind inzsavn
        read(inzsavn,*)((avnzs(i+(j-1)*ix),i=1,ix),j=1,iy)
        write(6,*)"note that it runs north to south"
        call amap ( avnzs, ix, iy, 'gbl sfczs(m)', 0., 200. )
        !write(6,*)"interp. zsavn to output grid"
        !call sintp16(avnzs,ix,iy,zs)
        !write(6,*) 'findxn model sfc.height (m)'
        !call findxn(zs,npts,-1.e29,xa,kx,an,kn)
        write(6,*) 'close unit inzsavn=',inzsavn
        close(inzsavn)
      endif

      write(6,*)"avnzs"
      do j=jgd2,jgd1,-1
       write(6,'(i4,20f10.2)')j,(avnzs(i+(j-1)*ix),i=igd1,igd2)
      enddo
      write(6,'(i4,20f10.2)')j,(real(i),i=igd1,igd2)

      write(6,*)'read land-sea mask (0=ocean, 1=land) inlsavn=',inlsavn
      if(inlsavn.gt.0)then
        write(6,*)"AVN land sea mask"
        write(6,*)"open sfc. orog. for avn data"
        write(6,*)"lsavn=",inlsavn,lsavn
        open(inlsavn,file=lsavn,form='formatted',status='old')
        read(inlsavn,*)((avnlnd(i+(j-1)*ix),i=1,ix),j=1,iy)
        write(6,*)"note that it runs north to south"
        call amap ( avnlnd, ix, iy, 'gbl land/sea mask', 0., 1. )
        close(inlsavn)
      else
        write(6,*)"######################## WARNING!!!!!!!!!!!!!!!!"
        write(6,*)"######################## setting input lsm == 0!"
        do j=1,iy
         do i=1,ix
          avnlnd(i+(j-1)*ix)=0.
         enddo
        enddo
      endif
      close(inzs)

      write(6,*)"avnlnd"
      do j=jgd2,jgd1,-1
       write(6,'(i4,20f10.2)')j,(avnlnd(i+(j-1)*ix),i=igd1,igd2)
      enddo
      write(6,'(i4,20f10.2)')j,(real(i),i=igd1,igd2)

      write(6,*)' reading variables '

c***********************************************************************

c sfc data
      write(6,*)'gbldat=',gbldat
      open(igbldat,file=gbldat,form='formatted',status='old')
      write(6,*)"note that it runs north to south",ix,iy
      read(igbldat,*)((datavn(i+(j-1)*ix),i=1,ix),j=1,iy)
      call amap ( datavn, ix, iy, 'gbl input', 0., 10. )

      write(6,*)"input gdata"
      do j=jgd2,jgd1,-1
       write(6,'(i4,20f10.2)')j,(datavn(i+(j-1)*ix),i=igd1,igd2)
      enddo
      write(6,'(i4,20f10.2)')j,(real(i),i=igd1,igd2)

      spval=-1.e10
      write(6,*)"spval=",spval
      write(6,*)"inlsavn=",inlsavn
      ijgd=igd+ix*(jgd-1)
      write(6,*)"igd,jgd,ijgd=",igd,jgd,ijgd
      ijd=id+il*(jd-1)
      write(6,*)"id,jd,ijd=",id,jd,ijd

!-----------------------------------------------------------------------
      if ( ols ) then
!-----------------------------------------------------------------------

      write(6,*)"prepare to interp. dat for sea and land separately"
      write(6,*)"igd,jgd,gdat=",igd,jgd,datavn(ijgd)
      do j=1,iy
       do i=1,ix
        iq = i+(j-1)*ix
!       write(6,*)i,j,iq,datavn(iq) ,avnlnd(iq)
        datlnd(iq)=datavn(iq)
        if(inlsavn.gt.0)then
        if ( avnlnd(iq) .gt. .5 ) then
! land
          datavn(iq)=spval
        else !!!  ( avnlnd(iq) .gt. .5 ) then
! ocean
          datlnd(iq)=spval
        endif ! ( avnlnd(iq) .gt. .5 ) then
        endif!(inlsavn.gt.0)then
       enddo ! ix
      enddo ! iy
      write(6,*) "two global dat arrays with spval=", spval
      write(6,*)"igd,jgd,ogdat=",igd,jgd,datavn(ijgd)
      write(6,*)"igd,jgd,lgdat=",igd,jgd,datlnd(ijgd)

      write(6,*)"input water gdata"
      do j=jgd2,jgd1,-1
       write(6,'(i4,20f10.2)')j,(datavn(i+(j-1)*ix),i=igd1,igd2)
      enddo
      write(6,'(i4,20f10.2)')j,(real(i),i=igd1,igd2)

      write(6,*)"input land gdata"
      do j=jgd2,jgd1,-1
       write(6,'(i4,20f10.2)')j,(datlnd(i+(j-1)*ix),i=igd1,igd2)
      enddo
      write(6,'(i4,20f10.2)')j,(real(i),i=igd1,igd2)

      if(inlsavn.gt.0)then
         write(6,*)"fill in missing dat land values"
         call fill(datlnd,ix,iy,.1*spval,tmp)
      endif!(inlsavn.gt.0)then
      write(6,*)"fill in missing water values"
      call fill(datavn,ix,iy,.1*spval,tmp)

      write(6,*)"filled input water gdata"
      do j=jgd2,jgd1,-1
       write(6,'(i4,20f10.2)')j,(datavn(i+(j-1)*ix),i=igd1,igd2)
      enddo
      write(6,'(i4,20f10.2)')j,(real(i),i=igd1,igd2)

      write(6,*)"filled input land gdata"
      do j=jgd2,jgd1,-1
       write(6,'(i4,20f10.2)')j,(datlnd(i+(j-1)*ix),i=igd1,igd2)
      enddo
      write(6,'(i4,20f10.2)')j,(real(i),i=igd1,igd2)

      write(6,*)"igd,jgd,ogdat=",igd,jgd,datavn(ijgd)
      write(6,*)"igd,jgd,lgdat=",igd,jgd,datlnd(ijgd)

!-----------------------------------------------------------------------
      endif !( ols ) then
!-----------------------------------------------------------------------

      write(6,*)"now interp. land and ocean points separately"
      call sintp16(datavn,ix,iy,dat)
      write(6,*)"id,jd,odat=",id,jd,dat(ijd)

      if ( opos ) then
         do iq=1,ifull
           dat(iq)=max(0.,dat(iq))
         enddo ! i
      endif ! ( opos ) then

      write(6,*)"filled dat on model water grid"
      do j=jd1,jd2
       write(6,'(i4,20f10.2)')j,(dat (i+(j-1)*il),i=id1,id2)
      enddo
      write(6,'(i4,20f10.2)')j,(real(i),i=id1,id2)

!-----------------------------------------------------------------------
      if ( ols ) then
!-----------------------------------------------------------------------
      if(inlsavn.gt.0)then
         write(6,*)"interpolate dat land values"
         call sintp16(datlnd,ix,iy,datl)
      endif!(inlsavn.gt.0)then
      write(6,*)"id,jd,ldat=",id,jd,datl(ijd)

      write(6,*)"filled dat on model land grid"
      do j=jd1,jd2
       write(6,'(i4,20f10.2)')j,(datl(i+(j-1)*il),i=id1,id2)
      enddo
      write(6,'(i4,20f10.2)')j,(real(i),i=id1,id2)

      if(inlsavn.gt.0)then
        write(6,*)"now recombine two (land/ocean) fields"
        do j=1,jl
         do i=1,il
          iq=i+(j-1)*il
          if(lsmod(iq).gt..5) dat(iq)=datl(iq)
         enddo ! i
        enddo ! j
      endif!(inlsavn.gt.0)then
      write(6,*)"id,jd,dat=",id,jd,dat(ijd)
!-----------------------------------------------------------------------
      endif !( ols ) then
!-----------------------------------------------------------------------

      write(6,*)" findxn model dat"
      call findxn(dat,npts,-1.e29,xa,kx,an,kn)
      call amap ( dat, il, jl, 'model dat', 0., 0. )
      is = 1+(48*il)
      call amap ( dat(is), il, il, 'HR model dat', 0., 0. )

      write(6,*)"dat on model grid"
      do j=jd1,jd2
       write(6,'(i4,20f10.2)')j,(dat(i+(j-1)*il),i=id1,id2)
      enddo
      write(6,'(i4,20f10.2)')j,(real(i),i=id1,id2)

      open(imoddat,file=moddat,form='formatted',status='new')
      write(imoddat,*) dat

      write(6,*)'*****************************************'
      write(6,*)'successfully completed one'
      write(6,*)'*****************************************'

      stop
      end
c***************************************************************************
