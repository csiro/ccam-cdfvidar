      subroutine sintp16(gbl_in,ngxi,ngyi,model_dat,sdiag)

      include 'newmpar.h'
      include 'xyzinfo.h'  ! x,y,z,wts
      include 'latlong.h'  ! rlat,rlong arrays of lat and lon on model grid
                           ! slat,slon,dlat,dlon for input data
      include 'vecsuv.h'
      include 'gblparm.h' ! nnx,nny
      common/glonlat/glon(nnx),glat(nny)

c assumes model_dat grid dimensions come through 'newmpar.h'
c as il, jl
c also assumes that you have already called lconset(ds)
c needs lconset routine.

c the gbl_in grid is dimensioned ng1 by ngyi
c with first data point at 0 deg lon, 90 deg N lat   ******  LH coord system!!

c southern lats are negative

      parameter ( np2=nnx+2 )
      !parameter ( nxn=np2+2 )

      real gbl_in(ngxi,ngyi)

      real globex(-1:np2,ngyi)

      !real globex(-1:(nnx+2),ngyi)
      !real globex(nxn*ngyi) ! max extended array

      real model_dat(il,jl)

      logical sdiag,pdiag

!     data pdiag/.false./
      data pdiag/.true./

      if ( pdiag ) then
	      write(6,*)"sintp16 nnx,nny=",nnx,nny
	      write(6,*)"sintp16 ngxi,ngyi=",ngxi,ngyi
      endif

c     extend gbl_in array over Grenwich meridion
      dx=-1.e29
      dn=1.e29
      do j=1,ngyi
       do ii=1,ngxi
        globex(ii,j)=gbl_in(ii,j)
        dn=min(dn,gbl_in(ii,j))
        dx=max(dx,gbl_in(ii,j))
       enddo ! ii
       !globex(    -1,j)=gbl_in(nnx-1,j)
       !globex(     0,j)=gbl_in(nnx  ,j)
       globex(    -1,j)=gbl_in(ngxi-1,j)
       globex(     0,j)=gbl_in(  ngxi,j)
       globex(ngxi+1,j)=gbl_in(     1,j)
       globex(ngxi+2,j)=gbl_in(     2,j)
      enddo ! j

      if ( pdiag ) then

      write(6,*) 'sintp16 dn,dx=',dn,dx
      jdi = (ngyi*3)/4
      write(6,*)"jdi=",jdi
      write(6,"('gbl_in(i=1:2,jdi) =',33x,6f12.2)")
     &          (globex(i,jdi),i=1,2)
      write(6,"('globex(i=-1:2,jdi)=',9x,6f12.2)")
     &          (globex(i,jdi),i=-1,2)
      write(6,"('gbl_in(i=ngxi-2:ngxi  ,jdi)=',7f12.2)")
     &          (globex(i,jdi),i=ngxi-1,ngxi)
      write(6,"('globex(i=ngxi-2:ngxi+2,jdi)=',7f12.2)")
     &          (globex(i,jdi),i=ngxi-1,ngxi+2)

      write(6,*)"rlong(1,ifull),dlon=",rlong(1),rlong(ifull),dlon
      write(6,*)"rlat(1,ifull),dlat=",rlat(1), rlat(ifull),dlat

      write(6,*)"np2,ngxi,ngyi=",np2,ngxi,ngyi
      write(6,*)"slon,slat=",slon,slat
      write(6,*)"dlon,dlat=",dlon,dlat

! these are fractional grid index of this model point in input grid
      do i = 1,ifull,ifull-1
        xx = 1.+(rlong(i)-slon)/dlon
        if ( xx.gt.ngxi ) xx=xx-ngxi
        yy = 1.+( rlat(i)-slat)/dlat
        write(6,"('(',i6,')lon,lat=',2f10.2,' xx,yy,md=',3f10.2)")
     &            i,rlong(i),rlat(i),xx,yy,globex(int(xx),int(yy))
      enddo ! i

      endif ! ( pdiag ) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop over all model grid points
      !write(6,*)"loop over all model points ifull=",ifull
      dn=1.e29
      dx=-1.e29
      do iq=1,ifull
      !do iq=1,ifull,500
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       write(6,*)"iq,dlon,dlat=",iq,dlon,dlat
!       write(6,'("rlong,rlat,dlon,dlat=",4f8.2)')
!    &             rlong(iq),rlat(iq),dlon,dlat

! calc nearest fractional i,j of gbl_in array
        xx = 1.+(rlong(iq)-slon)/dlon
        if ( xx.gt.ngxi ) xx=xx-ngxi
        yy = 1.+( rlat(iq)-slat)/dlat
!       write(6,"('(',i6,')lon,lat=',2f10.2,' xx,yy=',3f10.2)")
!    &            iq,rlong(iq),rlat(iq),xx,yy,globex(int(xx),int(yy))

        !if(mod(iq-1,48).eq.0)then
        if(pdiag .and. mod(iq-1,2304).eq.0)then
              write(6,'("iq,lo,la,x,y=",i8,4f8.1)')
     &                   iq,rlong(iq),rlat(iq),xx,yy
        endif

        !write(6,*)"call intp16 to interpolate from input grid to model"

        !write(6,*)"intp16 np2,ngxi,ngyi,xx,yy=",np2,ngxi,ngyi,xx,yy

! actual interpolation done here

        call intp16(globex,np2,ngxi,ngyi,xx,yy,model_dat(iq,1),.false.)

        dn=min(dn,model_dat(iq,1))
        dx=max(dx,model_dat(iq,1))

c       model_dat(i,j)=max(xmin,model_dat(i,j))

        !write(6,*)"output model_dat=",model_dat(iq,1)

!     if(ngxi.lt.3.and.abs(model_dat(iq,1)).gt.100.)then
!       write(6,*) 'iq,rlong,rlat ',iq,rlong(iq),rlat(iq)
!     endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo ! iq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !write(6,*) 'done sintp16 dn,dx=',dn,dx

      !pdiag=.false.

      return
      end
!#######################################################################

      subroutine intp16(globex,np2,ngxi,ngyi,x,y,model_dat,ofirst)

c**********************************************************************
c
c subroutine intp16 interpolates from gbl_in gridpoints to lat-lon
c     location specified by the main program, via a 16- point bessel
c     interpolation scheme.  interpolated
c     lat/lon grid values are contained in the variable "model_dat"
c     jlm: this is a quadratic fitting pts 2 & 3 & cubic value of centre point
c**********************************************************************

      logical prnt,ofirst,odiag
      real globex(-1:np2,ngyi)
      real model_dat

      data prnt/.false./
      data p25/.25/, c1/1./

      !write(6,*)"intp16 np2,ngxi,ngyi,x,y=",np2,ngxi,ngyi,x,y

c determine the interpolation distances.

      ix=x
      ixm1=ix-1
      ixp1=ix+1
      ixp2=ix+2
      !write(6,*)"ixm1,ix,ixp1,ixp2=",ixm1,ix,ixp1,ixp2

      iy=y
      iym1=iy-1
      iyp1=iy+1
      iyp2=iy+2
      !write(6,*)"iym1,iy,iyp1,iyp2=",iym1,iy,iyp1,iyp2

      !write(6,*)"ix,ngxi-2,iy,ngyi/2-1=",ix,ngxi-2,iy,ngyi/2-1

      odiag = (ix.ge.ngxi-2 .or. ix.le.1) .and.
     &        (iy.eq.ngyi/2)

      !write(6,*)"odiag=",odiag

      if(odiag)write(6,'("ngxi, ngyi, x, y, ofirst",2i8,2f8.1,l2)')
     &           ngxi, ngyi, x, y, ofirst

      odiag = .false.
      if(iym1.gt.ngyi)then
        iym1=ngyi
        iy  =ngyi
        iyp1=ngyi
        iyp2=ngyi
        odiag = .false.
      endif
      if(iy.gt.ngyi)then
        iy  =ngyi
        iyp1=ngyi
        iyp2=ngyi
        odiag = .false.
      endif
      if(iyp2.lt.1)then
        iym1=1
        iy  =1
        iyp1=1
        iyp2=1
        odiag = .true.
      endif
      if(iyp1.lt.1)then
        iym1=1
        iy  =1
        iyp1=1
        odiag = .true.
      endif
      if(iy  .lt.1)then
        iym1=1
        iy  =1
        odiag = .true.
      endif
      if(iym1.lt.1)then
        iym1=1
        odiag = .true.
      endif
      if ( iym1.lt.1 .or. iyp2.gt.ngyi ) then
        if ( ofirst) then
          write(6,'("ix,iy,iym1,iyp2,ngyi=",5i6)')ix,iy,iym1,iyp2,ngyi
          write(6,*) '*********** n out of bounds in intp16 **********'
        endif ! ofirst
        if ( iym1.lt.  1 ) iym1=1
        if ( iyp1.gt.ngyi ) iyp1=ngyi
        if ( iyp2.gt.ngyi ) iyp2=ngyi
        odiag = .false.
      endif
      odiag = .false.

      dx=x-ix
      dy=y-iy
      dxx=p25*(dx-c1)
      dyy=p25*(dy-c1)

c determine the 16 x-y gridpoints, i.e. top11-top44, to be used
c        for the interpolation.

      top11 = globex(ixm1,iym1)
      top12 = globex(ixm1,iy  )
      top13 = globex(ixm1,iyp1)
      top14 = globex(ixm1,iyp2)
      top21 = globex(ix  ,iym1)
      top22 = globex(ix  ,iy  )
      top23 = globex(ix  ,iyp1)
      top24 = globex(ix  ,iyp2)
      top31 = globex(ixp1,iym1)
      top32 = globex(ixp1,iy  )
      top33 = globex(ixp1,iyp1)
      top34 = globex(ixp1,iyp2)
      top41 = globex(ixp2,iym1)
      top42 = globex(ixp2,iy  )
      top43 = globex(ixp2,iyp1)
      top44 = globex(ixp2,iyp2)

c perform the 16-point interpolation.

      aa = top21 + dx*(top31-top21 + dxx*(top11-top21-top31+top41) )
      ab = top22 + dx*(top32-top22 + dxx*(top12-top22-top32+top42) )
      ac = top23 + dx*(top33-top23 + dxx*(top13-top23-top33+top43) )
      ad = top24 + dx*(top34-top24 + dxx*(top14-top24-top34+top44) )
      model_dat = ab + dy*(ac-ab + dyy*(aa-ab-ac+ad) )

      !if(x.gt.-1.e29)then

      odiag=odiag.or.(abs(model_dat).gt.1.e29)
      !odiag=odiag.or.(abs(model_dat).lt.1.e-10)
      !odiag=odiag .or. abs(top44-top11).lt.1.e-10
      if(abs(model_dat).gt.1.e29)write(6,*)"model_dat.gt.1.e29"
      !if(abs(model_dat).lt.1.e-10)write(6,*)"model_dat.lt.1.e-10"
      !if(abs(top44-top11).lt.1.e-10)write(6,*)"top11=top44*************"

      if(odiag)then
         write(6,'("intp16 np2,ngxi,ngyi,x,y=",3i6,2f10.3)')
     &                     np2,ngxi,ngyi,x,y
	      jdi = (ngyi*3)/4
	      write(6,*)"jdi=",jdi
	      write(6,"('gbl_in(i=1:2,jdi) =',33x,6f12.2)")
     &          (globex(i,jdi),i=1,2)
	      write(6,"('globex(i=-1:2,jdi)=',9x,6f12.2)")
     &          (globex(i,jdi),i=-1,2)
	      write(6,"('gbl_in(i=ngxi-2:ngxi  ,jdi)=',7f12.2)")
     &          (globex(i,jdi),i=ngxi-1,ngxi)
	      write(6,"('globex(i=ngxi-2:ngxi+2,jdi)=',7f12.2)")
     &          (globex(i,jdi),i=ngxi-1,ngxi+2)
        write(6,*) 'x,y: ',x,y
        write(6,111) ix,iy,dx,dy,dxx,dyy
111     format (5x,'(ix,iy)=',i3,',',i3,3x,'(dx,dy)=',f6.3,',',f6.3,
     &          3x,'(dxx,dyy)=',f6.3,',',f6.3)
	write(6,112)ixm1,ix,ixp1,ixp2
112     format (2x,'top11-top44=....',/,(6x,4i12) )
        write(6,122) iyp2,top14,top24,top34,top44
     &		    ,iyp1,top13,top23,top33,top43
     &              ,iy  ,top12,top22,top32,top42
     &		    ,iym1,top11,top21,top31,top41
122     format (i6,4f12.3)
        write(6,133) aa,ab,ac,ad,model_dat
133     format (1h ,2x,'aa=',f10.3,3x,'ab=',f10.3,3x,'ac=',f10.3,3x,
     &             'ad=',f10.3,7x,'model_dat=',f10.3)
        write(6,*)"+++++++++++++++++++++++++++++++++++++++++++++"
        write(6,*)"+++++++++++++++++++++++++++++++++++++++++++++"
      endif

      return
      end
