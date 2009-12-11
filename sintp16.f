      subroutine sintp16(gdat,ng1,ng2,rdat,glon,glat,odiag,il)

      use latlong_m

      logical odiag,opdiag

! assumes that you have already called lconset(ds)
! needs lconset routine.
! the gdat grid is dimensioned ng1 by ng2
! with first data point at 0 deg lon, 90 deg N lat ****** LH coord system!!
! southern lats are negative

      !include 'newmpar.h' ! il,jl,ifull

      integer il,jl,ifull

      real gdat(ng1*ng2),glon(ng1),glat(ng2)
      real rdat(6*il*il)

      jl=6*il
      ifull=il*jl

      !include 'latlong.h'  ! rlat,rlong

      if(odiag)write(6,'(72("="))')

      if(odiag)write(6,*)"sintp16 g(1),ng1,ng2,il,jl="
      if(odiag)write(6,*) gdat(1),ng1,ng2,il,jl
      if(odiag)write(6,*)"glon(1,ng1/2,ng1)="
     &                   ,glon(i),glon(ng1/2),glon(ng1)
      if(odiag)write(6,*)"glat(1,ng2/2,ng2)="
     &                   ,glat(i),glat(ng2/2),glat(ng2)

      if(odiag)write(6,*) "il,jl,ifull=",il,jl,ifull
      if(odiag)write(6,*) "rlong(1/ifull)=",rlong(1),rlong(ifull)
      if(odiag)write(6,*) "rlat(1/ifull)=",rlat(1),rlat(ifull)

      gdx=-1.e29
      gdn=1.e29
      do j=1,ng2
       do i=1,ng1
        iq=i+(j-1)*ng1
        gdn=min(gdn,gdat(iq))
        gdx=max(gdx,gdat(iq))
       enddo ! i=1,ng1
      enddo ! j=1,ng2

      write(6,*)'sintp16 gbl grid:gdn,gdx,ifull=',gdn,gdx,ifull

      dn=1.e29
      dx=-1.e29

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! loop over all model grid points
! to interpolate from input grid to model grid
      do iq=1,ifull
        jmod=1+((iq-1)/il)
        imod=iq-(jmod-1)*il
        !if(mod(iq,il).eq.0)write(6,*)"iq,im,jm=",iq,imod,jmod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

! compute x and y of model point in input grid space

        xlon=rlong(iq)
        if(rlong(iq).lt.0.)xlon=360.+rlong(iq)

        x=(xlon*ng1/360.)+1. ! assumes x=1 at 0 deg
        y=-999 ! default (no value) case

        do jg=1,ng2-1
         jgp=jg+1

         if(jgp.eq.ng2)then  ! special case if at end of y grid

          inc=0
          y=float(ng2-1)+(rlat(iq)-glat(ng2-1))
     &                   /(glat(ng2)-glat(ng2-1))

         elseif(glat(1).lt.glat(ng2))then  ! glat incr./w j

          inc=1
          if(rlat(iq).ge.glat(jg).and.rlat(iq).lt.glat(jgp))then
           y=float(jg)+(rlat(iq)-glat(jg))/(glat(jgp)-glat(jg))
           go to 88
          endif!(rlat(iq).gt.glat(jg))then

         else ! (glat(1).gt.glat(ng2))then  ! glat decr./w j

          inc=-1
          if(rlat(iq).ge.glat(jgp).and.rlat(iq).lt.glat(jg))then
           y=float(jg)+(rlat(iq)-glat(jg))/(glat(jgp)-glat(jg))
           go to 88
          endif!(rlat(iq).gt.glat(jg))then

         endif!(glat(1).gt.glat(ng2)then

        enddo ! jg=1,ng2-1

 88     continue

        if ( y.lt.-900. ) stop

!       if(i.eq.20.and.j.gt.45.and.j.lt.55) then
!       write(30,'("iq,i,j,xlon,rlat(iq),x,y=",3i5,4f10.2)')
!    &              iq,i,j,xlon,rlat(iq),x,y
!       endif!(i.eq.20.and.j.gt.45.and.j.lt.55) then

!***********************************************************************
! acutal interpolation to a point

        opdiag=imod.eq.il/2.and.il.lt.jmod.and.jmod.lt.2*il
        opdiag=(jmod.eq.1.5*il)
        opdiag=F

        call intp16(gdat,ng1,ng2,x,y,rdat(iq),opdiag)

!***********************************************************************

        !if(odiag.and.iq.eq.20+(20-1)*il)then
        if (odiag.and.
     &     ((imod.eq.il/2.and.il.lt.jmod.and.jmod.lt.2*il).or.
     &      (jmod.eq.1.5*il)))then
          write(6,*)"iq,imod,jmod=",iq,imod,jmod
          write(6,*)"xlon,rlong,rlat=",xlon,rlong(iq),rlat(iq)
          write(6,*)"jg,glon(jg),jgp",jg,glon(jg),glon(jgp)
          write(6,*)"inc,glat(jg),jgp",inc,glat(jg),glat(jgp)
        endif

! bilinear interp
!       m=int(x)
!       n=int(y)
!       dx=x-m
!       dy=y-n
!       mn=m+(n-1)*ng1
!       mpn=m+1+(n-1)*ng1
!       mnp=m+(n-1+1)*ng1
!       mpnp=m+1+(n-1+1)*ng1
!       d1 = dx*gdat(mpn)+(1-dx)*gdat(mn)
!       d2 = dx*gdat(mpnp)+(1-dx)*gdat(mnp)
!       rdat(iq) = dy*d2+(1-dy)*d1
!       if(iq.eq.20+(70-1)*il)then
!        write(6,*)"dx,dy,rdat(iq)",dx,dy,rdat(iq)
!        write(6,*)"m,n,mn,mpn,mnp,mpnp",m,n,mn,mpn,mnp,mpnp
!        write(6,*)"gdat(mpn),gdat(mn),gdat(mpnp),gdat(mnp)"
!        write(6,*)gdat(mpn),gdat(mn),gdat(mpnp),gdat(mnp)
!       endif
!         ig=x
!         jg=y
!         ig1=ig+(jg-1)*ng1
!         ig2=ig+1+(jg-1)*ng1
!         ig3=ig+(jg-1+1)*ng1
!         ig4=ig+1+(jg-1+1)*ng1
!         write(30,'("ig,jg,gdat=",2i4,5f10.2)') ig,jg,gdat(ig1)
!    &      ,gdat(ig2),gdat(ig3),gdat(ig4),rdat(iq)
!       endif
!       call intp16(globex,ng1+2,ng2,rlong(iq),91.-rlat(iq),
!    &              rdat(iq),.true.)

!***********************************************************************

! prevent extremes greater than global data ???
        rdat(iq)=min(gdx,max(gdn,rdat(iq)))
        dn=min(dn,rdat(iq))
        dx=max(dx,rdat(iq))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      enddo ! iq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

      print *,'sintp16 new grid:dn,dx=',dn,dx
      write(6,'(72("="))')

      return ! subroutine sintp16(gdat,ng1,ng2,rdat,odiag)
      end
!========================================================================
      subroutine intp16(globex,ng1,ng2,x,y,fine,ofirst)

!**********************************************************************
!
! subroutine intp16 interpolates from globex gridpoints to lat-lon
!     location specified by the main program, via a 16- point bessel
!     interpolation scheme.  interpolated
!     lat/lon grid values are contained in the variable "fine"
!     jlm: this is a quadratic fitting pts 2 & 3 & cubic value of centre point
!**********************************************************************

      logical prnt,ofirst
      real globex(ng1,ng2)

      data p25/.25/, c1/1./

! determine the interpolation distances.

      m=x
      mm1=m-1
      mp1=m+1
      mp2=m+2
! this code added 20 Nov 1998
! assumes data has wrap around in e-w direction
      if(mm1.lt.1)mm1=ng1+mm1
      if(mp1.gt.ng1)mp1=mp1-ng1
      if(mp2.gt.ng1)mp2=mp2-ng1

      n=y
      if ( n.lt.1 ) n=1
      if ( n.gt.ng2 ) n=ng2
      nm1=n-1
      np1=n+1
      np2=n+2

      if ( ofirst ) then
        print *,'m,nm1,mp2=',m,mm1,mp2
        print *,'n,nm1,np2=',n,nm1,np2
      endif
      if ( nm1.lt.1 .or. np2.gt.ng2 ) then
       if ( ofirst) then
        print *,'n,m,nm1,np2,ng2=',n,m,nm1,np2,ng2
        print *,'*********** n out of bounds in intp16 ***********'
       endif ! ofirst
       if ( nm1.lt.1 ) nm1=1
       if ( np1.gt.ng2 ) np1=ng2
       if ( np2.gt.ng2 ) np2=ng2
      endif
      dx=x-m
      dy=y-n
      dxx=p25*(dx-c1)
      dyy=p25*(dy-c1)
!     print *,dx,dy,dxx,dyy
!     do j=nm1,np2
!       print *,j,(globex(i,j),i=mm1,mp2,1)
!     enddo

! determine the 16 x-y gridpoints, i.e. top11-top44, to be used
!        for the interpolation.

      top11 = globex(mm1,nm1)
      top12 = globex(mm1,n  )
      top13 = globex(mm1,np1)
      top14 = globex(mm1,np2)
      top21 = globex(m  ,nm1)
      top22 = globex(m  ,n  )
      top23 = globex(m  ,np1)
      top24 = globex(m  ,np2)
      top31 = globex(mp1,nm1)
      top32 = globex(mp1,n  )
      top33 = globex(mp1,np1)
      top34 = globex(mp1,np2)
      top41 = globex(mp2,nm1)
      top42 = globex(mp2,n  )
      top43 = globex(mp2,np1)
      top44 = globex(mp2,np2)

! perform the 16-point interpolation.

      aa = top21 + dx*(top31-top21 + dxx*(top11-top21-top31+top41) )
      ab = top22 + dx*(top32-top22 + dxx*(top12-top22-top32+top42) )
      ac = top23 + dx*(top33-top23 + dxx*(top13-top23-top33+top43) )
      ad = top24 + dx*(top34-top24 + dxx*(top14-top24-top34+top44) )

      fine = ab + dy*(ac-ab + dyy*(aa-ab-ac+ad) )

      !if(ng1.lt.3.and.abs(fine).eq.100.)then
      if(ofirst)then
         write(6,111) m,n,dx,dy,dxx,dyy
111      format (1h ,5x,'(m,n)=',i3,',',i3,3x,'(dx,dy)=',f6.3,',',f6.3,
     .           3x,'(dxx,dyy)=',f6.3,',',f6.3)
         print *,'x,y: ',x,y
         write(6,122) top14,top24,top34,top44,top13,top23,top33,top43,
     .                top12,top22,top32,top42,top11,top21,top31,top41
122      format (1h ,2x,'top11-top44=....',(/,15x,4f12.3) )
         write(6,133) aa,ab,ac,ad,fine
133      format (1h ,2x,'aa=',f10.3,3x,'ab=',f10.3,3x,'ac=',f10.3,3x,
     .             'ad=',f10.3,7x,'fine=',f10.3,/)
      endif

      return
      end
