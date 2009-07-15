      subroutine sintp16(global,ndiag,ng2,fine,gslon,gdx,gslat,gdy,inu)

! NOTE: this is set up for 360,180 grid

      parameter (nx=2,ii0=195,jj0=95,nii=151,njj=111)

!     where nx  = number of intervals to divide 1 degree into (1, 2, or 6),
!           ii0 = integer value of (starting longitude)*nx
!           jj0 = integer value of (starting latitude)*nx
!              N.B. because grib data read in from north to south, must
!              provide ii0, jj0 as coordinates of NW point of subgrid
!           nii = i-dimension of substitution array being read in (0 turns off)
!           njj = j-dimension of substitution array being read in
      include 'gblparm.h'
      parameter (imax=nnx*nx+1,jmax=nny*nx+1)
      include 'newmpar.h'
      include 'xyzinfo.h'  ! x,y,z,rlat,rlong,wts
      include 'vecsuv.h'

c assumes c-c grid dimensions come through 'newmpar.h' as il, jl
c with first data point at 0 deg lon, 90 deg N lat; grib uses LH coord system!!
c southern lats are negative

      !real global(0:359,181),globalex(-1:361,181)
      real global(0:nnx-1,nny),globalex(-1:nnx+1,nny)
      real globex(-1:nnx*nx+1,(nny-1)*nx+1)      ! e.g. (-1:361,181) for nx=1
      real fine(il,jl)

      write(6,'("sintp16x nx,ndiag,ng2,imax,jmax,inu,gdx,gdy=",6i6
     &           ,2f7.2)')nx,ndiag,ng2,imax,jmax,inu,gdx,gdy

      write(6,*) "global ",global(0,1),global(nnx-1,nny)

! extend 1 degree global array over Grenwich meridian
      do j=1,nny
       do i=0,nnx-1
        globalex(i,j)=global(i,j)
       enddo ! i
       globalex(   -1,j)=global(nnx-1,j)
       globalex(  nnx,j)=global(    0,j)
       globalex(nnx+1,j)=global(    1,j)
      enddo ! j

      write(6,*) "globalex ",globalex(-1,1),globalex(nnx+1,nny)

! interpolate into higher resolution (*nx) log-lat array
      if(nx.eq.1)then
        do j=1,nny
         do i=-1,nnx+1
          globex(i,j)=globalex(i,j)
         enddo
        enddo
      else   !  assume larger nx is even at present (in fact 2 or 6)
        nx2=nx/2
!       first copy valid values to corresponding high-res grid points
        do j=1,nny
         jfine=(j-1)*nx+1
         do i=0,nnx
          globex(i*nx,jfine)=globalex(i,j)
         enddo
!        use cubic interp to half-way points in i direction first
         do i=0,359
          globex(i*nx+nx2,jfine)=( 9.*(globalex(i,j)+globalex(i+1,j))
     &          -globalex(i-1,j)-globalex(i+2,j) )/16.
         enddo
         if(nx.eq.6)then
!          linear interp in i direction for 4 remaining intermediate points
           do i=0,359
            globex(i*nx+1,jfine)=
     &          (2.*globex(i*nx,jfine)+globex(i*nx+3,jfine))/3.
            globex(i*nx+2,jfine)=
     &          (globex(i*nx,jfine)+2.*globex(i*nx+3,jfine))/3.
            globex(i*nx+4,jfine)=
     &          (2.*globex(i*nx+3,jfine)+globex(i*nx+6,jfine))/3.
            globex(i*nx+5,jfine)=
     &          (globex(i*nx+3,jfine)+2.*globex(i*nx+6,jfine))/3.
           enddo
         endif     ! (nx.eq.6)
! extend across Grenwich meridian
         globex(-1,jfine)=globex(nnx*nx-1,jfine)
         globex(nnx*nx+1,jfine)=globex(1,jfine)
        enddo  !  j loop
        jj=187
        print *,'A nx & bdy globex ',nx,globex(-1,jj),globex(0,jj),
     &     globex(720,jj),globex(721,jj),

! now interpolate in j direction
! use cubic interp for available half-way points
        do ii=-1,nnx*nx+1
         do j=2,nny-2
          jfine=(j-1)*nx+1
          globex(ii,jfine+nx2)=
     &        ( 9.*(globex(ii,jfine)+globex(ii,jfine+nx))
     &          -globex(ii,jfine-nx)-globex(ii,jfine+2*nx) )/16.
         enddo   ! j loop
! quadratic interp for half-way points near poles
         jfine=1
         globex(ii,jfine+nx2)=.375*globex(ii,jfine)
     &      +.75*globex(ii,jfine+nx) -.125*globex(ii,jfine+2*nx)
         jfine=nny*nx+1
         globex(ii,jfine-nx2)=.375*globex(ii,jfine)
     &      +.75*globex(ii,jfine-nx) -.125*globex(ii,jfine-2*nx)
         if(nx.eq.6)then
! then linear interp in j direction for 4 remaining intermediate points
           do j=1,180
            jfine=(j-1)*nx+1
            globex(ii,jfine+1)=
     &            (2.*globex(ii,jfine)+globex(ii,jfine+3))/3.
            globex(ii,jfine+2)=
     &            (globex(ii,jfine)+2.*globex(ii,jfine+3))/3.
            globex(ii,jfine+4)=
     &            (2.*globex(ii,jfine+3)+globex(ii,jfine+6))/3.
            globex(ii,jfine+5)=
     &            (globex(ii,jfine+3)+2.*globex(ii,jfine+6))/3.
           enddo   ! j loop
         endif     ! (nx.eq.6)
        enddo  ! ii loop
        if(nx.eq.6)then
        print *,'for nx=6 problem ii=950,980 and jj=420,430'
        do ii=950,980
         print 9,ii,(globex(ii,jj),jj=420,430)
9        format(i4,11f7.2)
        enddo
        endif   ! (nx.eq.6)
      endif     !  (nx.eq.1)  .. else ..

      write(6,*) "globex ",globex(-1,1),globex(nnx+1,nny)

      if(nii.gt.0.and.inu.ne.0)then
! now substitute interior grid points written during previous grib run
! calculate fine array indices of input data
        iia=ii0
        iib=iia+nii-1
        jja=90*nx -jj0 +1 ! -ve because jj0=90*nx is North pole (not South pole)
        jjb=jja+njj-1
        if ( inu.gt.0 ) then
          read(inu,err=20,end=20)((globex(ii,jj),ii=iia,iib),jj=jja,jjb)
          go to 30
 20      write(6,*)"setting input data to zero"
         do jj=jja,jjb
         do ii=iia,iib
          globex(ii,jj)=0.
         enddo
         enddo
        else
          read (abs(inu),*) ((globex(ii,jj),ii=iia,iib),jj=jja,jjb)
        endif
 30     dx=-1.e29
        dn= 1.e29
        do jj=jja,jjb
         do ii=iia,iib
          dx=max(dx,globex(ii,jj))
          dn=min(dn,globex(ii,jj))
         enddo
        enddo
        write(6,'("input ",i4," data x/n=",2f12.4)')inu,dx,dn
      endif   ! (nii.gt.0.and.inu.ne.0)then

! loop over all c-c grid points

      write(6,*) " now call intp16",imax,jmax

      do iq=1,ifull
! interpolate from input grid to c-c grid

!      write(6,*) iq,rlong(iq)*nx,(90.-rlat(iq))*nx+1.

       call intp16(globex,imax,jmax,ndiag,ng2,
     &             rlong(iq)*nx, (90.-rlat(iq))*nx +1., fine(iq,1))

       if(iq.lt.3) then
         print *,'iq,rlong,rlat,fine ',iq,rlong(iq),rlat(iq),fine(iq,1)
       endif

      enddo

      return
      end
      subroutine intp16 ( globex,imax,jmax,ndiag, ng2, x, y, fine)

c**********************************************************************
c
c subroutine intp16 interpolates from global gridpoints to lat-lon
c     location specified by the main program, via a 16- point bessel
c     interpolation scheme.  interpolated
c     lat/lon grid values are contained in the variable "fine"
c     jlm: this is a quadratic fitting pts 2 & 3 & cubic value of centre point
c**********************************************************************

      real globex(-1:imax,jmax)   !  (-1:361,181) for nx=1

      data p25/.25/, c1/1./

c determine the interpolation distances.

      m=x
      n=y
      mm1=m-1
      nm1=n-1
      mp1=m+1
      np1=n+1
      mp2=m+2
      np2=n+2
      if ( nm1.lt.0 .or. np2.gt.jmax+1 ) then
        print *,'n,m,nm1,np2,ng2=',n,m,nm1,np2,ng2
        stop '*********** n out of bounds in intp16 ***********'
      endif
      if ( nm1.lt.1 ) nmi=1
      if ( np2.gt.ng2 ) np2=ng2
      dx=x-m
      dy=y-n
      dxx=p25*(dx-c1)
      dyy=p25*(dy-c1)

c     write(6,*) x,y,m,n

c determine the 16 x-y gridpoints, i.e. top11-top44, to be used
c        for the interpolation.

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

c perform the 16-point interpolation.

      aa = top21 + dx*(top31-top21 + dxx*(top11-top21-top31+top41) )
      ab = top22 + dx*(top32-top22 + dxx*(top12-top22-top32+top42) )
      ac = top23 + dx*(top33-top23 + dxx*(top13-top23-top33+top43) )
      ad = top24 + dx*(top34-top24 + dxx*(top14-top24-top34+top44) )
      fine = ab + dy*(ac-ab + dyy*(aa-ab-ac+ad) )

      if(ndiag.lt.3.and.abs(fine).gt.100.)then
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
