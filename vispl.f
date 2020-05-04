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
      
      subroutine vispl ( gridp, grids, prelvs, press, ix, len3
     &                 , sgml, lm, ptop ,nplevs)

      logical olinex
      parameter ( olinex=.false. )

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     this routine interpolates a field 'gridp' on pressure surfaces
c     to 'grids' on sigma surfaces.
c     the interpolation is performed using cubic splines
c     prelvs(nplevs) are the pressure levels and sgml(lm) those
c     of sigma
c     'press' is the surface pressure and 'ix' and 'iy' specify the
c     size of the horizontal grid
c
c     all features of the interpolating spline are determined except
c     for the impostion of one condition at each end
c     these are prescribed through the quantities
c      'cmu1','clmdam','c1' and 'cm'
c
c     for specified slopes at each end - gpd(1) and gpd(nplevs) - have
c      cmu1=clmdam=0.0 , c1=gpd(1) , cm=gpd(nplevs)
c
c     for specified second derivative - gpdd(1) and gpdd(nplevs) -  have
c      cmu1=clmdam=0.5 ,
c     c1=1.5*(gp(2     )-gp(1       ))/h(1       )-
c        h(1       )*gpdd(1     )*0.25
c     cm=1.5*(gp(nplevs)-gp(nplevs-1))/h(nplevs-1)+
c        h(nplevs-1)*gpdd(nplevs)*0.25
c
c     note that the case gpdd(1)=gpdd(nplevs)=0.0 is a particular case
c     of the above and is referred to as the 'natural' spline
c
c     the theory upon which the routine is based may be found,for
c     example , in
c     ahlberg,j.h.,e.n.nilson and j.l.walsh,1967 : the theory of splines
c     and their applications. new york , academic press , 284 pp.
c     (pages 9-15)
c
c**********************************************************************c

      include 'nplevs.h' ! maxplev
      include 'lmax.h'
c
c**********************************************************************c
c
!      common / work / x1(maxplev), x2(maxplev), x3(maxplev)
!     .              , a (maxplev), c (maxplev), d (maxplev)
!     .              , h (maxplev), gd(maxplev), gp(maxplev)
!     .              , pre(lmax), gs(lmax)
      real, dimension(maxplev) :: x1, x2, x3
      real, dimension(maxplev) :: a, c, d, h, gd, gp
      real, dimension(lmax) :: pre, gs
c
      dimension gridp(len3,maxplev), grids(len3,lm), press(len3)
     .        , prelvs(maxplev), sgml(lm)
c
c-----------------------------------------------------------------------
c
      if ( lm.gt.lmax ) then
         print *,'**** lm gt lmax in vispl ****'
         stop 'lmax'
      endif
      if ( nplevs.gt.maxplev ) then
         print *,'**** nplevs gt maxplev in vispl ****'
         stop 'maxplev'
      endif
c
      cmu1=0.5
      clmdam=0.5
      c(1)=cmu1
      a(nplevs)=clmdam
c
c precompute some quantities
c
      nprm1=nplevs-1
c
      do i=1,nprm1
        h(i)=prelvs(i+1)-prelvs(i)
      enddo ! i=1,nprm1
c
      do i=2,nprm1
        x1(i)=0.5/(h(i)+h(i-1))
        x2(i)=h(i)/h(i-1)
        a(i)=h(i)*x1(i)
        c(i)=h(i-1)*x1(i)
      enddo ! i=2,nprm1
c
      c(nplevs)=0.0
c
c**********************************************************************c
c
c begin main horizontal loop
c
      do i=1,ix
c
c**********************************************************************c
c
c sfc pressure
        pres=press(i)
c
c setup press of sigma layers
        do 150 lev=1,lm
  150     pre(lev)=(pres-ptop)*sgml(lev)+ptop
c
        do 40 np=1,nplevs
   40     gp(np)=gridp(i,  np)
c
        c1  =1.5*(gp(2)-gp(1  )) / h(1)
        cm  =1.5*(gp(nplevs)-gp(nplevs-1)) / h(nplevs-1)
        d(1)=c1
        d(nplevs)=cm
c
        do 50 np=2,nprm1
   50     d(np)=3.0*x1(np)*(   (gp(np)-gp(np-1))*x2(np)+
     .                         (gp(np+1)-gp(np))/x2(np)    )
c
c calculate spline derivatives at original points
c
        gd(1)=-c(1)
c
        do 140 np=2,nplevs
          x3(np)=1.0 / (1.0+a(np)*gd(np-1))
          gd(np)=-c(np)*x3(np)
  140     d(np)=(d(np)-a(np)*d(np-1)) * x3(np)
c
        gd(nplevs)=d(nplevs)
c
        do 160 np=1,nprm1
          lev=nplevs-np
  160     gd(lev)=gd(lev)*gd(lev+1)+d(lev)
c
c find sigma levels lying between the limiting pressure levels
        mi1=1
        mim=lm

        do lev=1,lm
          iii=lm+1-lev
          if ( pre(lev).lt.prelvs(1) ) mi1=lev+1
          if ( pre(iii).gt.prelvs(nplevs) ) mim=iii-1
        enddo ! lev=1,lm

c do interpolation
        jj=2
        do ii=mi1,mim
          do 80 ja=jj,nplevs
            if ( pre(ii).le.prelvs(ja) ) go to 90
   80     continue
   90     j1=ja-1
          jj=ja
          v1=prelvs(ja)-pre(ii)
          w1=pre(ii)-prelvs(j1)
          v2=v1*v1
          w2=w1*w1
          h1=h(j1)
          h2=1.0/(h1*h1)
          h3=h2/h1
          z=(gd(j1)*v2*w1-gd(ja)*w2*v1)*h2
          gs(ii)=z+(gp(j1)*v2*(w1+w1+h1)+gp(ja)*w2*(v1+v1+h1))*h3
        enddo ! ii=mi1,mim

        if ( olinex ) then
c linear extrapolation
          if ( mi1.eq.1 ) go to 100
          mlim=mi1-1
          do ii=1,mlim
            gs(ii)=gp(1)+gd(1)*(pre(ii)-prelvs(1))
          enddo ! ii=1,mlim
  100     continue
          if ( mim.eq.lm ) go to 120
          mlim=mim+1
          do ii=mlim,lm
            gs(ii)=gp(nplevs)+gd(nplevs)*(pre(ii)-prelvs(nplevs))
          enddo ! ii=mlim,lm
  120     continue
        else ! ( .not.olinex ) then
c constant extrapolation
          if ( mi1.ne.1 ) then
            mlim=mi1-1
            do ii=1,mlim
              gs(ii)=gp(1)
            enddo ! ii=1,mlim
          endif ! ( mi1.eq.1 ) then
          if ( mim.ne.lm ) then
            mlim=mim+1
            do ii=mlim,lm
              gs(ii)=gp(nplevs)
            enddo ! ii=mlim,lm
          endif ! ( mim.eq.lm ) then
        endif ! ( .not.olinex ) then
c stuff values into output array
        do ii=1,lm
          grids(i,ii)=gs(ii)
        enddo ! ii=1,lm
c**********************************************************************c
c end of ix loop
      enddo ! i=1,ix
c**********************************************************************c
      return
      end
