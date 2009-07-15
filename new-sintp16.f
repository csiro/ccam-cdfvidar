      subroutine sintp16(global,ndiag,ng2,fine,globex)

      include 'newmpar.h'
      include 'xyzinfo.h'  ! x,y,z,wts
      include 'latlong.h'  ! rlat,rlong
      include 'vecsuv.h'

c assumes fine grid dimensions come through 'newmpar.h'
c as il, jl
c also assumes that you have already called lconset(ds)
c needs lconset routine.

c the global grid is dimensioned ng1 by ng2
c with first data point at 0 deg lon, 90 deg N lat   ******  LH coord system!!

c southern lats are negative

      real global(ndiag,ng2)
      !real globex(-1:363*ng2)
      real globex(-1:(ndiag+3)*ng2)
      real fine(il,jl)

      ng1p=ndiag+3

c     extend global array over Grenwich meridion
      dx=-1.e29
      dn= 1.e29
      do j=1,ng2
       do ii=1,ndiag
        iq=ii-1+(j-1)*ng1p
        globex(iq)=global(ii,j)
        dn=min(dn,global(ii,j))
        dx=max(dx,global(ii,j))
       enddo
       globex(-1+(j-1)*ng1p)=global(ndiag,j)
       globex(ndiag+(j-1)*ng1p)=global(1,j)
       globex(ndiag+1+(j-1)*ng1p)=global(2,j)
      enddo
      print *,'sintp16 dn,dx=',dn,dx

c loop over all lam grid points
      nngx = ndiag+1
      write(6,*)"nngx,ng2=",nngx,ng2
      xx = rlong(1)/dlon
      yy = 90-rlat(1)/dlat+1.
      write(6,*)"lon,lat(1),xx,yy=",rlong(1),rlat(1),xx,yy
      xx = rlong(ifull)/dlon
      yy = 90-rlat(ifull)/dlat+1.
      write(6,*)"lon,lat(ifull),xx,yy=",rlong(ifull),rlat(ifull),xx,yy
      xx = 360./dlon
      yy = 180/dlat+1.
      write(6,*)"lon,lat(ifull),xx,yy=",rlong(ifull),rlat(ifull),xx,yy

      dn=1.e29
      dx=-1.e29
      do iq=1,ifull

        xx = 1.+(rlong(iq)-slon)/dlon
        yy = 1.+(rlat(iq)-slat)/dlat
        if(ndiag.gt.300.)then
          xx = 1.+(rlong(iq)-slon)
          yy = (90.-rlat(iq))+1.
        else
          xx = rlong(iq)/1.25
          yy = (90.-rlat(iq))/1.25+1.
        endif
        !write(6,*)iq,rlong(iq),rlat(iq),xx,yy,ndiag

c interpolate from input grid to lam grid

        call intp16(globex,nngx,ng2,xx,yy,fine(iq,1),.true.)

        dn=min(dn,fine(iq,1))
        dx=max(dx,fine(iq,1))

c       fine(i,j)=max(xmin,fine(i,j))

c       print *,fine(i,j)
!     if(ndiag.lt.3.and.abs(fine(iq,1)).gt.100.)then
!       print *,'iq,rlong,rlat ',iq,rlong(iq),rlat(iq)
!     endif

      enddo ! iq
      print *,'sintp16 dn,dx=',dn,dx

      return
      end
      subroutine intp16 ( globex, ndiag, ng2, x, y, fine, ofirst )

c**********************************************************************
c
c subroutine intp16 interpolates from global gridpoints to lat-lon
c     location specified by the main program, via a 16- point bessel
c     interpolation scheme.  interpolated
c     lat/lon grid values are contained in the variable "fine"
c     jlm: this is a quadratic fitting pts 2 & 3 & cubic value of centre point
c**********************************************************************

      logical prnt,ofirst
      real globex(-1:ndiag,ng2)

      data prnt/.false./
      data p25/.25/, c1/1./

c determine the interpolation distances.

      m=x
      mm1=m-1
      mp1=m+1
      mp2=m+2
      n=y
      nm1=n-1
      np1=n+1
      np2=n+2
      if ( nm1.lt.1 .or. np2.gt.ng2 ) then
        if ( ofirst) then
          print *,'m,n,nm1,np2,ng2=',n,m,nm1,np2,ng2
          print *,'*********** n out of bounds in intp16 ***********'
        endif ! ofirst
        if ( nm1.lt.  1 ) nm1=1
        if ( np1.gt.ng2 ) np1=ng2
        if ( np2.gt.ng2 ) np2=ng2
      endif
      dx=x-m
      dy=y-n
      dxx=p25*(dx-c1)
      dyy=p25*(dy-c1)

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
