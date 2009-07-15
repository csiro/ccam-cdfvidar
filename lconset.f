      subroutine lconset(ds)
      common/mapproj/du,tanl,rnml,stl1,stl2
      common/lconther/ther
c rotate added to do S. Pole properly Wed  11-30-1994
c   this is because usually increasing j goes northwards, but not so for
c   S. Pole polar stereographic
c rewritten completely by jlm August 1993 & May 1995 for SH & NH  (lamb.for)
c it is cleaned up, & expects -ve latitudes for S. Hem.
c tanl with proper sign is latitude where j=1
c du is i value along standard longitude rnml
c all colatitudes from pi/2, so need sign=-1. for SH

c this subroutine converts the grid coordinates ri,rj of a point to its
c latitude (rlat) and longitude (rlong) or vice versa. 
c  Lambert conformal projection is used.
c     stl1,stl2  standard latitudes   -ve for S. Hem.
c     ds         length of grid unit(m) at standard latitudes
c     rearth     radius of earth(m)

      data pi/3.1415926536/,rearth/6371.22e3/,cfact/1.45444e-4/,sign/1./
      save
c     e.g. gu=137.057, du=19.76, tanl=-60.2, stl1=-40., stl2=-10., rnml=130.
cjlm  use chi for colatitudes in radians
      print *,'lconset ds,du,tanl,rnml,stl1,stl2'
      print *,ds,du,tanl,rnml,stl1,stl2
      pole=.5*pi
      if(.5*(stl1+stl2).lt.0.)sign=-1.       ! decides which hemisphere
      rotate=1.
      if(stl1.lt.0..and.stl1.eq.stl2)rotate=-1.  ! for grid over S. pole
      chi1=pole-sign*stl1*pi/180.
      chi2=pole-sign*stl2*pi/180.
      if(stl1.eq.stl2)then
        rk=1.
        rmult=(1.+cos(chi1))  *rearth/ds
      else
        rk=log(sin(chi1)/sin(chi2))/log(tan(.5*chi1)/tan(.5*chi2))    ! i.e. K
        rmult=sin(chi1)*tan(.5*chi1)**(-rk)  *rearth/(ds*rk)
      endif
      chiref2=.5*(pole-sign*tanl*pi/180.)
      refj=1.+rotate*sign*rmult*tan(chiref2)**rk
c     print *,'chi1,chi2,chiref2,rmult,refj,rk '
c    .        ,chi1,chi2,chiref2,rmult,refj,rk
      return

      entry lconij(rlong,rlat,ri,rj)
      theta=rk*(rlong-rnml)
      ther=theta*pi/180.
      chion2=.5*(pole-sign*rlat*pi/180.)
c     print *,'chion2,tan(chion2) ',chion2,tan(chion2)
      r=rmult*tan(chion2)**rk
      rj=refj-sign*r*cos(ther)*rotate   ! Wed  11-30-1994
c     rj=refj-sign*r*cos(ther)
      ri=du+r*sin(ther)*rotate          ! Wed  11-30-1994
      return

      entry lconll(rlong,rlat,ri,rj)
      ther=atan2( ri-du,rotate*sign*(refj-rj) )  ! Sun  05-07-1995
      rlong=rnml+rotate*ther*180./(pi*rk)
      if(rlong.gt.360.)rlong=rlong-360.
c     print *,'rmult,rk,ther,refj-rj ',rmult,rk,ther,refj-rj
c     print *,'cos(ther),rlong ',cos(ther),rlong
c     print *,'ratio ',(rj-refj)/(rmult*cos(ther))
      chion2= sign*atan( (-rotate*sign*(rj-refj)/
     .        (rmult*cos(ther)))**(1./rk) ) !OK 1,3
      rlat=(sign*pole-2.*chion2)*180./pi  ! Sun  05-07-1995
      return

      entry mapff ( rlat, facmap, coriol )
      chi=pole-sign*rlat*pi/180.
c     print *,'rk = ',rk
      if(rk.eq.1.)then
        facmap=(1.+cos(chi1))/(1.+cos(chi))  ! i.e. for polar stereographic
      else
        term1=sin(chi1)/sin(chi)
        term2=(tan(.5*chi)/tan(.5*chi1))**rk
        facmap=term1*term2
      endif
      coriol = sign*cfact*cos(chi)
      return
      end
c       include 'maxmin.f'
c       include 'setxyz.f'
c       include 'staguv.f'
c       include 'jimcc.f'
