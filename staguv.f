      subroutine staguv(u,v,uout,vout)   !  unstaguv as an entry
c     N.B. staguv & unstaguv require genuine 2D arrays as input
c     and they need to be contiguous (safest to have in common!)
c     N.B. for 3D input, can use staguv3, unstaguv3 (see entries below)
c     assume k level already given
      include 'newmpar.h'
      include 'parm.h'
      include 'indices.h' ! in,is,iw,ie,inn,iss,iww,iee
      real u(ifull),v(ifull),uout(ifull),vout(ifull)
      real ua(2*ifull),va(ifull),ud(2*ifull),vd(ifull) ! work arrays
      equivalence (ua(ifull+1),va),(ud(ifull+1),vd) ! to ensure contiguous

c     this is staguv
c     unstaggered u & v as input; staggered as output
      do iq=1,ifull
        ua(iq)=.5*(u(ieu2(iq))+u(iq))
        ud(iq)=    u(ieu2(iq))-u(iq)
        va(iq)=.5*(v(inv2(iq))+v(iq))
        vd(iq)=    v(inv2(iq))-v(iq)
      enddo   ! iq loop
      go to 5   ! for shared code

      entry staguv3(u,v,uout,vout)
c     unstaggered u & v as input; staggered as output
c     N.B. staguv3 & unstaguv3 require 3D arrays as input
      do iq=1,ifull
        ua(iq)=.5*(u(ieu(iq))+u(iq))
        ud(iq)=    u(ieu(iq))-u(iq)
        va(iq)=.5*(v(inv(iq))+v(iq))
        vd(iq)=    v(inv(iq))-v(iq)
      enddo   ! iq loop
      go to 5   ! for shared code

      entry unstaguv(u,v,uout,vout)
c     staggered u & v as input; unstaggered as output
c     N.B. staguv & unstaguv require genuine 2D arrays as input
      do iq=1,ifull
        ua(iq)=.5*(u(iq)+u(iwu2(iq)))
        ud(iq)=    u(iq)-u(iwu2(iq))
        va(iq)=.5*(v(iq)+v(isv2(iq)))
        vd(iq)=    v(iq)-v(isv2(iq))
      enddo   ! iq loop
      go to 5   ! for shared code

      entry unstaguv3(u,v,uout,vout)
c     staggered u & v as input; unstaggered as output
c     N.B. staguv3 & unstaguv3 require 3D arrays as input
      do iq=1,ifull
        ua(iq)=.5*(u(iq)+u(iwu(iq)))
        ud(iq)=    u(iq)-u(iwu(iq))
        va(iq)=.5*(v(iq)+v(isv(iq)))
        vd(iq)=    v(iq)-v(isv(iq))
      enddo   ! iq loop

c      now combine the above to give cubic interpolated values
c      code from here on identical in staguv & unstaguv
5      do iq=1,ifull
         uout(iq)=ua(iq)-(ud(ieu2(iq))-ud(iwu2(iq)))/16.
         vout(iq)=va(iq)-(vd(inv2(iq))-vd(isv2(iq)))/16.
       enddo   ! iq loop
c     if(diag)then
c       iq=id+il*(jd-1)
c       print *,'iq,ieu2(iq),iwu2(iq) ',iq,ieu2(iq),iwu2(iq)
c       write (6,'(a,7e20.9)'),'u,ue, ',u(iq),u(ieu2(iq))
c       write (6,'(a,7e20.9)'),'ua, ude, udw '
c    .   ,ua(iq),ud(ieu2(iq)),ud(iwu2(iq))
c     endif
      return
      end
