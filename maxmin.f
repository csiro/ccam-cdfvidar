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
      
      subroutine maxmin(u,char,ktau,fact,il,kl)

! has more general format & scaling factor
! has entry maxmin0 for single level fields

      !include 'newmpar.h'
      integer il,jl,kl,ifull

      character*2 char
      dimension u(il,6*il,kl)
     & ,umin(kl),umax(kl),iumax(kl),jumax(kl),iumin(kl),jumin(kl)

      jl=6*il
      ifull=il*jl

      kup=kl

2     do k=1,kup
       umax(k)=u(1,1,k)
       umin(k)=u(1,1,k)
       jumax(k)=1
       iumax(k)=1
       jumin(k)=1
       iumin(k)=1
       do j=1,jl
        do i=1,il
         if(umax(k).lt.u(i,j,k))then
           umax(k)=u(i,j,k)
           jumax(k)=j
           iumax(k)=i
         endif
         if(umin(k).gt.u(i,j,k))then
c          print *,'k,i,j,iu,ju,umin,u: ',
c    .              k,i,j,iumin(k),jumin(k),umin(k),u(i,j,k)
           umin(k)=u(i,j,k)
           jumin(k)=j
           iumin(k)=i
         endif
        enddo   ! i loop
       enddo   ! j loop
       umax(k)=fact*umax(k)
       umin(k)=fact*umin(k)
      enddo   ! k loop
      if(umax(kup).gt.30.)then  ! format for T, & usually u,v
        print 971,ktau,char,(umax(k),k=1,kup)
971     format(i7,1x,a2,'max ',10f8.2/(14x,10f8.2)/(14x,10f8.2))
        print 977,ktau,(iumax(k),jumax(k),k=1,kup)
        print 972,ktau,char,(umin(k),k=1,kup)
972     format(i7,1x,a2,'min ',10f8.2/(14x,10f8.2)/(14x,10f8.2))
        print 977,ktau,(iumin(k),jumin(k),k=1,kup)
977     format(i7,'  posij',10(i3,i5)/(14x,10(i3,i5))/(14x,10(i3,i5)))
      else  ! for qg & sd
        print 981,ktau,char,(umax(k),k=1,kup)
981     format(i7,1x,a2,'max ',10f8.3/(14x,10f8.3)/(14x,10f8.3))
        print 977,ktau,(iumax(k),jumax(k),k=1,kup)
        print 982,ktau,char,(umin(k),k=1,kup)
982     format(i7,1x,a2,'min ',10f8.3/(14x,10f8.3)/(14x,10f8.3))
        print 977,ktau,(iumin(k),jumin(k),k=1,kup)
      endif
      return

      entry maxmin0(u,char,ktau,fact,il,kl)
      kup=1
      jl=6*il
      ifull=il*jl      
      go to 2

      end
