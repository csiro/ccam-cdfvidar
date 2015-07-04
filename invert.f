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
      
      subroutine invert3(temp,il,kl)

c invert vertical levels in temp() like darlam (bottom up)

      !include 'newmpar.h' ! il,jl,kl,ifull,ijk

      integer il,jl,ifull,kl

      real temp(6*il*il,kl)
      real tmp(6*il*il,kl)

      jl=6*il
      ifull=il*jl

      do k=1,kl
       do iq=1,ifull
        tmp(iq,k)=temp(iq,kl+1-k)
       enddo ! n=1,ifull
      enddo ! k

      do k=1,kl
       do iq=1,ifull
        temp(iq,k)=tmp(iq,k)
       enddo ! n=1,ifull
      enddo ! k

      return ! invert3
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine invert1(sig,kl)

      !include 'newmpar.h' ! il,jl,kl,ifull,ijk

      integer kl

      real sig(kl)
      real tmp(kl)

      do k=1,kl
        tmp(k)=sig(kl+1-k)
!       write(6,*)"k=",kl+1-k," sig=",sig(kl+1-k)
      enddo ! k

      do k=1,kl
        sig(k)=tmp(k)
!       write(6,*)"k=",k," sig=",sig(k)
      enddo ! k
      return ! invert1
      end
