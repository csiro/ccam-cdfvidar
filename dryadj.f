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
      
      subroutine dryadj ( t, pd, pt, len, icnt, len3, dsg, sgml, lm )
c
c**********************************************************************c
c                                                                      c
c ********** dry convective adjustment **********                      c
c                                                                      c
c input: dsg(l), pd(i), sgml(l), pt, t(i,l), rn                        c
c                                                                      c
c***********************************************************************
c                                                                      c
c dry adiabatic adjustment of unstable atmospheric layers in the model c
c               original non-vectorized version                        c
c                                                                      c
c**********************************************************************c
c
      include 'lmax.h'
c
      dimension dsg(lm), sgml(lm)
      dimension papai(lmax), rdsd(lmax)
      dimension t(len3,lm), pd(len)
c
      logical drca
c
      data rn/3.e-3/
c
c***********************************************************************
c
      if ( lm.gt.lmax ) then
         print *,'**** lm gt lmax in dryadj ****'
         stop 'lmax'
      endif
c
      lmm1=lm-1
      do 100 l=1,lmm1
        rdsd(l)=dsg(l)/dsg(l+1)
 100  continue
c
c.......................................................................
c
c main horizontal loop
c
      do 230 i=1,len
c
c.......................................................................
c
          pdij=pd(i)
c
          do 451 l=1,lm
 451        papai(l)=(pt+sgml(l)*pdij)**(-.2858964143)
c
 232      drca=.false.
          pplp1=papai(1)
c
          do 231 l=1,lmm1
c
            ppl=pplp1
            pplp1=papai(l+1)
            thl=t(i,l)*ppl
            thlp1=t(i,l+1)*pplp1
            dth=thlp1-thl
c
            if ( dth.gt.rn ) then
c
               dtl=dth/(ppl+pplp1*rdsd(l))
               dtlp1=-dtl*rdsd(l)
               t(i,l)=t(i,l)+dtl
               t(i,l+1)=t(i,l+1)+dtlp1
               drca=.true.
               icnt=icnt+1
c
            endif
c
 231      continue
c
          if (drca) go to 232
c
c.......................................................................
c
c end of horizontal loop
c
 230  continue
c
c.......................................................................
c
      return
      end
