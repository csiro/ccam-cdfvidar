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

      logical spline,oesig,debug,notop,opre,calout,oform,have_gp

      logical splineu,splinev,splinet,zerowinds,osig_in

      character*1024 zsfil,tsfil,smfil,vfil

      common / vi / ntimes,spline,mxcyc,nvsig,nrh,                       &
     &              oesig,ptop,debug,notop,opre,have_gp,                 &
     &              in,calout,                                           &
     &              iout,oform,osig_in,                                  &
     &              inzs,zsfil,ints,tsfil,insm,smfil,                    &
     &              vfil,                                                &
     &              splineu,splinev,splinet,zerowinds
