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
      
      subroutine vidar(nplevs,zp,tp,up,vp,hp,validlevcc
     &                ,iyr,imon,idy,ihr,nt,time,mtimer,pm,io_out,il,kl)
     
      use comsig_m, dsg => dsgx, sgml => sgmlx, sg => sgx
      use sigdata_m
! version for globpe
c**********************************************************************
c   interpolation from nplevs pressure levels to lm sigma levels       
c     for imfjmf grid, compute pstar, check temps for dry adiabat      
c     lapse rates & adjust                                             
c**********************************************************************
c  zs (in sigdata.h) : output grid zs (m)
c  ints : disc with t*
c  insm : disc with soilm
c      
c   spline : turns on cubic spline vertical interpolation (for u,v,t
c            only)  r,h are always linearly interpolated

c   nt : time period to process
c   ptop : pressure at top sigma level (0.0 mb)

c***********************************************************************
      !include 'newmpar.h'

      integer il,jl,kl,ifull
      integer imf,jmf,lm,im,imp1,imid
      !parameter ( imf=il, jmf=jl, lm=kl )

      include 'nplevs.h' ! maxplev

      !parameter ( im=imf )

      !parameter ( lmp1 = lm+1 )
      !parameter ( imid = il/2+il*(jl/2-1) )

      !common / comsig / dsg(kl), sgml(kl), sg(kl+1)

c**********************************************************************

      include 'vidar.h'

!        logical spline,oesig,debug,notop,opre,calout,oform,have_gp
!        logical splineu,splinev,splinet,zerowinds,osig_in
!        character*80 zsfil,tsfil,smfil,vfil
!        common / vi / ntimes,spline,mxcyc,nvsig,nrh
!       &             ,oesig,ptop,debug,notop,opre,have_gp
!       &             ,in,calout
!       &             ,iout,oform,osig_in
!       &             ,inzs,zsfil,ints,tsfil,insm,smfil
!       &             ,vfil
!       &             ,splineu,splinev,splinet,zerowinds
c**********************************************************************

      real zp(6*il*il,maxplev),tp(6*il*il,maxplev),rp(6*il*il,maxplev)
     &    ,hp(6*il*il,maxplev),up(6*il*il,maxplev),vp(6*il*il,maxplev)
      real validlevcc(6*il*il)

      !include 'sigdata.h'
!n    common/sigdata/pmsl(ifull),sfct(ifull),zs(ifull),ps(ifull)
!n   &             ,us(ifull,kl)    ,vs(ifull,kl)    ,ts(ifull,kl)
!n   &             ,rs(ifull,kl)    ,hs(ifull,kl)    ,psg_m(ifull)
!n   &             ,zsi_m(ifull)
!o    common/sigdata/pmsl(ifull),sfct(ifull),zs(ifull),ps(ifull)
!o   &             ,us(ifull,lm),vs(ifull,lm),ts(ifull,lm)
!o   &             ,rs(ifull,lm),hs(ifull,lm),psg_m(ifull)
!o   &             ,zsi_m(ifull)

      common/datatype/moist_var,in_type ! set in cdfvidar
      character*1 in_type
      character*2 moist_var
      character*10 header

      real alm(6*il*il),prein(6*il*il)

      common / labcom / lab(17)
                        character*4 lab

      common / mapproj / du,tanl,rnml,stl1,stl2

      real pm(maxplev), alpm(maxplev)
      real ac(kl+1), bc(kl+1), cc(kl+1), dc(kl+1)
      real tpold(6*il*il,maxplev) ! MJT suggestion

      integer pcoun
      integer dypmo(12)

c**********************************************************************

      parameter (epsiln=0.05, rrr=287.04, grav=9.80665, ccp=1.00464e3)
      parameter (c0=0., c100=100.)

      data dypmo/31,28,31,30,31,30,31,31,30,31,30,31/

c***********************************************************************
      save

      jl=6*il
      ifull=il*jl
      imf=il
      jmf=jl
      lm=kl
      im=imf
      lmp1 = lm+1
      imid = il/2+il*(jl/2-1)

      write(6,*)"#####################################################"

      print *,'vidar *** nplevs=',nplevs," nt=",nt," time=",time

c initialize some constants
      grds=grav/rrr
      cuppa=rrr/ccp
      lmm1 = lm-1
      nplevsm=nplevs-1

      imt=il   !  jjk's was +1
      jmt=jl   !  jjk's was +1
      lmt=kl
      npts=imt*jmt
      if ( npts .gt. ifull ) stop 'ifull'

c pressure coordinates (pa)
c top down

       write(6,*)"convert pm from hPa to Pa  osig_in=",osig_in
       write(6,*)"then precompute alog(pm)"

       if (osig_in) then
         idiag = il/2+(jl/2-1)*il
         do k=1,nplevs
            pm(k)=pm(k)*1.e2 ! pm now Pa
            pr = calc_p(psg_m(idiag),pm(k)) ! psg_m(Pa), pm(sig*1.e5)-now
            alpm(k)=alog(pr)
         end do ! m=1,nplevs
       else
         do k=1,nplevs
           pm(k)=pm(k)*1.e2 ! pm now Pa
           alpm(k)=alog(pm(k))
         end do ! m=1,nplevs
       end if

c print out values determined in sigma

           write(6,9010) (sgml(l),l=1,lm)
 9010      format(1x,'sgml  ',9e12.4/(7x,9e12.4))
           write(6,9030) (sg(l),l=1,lmp1)
 9030      format(1x,'sg    ',10e12.4/(7x,10e12.4))
           write(6,9040) (pm(l),l=1,nplevs)
 9040      format(1x,'pm(Pa)',12e10.3/(7x,12e10.3))
           write(6,9050) (alpm(l),l=1,nplevs)
 9050      format(1x,'alpm  ',12e10.3/(7x,12e10.3))

c determine some additional values

           do n = 1, lmm1
             dc(n)=cuppa*(sgml(n+1)-sgml(n))/(sgml(n+1)+sgml(n))
             ac(n)=sg(n+2)-sg(n)
             bc(n)=sg(n+2)-sg(n+1)
             cc(n)=sg(n+1)-sg(n)
           end do ! n = 1, lmm1

!########## if nt=1 #####################
      if ( nt.eq.1 ) then
!########## if nt=1 #####################

          print *,'**** imt,jmt,lmt,npts= ',imt,jmt,lmt,npts,' *****'
          idiag = il/2+(jk+jl/2-1)*il

!########## end nt=1 #####################
      endif ! nt.eq.1
!########## end nt=1 #####################

      if ( ints.gt.0 .and. ints.ne.inzs ) then
          print *,'open ',ints,' tsfil=',tsfil
          open (ints,file=tsfil,form='unformatted',status='unknown')
      elseif ( ints.lt.0 .and. nt.eq.1 ) then
          print *,'open ',ints,' tsfil=',tsfil
          open(abs(ints),file=tsfil,form='formatted',status='unknown')
      endif

      if ( insm.gt.0 .and. insm.ne.inzs ) then
          print *,'open ',insm,' smfil=',smfil
          open (insm,file=smfil,form='unformatted',status='unknown')
      endif

      print *,'ints,insm=',ints,insm
!     if ( insm.gt.0 ) rewind insm

c***********************************************************************
      write(6,5050) nt
 5050 format(1x,'*** vidar starts *** nt=',i4)
c***********************************************************************
      write(6,*)"input iyr=",iyr," imon=",imon
     &              ," idy=",idy," ihr=",ihr
 51   if ( ihr.gt.24 ) then
        ihr=ihr-24
        idy=idy+1
        if ( mod(iyr,4).eq.0 ) dypmo(2)=29
        if ( idy.gt.dypmo(imon) ) then
          idy=idy-dypmo(imon)
          imon=imon+1
          if ( imon.gt.12 ) then
            imon=imon-12
            iyr=iyr+1
          endif ! ( imon.gt.12 ) then
        endif ! idy.gt.dypmo(imon) ) then
        go to 51
      endif
      write(6,*)"fixed input iyr=",iyr," imon=",imon
     &     ," idy=",idy," ihr=",ihr
!     if ( nt.eq.1 .and. calout ) then
!        ih1=ihr/10
!        ih2=ihr-10*ih1
!        id1=idy/10
!        id2=idy-10*id1
!        im1=imon/10
!        im2=imon-10*im1
!        write(vfil,'("a",i2,6i1)') iyr,im1,im2,id1,id2,ih1,ih2
!     endif ! ( nt.eq.1 ) then

      write(6,*)"pmsl (1)=",pmsl(1)," npts=",pmsl(npts)
      call prt_pan(pmsl,il,jl,2,'pmsl(Pa)')

c***********************************************************************

! hp comes from cdfvidar.  Can be either rh or mixrat.

      if ( moist_var .eq. "rh" ) then

        write(6,*)"input rh (0-100), output mix.rat(g/g) with kl=",kl

        do k=1,nplevs

! rp is returned with h2o saturation partial pressures
          call esmtrv ( tp(:,k), rp(:,k), npts, rs(:,1),
     &                  ts(:,1), us(:,1), vs(:,1) )

          pr = .01*pm(k) ! hPa
          if ( osig_in ) pr = .01*calc_p(psg_m(imid),pm(k)) ! hPa

          write(6,'("k,p(hPa),tp,h,esp=",i3,3f10.2,1p,e12.4)')
     &               k,pr,tp(imid,k),hp(imid,k),rp(imid,k)

! calculate actual.mix.ratio (input hp = rh, rp = sat vp )
          do iq=1, npts
            hp(iq,k)=max(min(hp(iq,k),100.),0.) ! keep 0 < rh < 100
            ! why was this here ????? hp(iq,2)=hp(iq,1) ! removed 07 dec 2005
            pr = pm(k) ! Pa
            if ( osig_in ) pr = calc_p(psg_m(iq),pm(k)) ! Pa
            rp(iq,k)=0.622*rp(iq,k)/(pr-rp(iq,k)) ! sat. mix.rat. from es
            rp(iq,k)=hp(iq,k)*rp(iq,k)*.01 ! actual mix.rat. = rh*qsat
          enddo ! iq=1, npts

          pr = .01*pm(k) ! hPa
          if ( osig_in ) pr = .01*calc_p(psg_m(imid),pm(k)) ! hPa
          write(6,'("k,p,tp,hp,rp=",i3,3f10.2,1p,e12.4)')
     &               k,pr,tp(imid,k),hp(imid,k),rp(imid,k)

        enddo ! k=1,nplevs

      else

        write(6,*)"input mix.rat(g/g), output rh (0-100) with kl=",kl

        do k=1,nplevs

! hp is returned with h2o saturation partial pressures
          call esmtrv ( tp(:,k), rp(:,k), npts, rs(:,1), ts(:,1),
     &                  us(:,1), vs(:,1) )

          pr = .01*pm(k)
          if ( osig_in ) pr = .01*calc_p(psg_m(imid),pm(k))
          write(6,'("k,p,tp,h,esp=",i3,3f10.2,1p,e12.4)')
     &               k,pr,tp(imid,k),hp(imid,k),rp(imid,k)

! calculate rel. hum. ( input hp = mixrat, rp = sat vp )
          do iq=1, npts
            pr = pm(k)
            if ( osig_in ) pr = calc_p(psg_m(iq),pm(k))
            q = hp(iq,k)
            satmr=0.622*rp(iq,k)/(pr-rp(iq,k)) ! sat. mix.rat. from es
            if ( satmr .gt. 1.e-10 ) then
              hp(iq,k)=100.*q/satmr              ! rh(%) = 100 * actual mix.rat. / qsat
            else
              hp(iq,k)=0.
            endif
            hp(iq,k)=max(min(hp(iq,k),100.),0.) ! keep 0 < rh < 100
            rp(iq,k)=q
          enddo ! iq=1, npts

          pr = .01*pm(k)
          if ( osig_in ) pr = .01*calc_p(psg_m(imid),pm(k))
          write(6,'("k,p,tp,hp,rp=",i3,3f10.2,1p,e12.4)')
     &               k,pr,tp(imid,k),hp(imid,k),rp(imid,k)

        enddo ! k=1,nplevs

      endif ! ( moist_var .eq. "rh" ) then
c***********************************************************************

      write(6,*)"inzs=",inzs

      if ( inzs.eq.0 ) then
c use input zs if inzs = 0
c determine land mask
         do i=1,npts
           if ( zs(i).gt.0.01 ) then
              alm(i)=1.
           else
              alm(i)=0.  ! here assuming 0 zs = water!
           endif
         enddo ! i=1,npts

      else ! inzs.ne.0

!J      write(6,*)"read in new lm;not using zs here,using input grid zs"
!J      rewind(inzs)

!J       read(inzs,*)ilx,jlx,rlong0,rlat0,schmidt,ds,header
!J       write(6,*)ilx,jlx,rlong0,rlat0,schmidt,ds,header
!J       write(6,*)"npts=",npts,il,jl

!J       if ( npts.ne.(ilx*jlx) ) stop 'zs has wrong # of points'

! dummy read to get past zs
!J       read(inzs,*)(alm(i),i=1,npts) ! added alm to read 26-apr-2001

!        read(inzs,*)(zs(iq),iq=1,ifull)  ! formatted zs
         zsmax=-1.e29
         zsmin=+1.e29
         do i=1,npts
           zsmax=max(zsmax,zs(i))
           zsmin=min(zsmin,zs(i))
         enddo !i=1,npts
         write(6,'(''zs(m) npts: max,min='',2f15.2)')zsmax,zsmin

! now read land-sea mask
!J       read(inzs,*)(alm(i),i=1,npts)
!J       if(nt.eq.1)print *,'alm read 1,',npts,'=',alm(1),alm(npts)

      endif ! inzs.eq.0

      if ( notop ) then
         write(6,*)"############ reset topog if notop=t notop=",notop
         do i=1,npts
           zs(i)=0.0
         enddo !i=1,npts
      endif

      if ( nt.eq.1 ) then
         !!call amap ( zs, imt, jmt, 'zs(m)', 100., 0. )
         !!call amap ( alm, imt, jmt, 'lmask', .5, 0. )
         call prt_pan(zs,il,jl,2,'zs(m)')
         call prt_pan(alm,il,jl,2,'zs(m)')
      endif

c***********************************************************************
c begin processing data
c  ltcoun = counts total nr of grid columns where lapse rate
c           exceeded dry adiab on first pass
c  pcoun  = counts total number of grid points which will not converge
      ltcoun = 0
      pcoun  = 0
c***********************************************************************
c print out data at i=1 and j=1,nplevs
      if ( debug ) then
         write(6,'(1x,"pressure-level data 1st and last point:")')
         write(6,98)(k,up(1,k),vp(1,k),tp(1,k),rp(1,k),hp(1,k)
     &            ,zp(1,k), k=1,nplevs)
         write(6,98)(k,up(npts,k),vp(npts,k),tp(npts,k),rp(npts,k)
     &            ,hp(npts,k),zp(npts,k), k=1,nplevs)
 98      format(10x,'up',11x,'vp',11x,'tp',11x,'rp',11x,'hp',11x,'zp'
     &          /(i3,1x,6e13.5))
      endif
c***********************************************************************
c convert sensible temp to virt. temp
      write(6,'(a3,6a12)') 
     &          "k","pm(hPa)","zp(m)","tp(K) (1)","(npts) "
     &                           ," rp (g/g) (1)","(npts)"
      tpold=tp ! MJT suggestion
      do k=1,nplevs
        do i=1,npts
          tp(i,k)=tp(i,k)*(rp(i,k)+.622)/(.622*(1.+rp(i,k)))
        enddo ! i=1,npts
        pr = .01*pm(k)  ! hPa
        if(osig_in)pr  = .01*calc_p(psg_m(1),pm(k)) ! hPa
        write(6,'(i3,4f12.4,1p,2e12.4)')
     &         k,pr,zp(1,k),tp(1,k),tp(npts,k),rp(1,k),rp(npts,k)
      enddo ! k=1,nplevs
c***********************************************************************
c start interpolating from pressure to sigma
      write(6,*)"=====> determine p* (pressure at surface)"
      write(6,*)"************************************> have_gp=",have_gp
      g=grav
c***********************************************************************
      if ( have_gp ) then
c***********************************************************************
       do i=1, npts          
c loop through pressure levels ( bottom -> up )
c loop through pressure levels ( top-down? )
        do kk=1,nplevsm

!         if ( i.eq.idiag .or. zs(i).gt.2000. ) then
!           write(6,*)zs(i),zp(i,kk+1),kk,nplevsm,pm(kk)
!         endif ! ( i.eq.idiag .or. zs(i).gt.2000. ) then

          if(zs(i).lt.zp(i,kk+1).and.kk.ne.nplevsm ) go to 172

          tem1=tp(i,kk)
          tem2=tp(i,kk+1)

          if ( osig_in ) then
            alp  = alog(calc_p(psg_m(i),pm(kk  )))
            alpp = alog(calc_p(psg_m(i),pm(kk+1)))
            if ( i .eq. idiag ) then
              write(6,*)"i,kk,psg_m,pm,calcp="
              write(6,*)i,kk,psg_m(i),pm(kk),calc_p(psg_m(i),pm(kk))
            endif ! ( i .eq. idiag ) then
          else
            alp  = alpm(kk)
            alpp = alpm(kk+1)
          endif ! ( osig_in ) then

          icase=0
          if ( zs(i) .lt. zp(i,nplevs) ) then
            !!write(6,*)"NEW CASE setting tem1=tem2=tp(i,nplevs)"
            tem1=tp(i,nplevs)
            tem2=tp(i,nplevs)
            icase=-1
          endif

          sem1= .5*(tem2-tem1)/(alpp-alp)
          sem3= grds*(zs(i)-zp(i,kk+1))

!         if ( i.eq.idiag .or. zs(i).gt.2000. ) then
          if ( i.eq.idiag ) then
            write(6,*)"i,kk,tem1,tem2,alp,alpp,pm(kk)="
            write(6,'(2i6,5f10.3)')i,kk,tem1,tem2,alp,alpp,pm(kk)
            write(6,*)"zs,zp(kk+1)=",zs(i),zp(i,kk+1)
            write(6,*)"sem1,sem3=",sem1,sem3
          endif ! ( i .eq. idiag ) then

          if ( sem1.ne.c0 ) then
c non-isothermal case (general)
             rem1 = tem2**2 - 4.*sem1*sem3
             if ( rem1.lt..00001 ) then
                write(6,142)i,kk,rem1,tem1,tem2,zs(i),zp(i,kk+1)
  142           format(6x,'** diagnostic in ps calc: rem1.le..00001 '
     &          ,' i=',i5,' kk=',i3,' rem1='
     &          ,e12.4,/,10x,'t(kk),t(kk+1),zs,z(kk+1)=',4f10.1)
                write(6,*)"i,kk,tem1,tem2,alp,alpp="
                write(6,*)i,kk,tem1,tem2,alp,alpp
                write(6,*)"zs,zp(kk+1)=",zs(i),zp(i,kk+1)
                write(6,*)"sem1,sem3=",sem1,sem3
                rem1=.001
                rem1=.0
                icase=10+icase
             endif
             rem2=2.0*sem3/(tem2+sqrt(rem1))
          else
c isothermal case
             rem2=sem3/tem2
             icase=20+icase
          endif

c calculate surface pressure
          if ( osig_in ) then
            pr=calc_p(psg_m(i),pm(kk+1))
          else
            pr=pm(kk+1)    
          end if
          
          rem2 = max( min( rem2, 60. ), -60. ) ! MJT suggestion for single precision

          ps(i)=pr*exp(-rem2)

!         if ( i.eq.idiag .or. zs(i).gt.2000. ) then
          if ( i.eq.idiag  ) then
            write(6,*)"i,kk,ps,rem1,rem2=",i,kk,ps(i),rem1,rem2,icase
          endif

          go to 173
c end of pressure loop
 172      continue
        enddo ! kk=1,nplevsm

        do kk=1,nplevsm ! (bottom - up )
          if ( zs(i).lt.zp(i,kk) ) then ! zs below this level
            if ( kk .eq. 1 ) then
              if ( osig_in ) then
                pr1=calc_p(psg_m(i),pm(1))
                pr2=calc_p(psg_m(i),pm(2))
              else
                pr1=pm(1)
                pr2=pm(2)
              endif ! ( osig_in ) then
              ps(i)=pr
     &             +(zs(i)-zp(i,1))*(pr2-pr1)
     &                             /(zp(i,2)-zp(i,1))
            else
              if ( osig_in ) then
                pr =calc_p(psg_m(i),pm(kk  ))
                prm=calc_p(psg_m(i),pm(kk-1))
              else
                pr =pm(kk)
                prm=pm(kk-1)
              endif ! ( osig_in ) then
              ps(i)=prm
     &             +(zs(i)-zp(i,kk-1))*(pr-prm)
     &                                /(zp(i,kk)-zp(i,kk-1))
            endif
            go to 173
          endif
        enddo ! kk
 173    continue
       enddo ! i        
      else ! if ( ! have_gp ) then

       do i=1, npts
         k=nplevs+1-nint(validlevcc(i))
         avgtmp=tp(i,k)
         
         ! Use lowest pressure level temp as to compute ps from pmsl
         !ps(i)=exp(log(100.*pmsl(i))-max(0.,grav*zs(i))/(rrr*avgtmp)) ! jjk
 	   
         tmsl=avgtmp +zs(i)*.0065                                     ! jlm
         ps(i)=100.*pmsl(i)*(1.-.0065*zs(i)/tmsl)**(grav/(.0065*rrr)) ! jlm

         if ( mod(i,100).eq.0 ) then
           write(6,*)"pmsl(i),grav,zs(i),rrr,tp(i,1),ps(i)"
           write(6,'(6f11.3)')pmsl(i),grav,zs(i),rrr,tp(i,1),ps(i)
         endif

c***********************************************************************
c end of npts (xy) loop
       enddo ! i
      endif ! if ( have_gp ) then       
      print *,'calc ps(1),(npts)=',ps(1),ps(npts)
      call amap ( ps, imt, jmt, 'calculated ps', 4., 1.e3 )
      call prt_pan(ps,il,jl,2,'ps(Pa)')
      call prt_pan(pmsl,il,jl,2,'pmsl(hPa)')
      call prt_pan(tp(1,1),il,jl,2,'temp(K)')
c***********************************************************************
c interpolate hs,rs,us,vs,ts
c   use linear or cubic spline interpolation in the vertical with lm
c     roundoff values occur in levels 1-5 of hs, rs when spline = t ,
c     consequently hs, rs are always interpolated linearly
      write(6,'(3x,"vertical interp.: linear (with heights of"
     &      ," levels expressed as log p): ps(1)=",f10.2)') ps(1)*.01
c***********************************************************************

c xy loop
      do i=1,npts

c sigma level loop ( top down here )
        do k=1,lm

          sigp=(ps(i)-ptop)*sgml(k)+ptop

c pressure level loop ( top down )
          do lev=1,nplevsm ! (=nplevs-1)
            pr1 = pm(1)
            prx = pm(nplevs)
            if ( osig_in ) then
              pr1 =calc_p(psg_m(i),pm(1))
              prx =calc_p(psg_m(i),pm(nplevs))
            endif ! ( osig_in ) then

            if ( sigp.lt.pr1 ) then ! sigp above top press
c assumes constant values above top input pressure
              hs(i,k)=hp(i,1) ! set to highest data value 
              rs(i,k)=rp(i,1)
              us(i,k)=up(i,1)
              vs(i,k)=vp(i,1)
              ts(i,k)=tp(i,1)
              if ( i.eq.1 ) then
                write(6,'("top index,pres=",i3,f12.3)')1,pr1
                write(6,'("sig index,sigp=",i3,f12.3)')k,sigp
              endif
              go to 361 ! go to next sigma
            endif ! ( sigp.lt.pr1 ) ! sigp above top press

            if ( sigp.gt.prx ) then ! sigp below bot press
              hs(i,k)=hp(i,nplevs) ! set to lowest data value 
              rs(i,k)=rp(i,nplevs)
c assumes 6.5 deg/km lapse below lowest input press.
              ts(i,k)=tp(i,nplevs)-6.5e-3*(alog(prx/sigp))
     .                  *(rrr*tp(i,nplevs))/grav
             if ( zerowinds ) then
c assumes 0. wind at sfc.
               fac=1.-(prx-sigp)/(prx-ps(i))
               us(i,k)=up(i,nplevs)*fac
               vs(i,k)=vp(i,nplevs)*fac
             else
c uses lowest input level winds
               us(i,k)=up(i,nplevs)
               vs(i,k)=vp(i,nplevs)
             endif ! ( zerowinds ) then
              if ( i.eq.1 ) then
                write(6,'("sig index,sigp=",i3,f12.3)')k,sigp
                write(6,'("bot index,pbot=",i3,f12.3)')nplevs,prx
              endif
              go to 361 ! go to next sigma
            endif ! ( sigp.gt.prx) ! sigp above top press

c set top/bot input pressures
            prest=pm(lev) ! top pressure
            presb=pm(lev+1) ! bottom pressure
            if ( osig_in ) then
              prest = calc_p(psg_m(i),pm(lev))
              presb = calc_p(psg_m(i),pm(lev+1))
            endif ! ( osig_in ) then
c test to see if sigp between input pressures
            if ( sigp.lt.prest .or. sigp.gt.presb ) go to 360 ! go to next press
            if ( i.eq.1 ) then
              write(6,'("ptop index, pres=",i3,f12.3)') lev,prest
              write(6,'("sigp index, sigp=",i3,f12.3)') k,sigp
              write(6,'("pbot index, pbot=",i3,f12.3)') lev+1,presb
            endif
c set up linear interpolation
            asigp=alog(sigp)
              if ( osig_in ) then
                alp =alog(calc_p(ps(i),pm(lev)))
                alpp=alog(calc_p(ps(i),pm(lev+1)))
              else
                alp =alpm(lev)
                alpp=alpm(lev+1)
              endif ! ( osig_in ) then
            !asp_abp=asigp-alpm(lev+1)
            !atp_abp=alpm(lev+1)-alpm(lev)
            asp_abp=asigp-alpp
            atp_abp=alpp-alp
            fap=asp_abp/atp_abp
c always linear interpolate rh and mix.ratio
            hs(i,k)=hp(i,lev+1)+(hp(i,lev+1)-hp(i,lev))*fap
            rs(i,k)=rp(i,lev+1)+(rp(i,lev+1)-rp(i,lev))*fap

            if(.not.splineu)then
                us(i,k)=up(i,lev+1)+(up(i,lev+1)-up(i,lev))*fap
            endif ! not splineu

            if(.not.splinev)then
                vs(i,k)=vp(i,lev+1)+(vp(i,lev+1)-vp(i,lev))*fap
            endif ! not splinev

            if ( .not. splinet ) then
              ts1=tp(i,lev+1)+(tp(i,lev+1)-tp(i,lev))*fap
c convert back to sensible temperature
              ts(i,k)=ts1*0.622*(1.0+rs(i,k))/(0.622+rs(i,k))
	      !-------------------------------------------------------
	      ! MJT suggestion
	      if (ts(i,k).gt.350.) then
               print *,"bad ts at ",i,k,ts(i,k)
               ts(i,k)=tpold(i,lev+1)+(tpold(i,lev+1)-tpold(i,lev))*fap
               print *,"new ts at ",i,k,ts(i,k)
	      end if
	      !-------------------------------------------------------
            endif ! not splinet

c end of pressure loop
 360        continue
          enddo ! lev=1,nplevsm ! (=nplevs-1)

c end of sigma loop
 361      continue
        enddo ! k=1,lm

c end of x/y loop
      enddo ! i=1,npts
c***********************************************************************
c  vertical interpolation with spline as linear p
        write(6,'(3x,"vertical interp.: spline (with heights of"
     &      ," levels expressed as linear p)")')

        if ( splinet ) then

c temperature ( virtual )
          call vispl ( tp,ts,pm,ps,npts,ifull,sgml,lm,ptop,nplevs )

c extrapolate temps for press > pm using 6.5/km lapse
c note that pm(nplevs)=bottom data level pressure
c and ts(lm) = bottom level sigma

          do k=1,lm
           do i=1,npts
            sigp=sgml(k)*ps(i)+ptop
            prx = pm(nplevs)
            if ( osig_in ) prx=calc_p(psg_m(i),pm(nplevs))
            if ( sigp.gt.prx ) then
              ts(i,k)=tp(i,nplevs)-6.5e-3*(alog(prx/sigp))
     .               *(rrr*tp(i,nplevs))/grav
            endif ! p>pm
           end do ! i=1,npts
          end do ! k=1,lm
c convert back to sensible temperature
          do k=1,lm
            do i=1,npts
              ts(i,k)=ts(i,k)*.622*(1.+rs(i,k))/(.622+rs(i,k))
            end do ! i=1,npts
          end do ! k=1,lm

        endif ! ( splinet ) then

c u component of the wind
        if(splineu)then
          call vispl(up,us,pm,ps,npts,ifull,sgml,lm,ptop,nplevs)
        endif!(splineu)then

c v component of the wind
        if(splinev)then
          call vispl(vp,vs,pm,ps,npts,ifull,sgml,lm,ptop,nplevs)
        endif!(splinev)then

c extrapolate winds for press > pm (assuming usfc=0 if zerowinds=t)
c note that pm(nplevs)=bottom data level pressure
c and ps(lm) = bottom level sigma
         do k=1,lm
          do i=1,npts
           sigp=sgml(k)*ps(i)+ptop
           prx = pm(nplevs)
           if ( osig_in ) prx=calc_p(psg_m(i),pm(nplevs))
           if ( sigp.gt.prx ) then
             if ( zerowinds ) then
               fac=1.-(prx-sigp)/(prx-ps(i))
               if(splineu)us(i,k)=up(i,nplevs)*fac
               if(splinev)vs(i,k)=vp(i,nplevs)*fac
             else ! constant winds
               if(splineu)us(i,k)=up(i,nplevs)
               if(splinev)vs(i,k)=vp(i,nplevs)
             endif ! ( zerowinds ) then
           endif ! p>pm
          end do ! i=1,npts
         end do ! k=1,lm

c***********************************************************************
c temperature adjustment for dry adiabat
c  lcoun  = counts # of grid columns  where temp lapse rate exceeded
c           dry adiab on first pass
c  ncoun  = counts # of grid points (actually layers) in vertical
c           at each horiz gridpoint where lapse exceeded dry adiab
c***********************************************************************
      if ( debug ) then
      write(6,96)
 96   format(3x,'start temperature adjustment for dry adiabats')
      endif

      print *,'adiabatic adjustment using new method (dryadj)'
      lcoun=0

      call dryadj ( ts, ps, ptop, npts, lcoun, ifull, dsg, sgml, lm )

      ltcoun = ltcoun+lcoun

c     write(6,5310)lcoun
c5310 format(3x,'gridpoints where temp lapse rate exceeded dry'
c    .  ,'adiabatic=',i6)

c***********************************************************************

      if ( nrh.ne.0 ) then
c top nrh layers are assumed dry
         do k=1,nrh
           do iq=1,npts
             hs(iq,k)=c0
             rs(iq,k)=c0
           enddo ! iq=1,npts
         enddo ! k=1,nrh
      endif ! nrh

c***********************************************************************

c limit humidity to 0 -> 100 %

      nrhp=nrh+1
      do 745 k=nrhp,lm
        do 745 iq=1,npts
          hs(iq,k) = min ( max ( hs(iq,k),c0 ), c100 )
 745  continue

c***********************************************************************

c recomputation of mixing ratio based on new temperature

c recalculate mix.rat. (rs)
c     input : hs=rel.hum.(%) , ts=temps.(k) , tp=press.(pa)
c     output: rs=mix.rat.(g/g)   all npts

c tp is temporarily used to hold the pressures (pa)
c rp,zp,up,vp are used as temporary arrays

      do 748 k=nrhp,lm

        do 746 iq=1,npts
 746      tp(iq,1)=(ps(iq)-ptop)*sgml(k)+ptop

c ts is returned with h2o saturation partial pressures
        call esmtrv ( ts(:,k), rs(:,k), npts, rp(:,1), zp(:,1),
     &                up(:,1), vp(:,1) )

c calculate actual.mix.ratio
        do 750 iq=1, npts
          qs=0.01*hs(iq,k)*rs(iq,k)
          rs(iq,k)=0.622*qs/(tp(iq,1)-rs(iq,k))
 750    continue

 748  continue

c***********************************************************************

c limit mix.rat. to 0.0 < mix.rat. < 50.e-3 g/g

      do 749 k=nrhp,lm
        do 749 iq=1,npts
          rs(iq,k) = min ( max ( rs(iq,k),c0 ) , 50.e-3 )
 749  continue

c***********************************************************************

      write(6,'(1x,"sigma-level data at 1st and last point :")')
      i=1
      write(6,9056)(k,us(i,k),vs(i,k),ts(i,k),rs(i,k),hs(i,k),k=1,lm)
      i=npts
      write(6,9056)(k,us(i,k),vs(i,k),ts(i,k),rs(i,k),hs(i,k),k=1,lm)
 9056 format(9x,'us',11x,'vs',11x,'ts',11x,'rs',11x,'hs'/(i3,1x,5e13.5))

c***********************************************************************

      write(6,99)ltcoun,pcoun
   99 format(/
     & 1x,'total nr of horiz gridpoints exceeded  dry adiabatic=',i10/
     & 1x,'total nr of horiz non-convergnent gridpoints for temp adjust'
     & ,'=',i10)

c***********************************************************************

c write sigma data to disk

c***********************************************************************
c read replacement sfct if ints>0
         if ( ints.ne.0 ) then
            print *,'read new ts from unit=',ints,' npts=',npts
            read ( abs(ints),* ) (sfct(i),i=1,npts)
            call amap ( sfct, imt, jmt, 'new sfct', 5., 0. )
         endif

        call maxmin(ts,'ts',0,1.,il,kl)
        call maxmin(rs,'rs',0,1.,il,kl)
        call maxmin(us,'us',0,1.,il,kl)
        call maxmin(vs,'vs',0,1.,il,kl)
        
        print *,"*** ",maxval(ts),minval(ts)

c write out file
        if ( io_out .eq. 3 ) then
           write(6,*)"io_out=3 no longer supported!!!!!!!!!!"
           stop
        else
           call invert1(sgml,kl)
           call invert3(ts,il,kl)
           call invert3(rs,il,kl)
           call invert3(us,il,kl)
           call invert3(vs,il,kl)

           call outcdf(ihr,idy,imon,iyr,iout,nt,time,mtimer
     &                 ,sgml,vfil,ds,il,kl)

           call invert1(sgml,kl)

        endif
c***********************************************************************
      return ! vidar
      end ! vidar
!=======================================================================
      subroutine calcsig(sigh,kl)

c these routines are written from top down
c nsig=50 option general cubic half-sigma levels (normal formula bottom up)

      integer kl

      real sigh(kl+1)

      sigt=.008916
      den=.75*( 2./kl-1. - (2./kl-1.)**3 )
      alf=(sigt-.5 -.5*(2./kl-1.)**3)/den

      do k=0,kl
        rk=(.5*(kl+1)-k-.5)/kl
        sigh(kl+1-k)=.5 + 1.5*alf*rk - 2.*(3.*alf-2.)*rk**3
      enddo

      return
      end
!=======================================================================
      function calc_p(ps,pl)
      real calc_p, ps, pl
! compute actual press (Pa) using ps(Pa) and sigma (pl=sigma*100000)
      calc_p = ps*(pl/1.e5) ! compute actual press usng sfc p and sigma
      return
      end
!=======================================================================
