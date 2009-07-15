      subroutine tomod ( alm, sgml, prein
     &                ,ihr,idy,imn,iyr,iout,oform
     &                ,outfil,ddss,opre)

! version for globpe   stagger winds here
c-----------------------------------------------------------------------

c  this program writes out the interpolated file in a format
c  consistent with the input to the nested model

c-----------------------------------------------------------------------

      include 'newmpar.h' ! il,jl,kl,ifull,ijk
      include 'map.h'
      include 'parm.h'
      parameter ( ip=il, jp=jl, ipjp=ifull )
      parameter ( klp=kl+1 )
      parameter ( cfact=1.45444e-4, pi=3.1415926536 )
      parameter ( dtr=pi/180. )

      real alm(ipjp),prein(ipjp)

      include 'sigdata.h'
!n    common/sigdata/pmsl(ifull),sfct(ifull),zs(ifull),ps(ifull)
!n   &             ,us(ifull,kl)    ,vs(ifull,kl)    ,ts(ifull,kl)
!n   &             ,rs(ifull,kl)    ,hs(ifull,kl)    ,psg_m(ifull)
!n   &             ,zsi_m(ifull)
!o    common/sigdata/pmsl(ifull),ts(ifull),zs(ifull),ps(ifull)
!o   &             ,u(ifull,kl)    ,v(ifull,kl)    ,t(ifull,kl)
!o   &             ,qg(ifull,kl)   ,temp(ifull,kl)

c     real swb(ipjp),swt(ipjp),stt(ipjp),stb(ipjp)
      real sgml(kl)

      real precip(ifull),preold(ifull)
      real psf(ifull)
      real zss(ifull)
      real tss(ifull),wt(ifull),wb(ifull)
      common/worktst/tsb(ifull),tst(ifull)

      character*80 outfil,oldout
      character*10 rundate

      logical tvirt,oform,ofirst,opre
      data oldout/'              '/
      data tvirt/.true./,ofirst/.true./
      data rundate/"ncepavnanl"/

      save

      print *,'**** subroutine tomod ****, oform=',oform
      print *,'ipjp,ifull,opre=',ipjp,ifull,opre

      call maxmin(t,'t',0,1.)

c***********************************************************************
      print *,'ds,ofirst=',ds,ofirst
c     if ( ofirst ) then     
c       call lconset ( ds )
        do i=1,ifull
          preold(i)=0.
        enddo
c     endif ! ofirst

c invert vertical levels in hs() like darlam (bottom up)
      do k=1,kl
       do n=1,ifull
        hs(n,k)=ts(n,kl+1-k)
       enddo ! n=1,ifull
       print *,'tsig(k)=',k
       call findxn ( hs(1,k), ifull, -1.e29, amax, k1, amin, k2 )
      enddo
      call maxmin(hs,'temp',0,1.)
c***********************************************************************
      print *,'zs'
      call findxn ( zs, ipjp, -1.e29, amax, k1, amin, k2 )
      print *,'ps'
      call findxn ( ps, ipjp, -1.e29, amax, k1, amin, k2 )
      print *,'ts'
      call findxn ( ts, ipjp, -1.e29, amax, k1, amin, k2 )
      print *,'alm'
      call findxn ( alm, ipjp, -1.e29, amax, k1, amin, k2 )
      n=0
      do j=1,jl
        do i=1,il
          iq = i+il*(j-1)
          n=n+1
c this assumes that precip is accumulated for each time period
          if ( opre ) then
            precip(n)= max(0.,prein(iq)-preold(n))
          else
            precip(n)= 0.
          endif
c set zs to g*sfc.topog
          zss(n)   = 9.80616*zs(iq)
c psf = log(1.e-5*sfc p)
          psf(n)   = log(1.e-5*ps(n))
          tss(n)   = sfct(iq)

c         tst(n)   = stt(iq)
c         tsb(n)   = stb(iq)
c         wt(n)    = swt(iq)
c         wb(n)    = swb(iq)

          tst(n)   = hs(n,2)
          tsb(n)   = hs(n,2)
          wt(n)    = .15
          wb(n)    = .15
c set tss neg if ocean point
          if ( alm(iq).lt..5 ) then
             tss(n)=-tss(n)
          endif ! alm(iq).lt..5 ) then
!         rj=float(j)
!         ri=float(i)
c calc. map fac=em and coriol fac = f
!         call lconll ( alon, alat, ri, rj )
!         call mapff  ( alat, facmap, coriol )
!          if(ofirst.and.i.eq.1.and.j.eq.1)then
!            print *,'j,alat,alon,facmap,coriol='
!            print *,j,alat,alon,facmap,coriol
!          endif
!          em(n) = facmap
!          f (n) = coriol
        end do ! i=1,il
      end do ! j=1,jl
c***********************************************************************

      print *,'ifull,n=',ifull,n
      if ( ifull.ne.n ) then
         print *,'********************* ifull<>n'
         stop 'ifull<>n'
      endif ! ifull.ne.n ) then
      ofirst=.false.

      do i=1,ifull
        preold(i)=precip(i)
      enddo

c***********************************************************************
c temporary hard-coding for control block
      kdate = iyr*10000+imn*100+idy
      ktime = ihr*100
      print *,'kdate,ktime=',kdate,ktime
      ktau  = 0
      ik    = il
      jk    = jl
      kk    = kl
      m     = 0
      nsd   = 0
      nbd   = 0
      nps   = 0
      mex   = 0
      mup   = 0
      nem   = 0
      nsi   = 0
      nmi   = 0
      ndt   = 0
      npsav = 0
      nhor  = 0
      nkuo  = 0
      khdif = 0
      kwt   = kl
      iaa   = 1
      jaa   = 1
      timeb = 0.
      timelb= 0.
      nvad  = 0
      ndum  = 0
      difk  = 0
      rhk   = 0
      nqg   = 3
      ntsur = 3

      write(6,*)"outfil=",outfil
      write(6,*)"oldout=",oldout
      if ( outfil.ne.oldout ) then
          meso=3
          write(6,'("open unformatted disk (",i3,") called ",a40)')
     &            iout,outfil
          open(iout,file=outfil,form='unformatted',status='unknown')
        oldout=outfil
      else ! outfil.eq.oldout
         print *,'writting to file=',outfil
      endif ! outfil.ne.oldout

c unformatted writes
            ivirt=0
            if ( tvirt ) ivirt=1

c write contrl block
            write(iout) kdate,ktime,ktau,ik,jk,kk,m,nsd,meso,nbd,nps
     .             ,mex,mup,nem,nsi,nmi,ndt,npsav
     .           ,rundate,nhor,nkuo,khdif,kwt,iaa,jaa,timeb,timelb
     .             ,ds,nvad,nqg
     .           ,(ndum,nn=1,10),ivirt,ntsur,(ndum,n2=1,8)
     .             ,difk,rhk,du,tanl,rlong0,rlat0,schmidt
            print *,kdate,ktime,ktau,ik,jk,kk,m,nsd,meso,nbd,nps
     .             ,mex,mup,nem,nsi,nmi,ndt,npsav
     .           ,rundate,nhor,nkuo,khdif,kwt,iaa,jaa,timeb,timelb
     .             ,ds,nvad,nqg
     .           ,(ndum,nn=1,10),tvirt,ntsur,(ndum,n2=1,8)
     .             ,difk,rhk,du,tanl,rlong0,rlat0,schmidt

            write(iout) (sgml(klp-k),k=1,kl)
            write(6,'(''sig='',(10f8.5))') (sgml(klp-k),k=1,kl)

! scaled sfc press
            write(iout) (psf(n),n=1,ifull)
            write(6,'(''psf='',2f8.5)') psf(1),psf(ifull)

c mean sea level pressure
            write(6,*)"sgml=",sgml
            write(6,*)"ps(5),zss,temp,ifull,kl=",
     &                ps(5),zss(5),hs(5,2),ifull,kl
            is=il/2+(il+il/2-1)*il
            write(6,*)"ps=",(ps(is+i),i=1,5)
            write(6,*)"zss=",(zss(is+i),i=1,5)
            write(6,*)"temp=",(hs(is+i,2),i=1,5)

            xpmsl=-1.
            do iq=1,il*jl
              xpmsl=max(xpmsl,pmsl(iq))
            enddo

            if ( xpmsl .lt. 100.e2 ) then
              write(6,*)"xpmsl=",xpmsl
              write(6,*)"mslp(ps,pmsl,zss,temp,sgml,ifull,ifull,kl)"
              call mslp(ps,pmsl,zss,hs,sgml,ifull,ifull,kl)
            endif

            write(6,*)"pmsl=",(pmsl(is+i),i=1,5)

            write(iout) (pmsl(n),n=1,ifull)
            write(6,'(''pmsl '',2f8.0)') pmsl(1),pmsl(ifull)

            write(iout) (zss(n),n=1,ifull)
            write(6,'(''zs '',2f8.0)') zss(1),zss(ifull)

            write(iout) (em(n,1),n=1,ifull)
            write(6,'(''em '',2f8.4)') em(1,1),em(il,jl)

            write(iout) (f(n,1),n=1,ifull)
            write(6,'(''f '',4p2f8.4)') f(1,1),f(il,jl)

            write(iout) (tss(n),n=1,ifull)
            write(6,'(''tss '')')
            call findxn(tss,ifull,-1.e29,ax,ix,an,in)

            if ( nqg.ge.2 ) then

              write(iout) (precip(n),n=1,ifull)
              write(6,'(''precip '')')
              call findxn(precip,ifull,-1.e29,ax,ix,an,in)

              if ( nqg.ge.3 ) then

                if ( ntsur.gt.2 ) then
                  write(iout) (tsb(n),n=1,ifull)
                  print *,'tsb'
                  call findxn(tsb,ifull,-1.e29,ax,ix,an,in)
                endif ! ( ntsur.gt.2 ) then

                write(iout) (tst(n),n=1,ifull)
                print *,'tst'
                call findxn(tst,ifull,-1.e29,ax,ix,an,in)

                write(iout) (wt(n),n=1,ifull)
                print *,'wt'
                call findxn(wt,ifull,-1.e29,ax,ix,an,in)

                write(iout) (wb(n),n=1,ifull)
                print *,'wb'
                call findxn(wb,ifull,-1.e29,ax,ix,an,in)
              end if
            end if

         write(iout) ((hs(n,k),n=1,ifull),k=1,kl)
         write(6,'("temp written")')
         do k=1,kl
           call findxn(hs(1,k),ifull,-1.e29,ax,ix,an,in)
         enddo ! k=1,kl

c        don't stagger the winds for globpea (for writing-out purposes only)
!        print *,'before staguv - not done because now globpea'
!        call maxmin(us,' us',0,1.)
!        call maxmin(vs,' vs',0,1.)
c        do k=1,kl
c          print *,'before staguv; k,u,v: ',k,us(1,k),vs(1,k)
c          do iq=1,ifull
c            tsb(iq)=us(iq,k)   ! tsb just as temporary storage
c            tst(iq)=vs(iq,k)   ! tst just as temporary storage
c          enddo      ! iq loop
c          call staguv(tsb,tst,us(1,k),vs(1,k))  ! *** need tsb & tst contiguous
c          print *,'after staguv; k,u,v: ',k,us(1,k),vs(1,k)
c        enddo  ! k loop
c        print *,'after staguv'
c        call maxmin(us,' us',0,1.)
c        call maxmin(vs,' vs',0,1.)

         write(iout) ((us(n,klp-k),n=1,ifull),k=1,kl)
         write(6,'("us written")')
         do k=1,kl
           call findxn(us(1,k),ifull,-1.e29,ax,ix,an,in)
         enddo ! k=1,kl

         write(iout) ((vs(n,klp-k),n=1,ifull),k=1,kl)
         write(6,'("vs written")')
         do k=1,kl
           call findxn(vs(1,k),ifull,-1.e29,ax,ix,an,in)
         enddo ! k=1,kl

c convert mix.rat. to split grid
         write(iout) ((rs(n,klp-k),n=1,ifull),k=1,kl)
         write(6,'("rs written")')
         do k=1,kl,1
           call findxn(rs(1,k),ifull,-1.e29,ax,ix,an,in)
         enddo ! k=1,kl

         print *,'initial conditions for ktau = ',ktau

         return
         end
