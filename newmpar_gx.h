!     plan to replace ij by ifull eventually
      parameter(il=48,npanels=5,jl=il+npanels*il,kl=18,ksl=3)
!     parameter(il=120,npanels=5,jl=il+npanels*il,kl=18,ksl=3)
!     parameter(il=200,npanels=5,jl=il+npanels*il,kl=18,ksl=3)
      parameter(ifull=il*jl,ij=il*jl,ijk=il*jl*kl)
      parameter( iquad=1+il*((8*npanels)/(npanels+4)) )
!     for     npanels:   0          5        13
!                  jl:   -         6*il     14*il
!               iquad:   1         4*il+1   6*il+1
