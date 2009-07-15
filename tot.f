      subroutine tot ( d, t, nl )

      include 'newmpar.h'
      parameter ( ijlp = (il+1)*(jl+1) )

      dimension t(ij,nl)
      dimension d(ijlp,nl)

         do l=1,nl
           do j=1,jl
             do i=1,il
               inj=i+(il+1)*(j-1)
               n=i+il*(j-1)
               t(n,l)=d(inj,l)
             end do ! i=1,il
           end do ! j=1,jl
         end do ! l=1,nl

      return

      entry tou ( d, t, nl )

         do l=1,nl
           do j=1,jl
             do i=1,il
               inj=i+(il+1)*(j-1)
               n=i+il*(j-1)
               t(n,l)=.5*(d(inj,l)+d(inj+1,l))
             end do ! i=1,il
           end do ! j=1,jl
         end do ! l=1,nl

      return

      entry tov ( d, t, nl )

         do l=1,nl
           do j=1,jl
             do i=1,il
               inj=i+(il+1)*(j-1)
               n=i+il*(j-1)
               t(n,l)=.5*(d(inj,l)+d(inj+il+1,l))
             end do ! i=1,il
           end do ! j=1,jl
         end do ! l=1,nl

      return
      end
