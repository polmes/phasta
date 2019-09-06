      subroutine AddEddyViscKomgq(blk,ith,yl, dwl, shdrv, shape, rho, rmu, ssq,
     &                               xl)
      use eblock
      use turbKW ! access to KW model constants
      include "common.h"
      type (LocalBlkData) blk

      real*8 yl(bsz,blk%s,ndof), xl(bsz,blk%n,nsd), dwl(bsz,blk%n),
     &       shape(blk%e,blk%s), rho(blk%e), ssq(blk%e), 
     &       shdrv(blk%e,nsd,blk%s), rmu(blk%e)
      real*8 kay, omega, dw, ssqInv, 
     &       nu, Rou, F1, F2, mut, gradK(blk%e,nsd),  
     &       gradW(blk%e,nsd), delKdelW, mutden
      integer e,n,ith

      rho = datmat(1,1,1)
      rmu = datmat(1,2,1)
      mut = zero
      call e3qvarkwSclr(blk,ith,yl, shdrv, xl, gradK, gradW)
      do e = 1, blk%e
         kay = zero
         omega = zero
         dw = zero
         do n = 1, blk%s
            kay = kay + shape(e,n)*yl(e,n,6)
            omega = omega + shape(e,n)*yl(e,n,7)
            ! This won't work for oder p>1
            dw = dw + shape(e,n)*dwl(e,n)
         enddo
         kay = max(kay,zero)
         omega = max(omega,zero)

         nu = rmu(e)/rho(e)
         Rou = rho(e)
         delKdelW = gradK(e,1)*gradW(e,1)  +  gradK(e,2)*gradW(e,2)
     &          +    gradK(e,3)*gradW(e,3)

         call getblendfunc (delKdelW, kay, omega, dw, Rou, nu, F1, F2)
         mutden = max(max(a1*omega,1.0d-10),max(F2*ssq(e),1.0d-10))
         mut = Rou * a1 *kay / mutden
         mut = max(mut,1.0d-4*rmu(e))
         rmu(e) = rmu(e)+mut
 
      enddo

      return
      end subroutine AddEddyViscKomgq
