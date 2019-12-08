c     ------------------------------------------------------------------
c                       SLIP BOUNDARY CONDITIONS
c                       by Pol Mesalles Ripoll
c                       CU Boulder, Fall 2019
c     ------------------------------------------------------------------

      subroutine getSlipVelocity1D(mu, rho, T, dudy, uslip)
c     ------------------------------------------------------------------
c        Computes Maxwell 1D slip velocity at the wall
c        Input:
c           mu   - dynamic viscosity
c           rho  - density at the integration point
c           T    - temperature at the integration point
c           dudy - du1/dx2 velocity gradient
c        Output:
c           uslip - slip velocity at the integration point
c        Note: pi, Rgas, npro, slipSigma, slipConst
c              are constants/parameters defined in common.h
c     ------------------------------------------------------------------

         include "common.h"

         real*8, intent(in) :: mu, rho(npro), T(npro), dudy(npro)
         real*8, intent(out) :: uslip(npro)
         real*8 :: mfp(npro) ! or dimension(npro)

         ! MFP = Mean Free Path, typically referred to as \lambda
         mfp = mu / rho * sqrt(pi / (2 * Rgas * T))

         ! Slip velocity at the wall
         uslip = slipConst * (2 - slipSigma) / slipSigma * mfp * dudy

         return
      end subroutine getSlipVelocity1D

      subroutine setSlipBC(iBC, iBCB, ienb)
c     ------------------------------------------------------------------
c        Corrects iBCB to account for Nitsche BC
c        Input:
c           iBC  - essential boundary condition codes
c           iBCB - natural boundary condition codes
c           ienb - nodal connectivity for boundary elements (in block)
c           dudy: du1/dx2 velocity gradient
c        Output:
c           iBCB - modified natural boundary condition codes
c        Note: nshg, npro are defined in common.h
c     ------------------------------------------------------------------

         include "common.h"

         integer, intent(in) :: iBC(nshg), ienb(npro, nshl)
         integer, intent(inout) :: iBCB(npro,ndiBCB)
         integer :: el

         do el = 1, npro
            if (any(btest(iBC(ienb(el,:)), 6))) then
               iBCB(el,1) = iBCB(el,1) + 32 ! enable sclr1 flux (2^5)
               ! iBCB(el,2) = 1 ! surfID = 1
            endif
         enddo

      end subroutine setSlipBC
