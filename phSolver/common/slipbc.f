c     ------------------------------------------------------------------
c                       SLIP BOUNDARY CONDITIONS
c                       by Pol Mesalles Ripoll
c                       CU Boulder, Fall 2019
c     ------------------------------------------------------------------

      subroutine getSlipVelocity1D(mu, rho, T, dudy, uslip)
c     ------------------------------------------------------------------
c        Computes Maxwell 1D slip velocity at the wall
c        Input:
c           mu:   dynamic viscosity
c           rho:  density at the integration point
c           T:    temperature at the integration point
c           dudy: du1/dx2 velocity gradient
c        Output:
c           uslip: slip velocity at the integration point
c        Note: pi, Rgas, npro are constants defined in common.h
c     ------------------------------------------------------------------

         include "common.h"

         real*8, intent(in) :: mu, rho(npro), T(npro), dudy(npro)
         real*8, intent(out) :: uslip(npro)
         real*8 :: mfp(npro) ! or dimension(npro)

         ! MFP = Mean Free Path, typically referred to as \lambda
         mfp = mu / rho * sqrt(pi / (2 * Rgas * T))

         ! Slip velocity at the wall
         uslip = mfp * dudy

         return
      end subroutine getSlipVelocity1D
