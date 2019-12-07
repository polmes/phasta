c ---------------------------------------------------------------------
c                       SLIP BOUNDARY CONDITIONS
c                       by Pol Mesalles Ripoll
c                       CU Boulder, Fall 2019
c ---------------------------------------------------------------------

      ! include "common.h"

      module slipbc
         implicit none

            ! Allocatable variables, constant parameters here
            ! (public / private)

         contains

            function testhere()
               character(len = 8) :: testhere
               testhere = "Hell-oh!"
               ! return ! just a quick test that we can load the module
            end function testhere

            subroutine getSlipVelocity1D(mu, rho, T, dudy, uslip)
c           -----------------------------------------------------
c           Computes Maxwell 1D slip velocity at the wall
c           Input:
c              mu:   dynamic viscosity
c              rho:  density at the integration point
c              T:    temperature at the integration point
c              dudy: du1/dx2 velocity gradient
c           Output:
c              uslip: slip velocity at the integration point
c            Note: pi, Rgas are constants defined in common.h
c           -----------------------------------------------------

               real*8, intent(in) :: mu
               real*8, dimension(2) :: mfp
               real*8, dimension(2) :: rho, T, dudy, uslip

               print *, "Address In", loc(mu)
               print *, "Address In", loc(uslip)

               ! MFP = Mean Free Path, typically referred to as \lambda
               mfp = mu / rho * sqrt(3 / (2 * 287 * T))

               ! Slip velocity at the wall
               uslip = mfp * dudy

               return
            end subroutine getSlipVelocity1D

      end module slipbc

      program slipbctest
         use slipbc

         implicit none

         real*8 mu, rho(2), T(2), dudy(2), uslip(2)
         mu = 1
         rho = (/1, 2/)
         T = (/200, 200/)
         dudy = (/0.5, 0.5/)
         uslip = (/0, 0/)

         print *, "Address Before", loc(mu)
         print *, "Address Before", loc(uslip)

         call getSlipVelocity1D(mu, rho, T, dudy, uslip)

         print *, "Address After", loc(mu)
         print *, "Address After", loc(uslip)

         print *, "This is a test", uslip
      end program slipbctest
