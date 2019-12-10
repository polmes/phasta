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

      subroutine setSlipiBC(iBC, iBCB, ienb)
c     ------------------------------------------------------------------
c        Corrects iBCB to account for Nitsche BC
c        Input:
c           iBC  - essential boundary condition codes
c           iBCB - natural boundary condition codes
c           ienb - nodal connectivity for boundary elements (in block)
c           dudy - du1/dx2 velocity gradient
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

      end subroutine setSlipiBC

      subroutine getNodalShapeFunctions(shpnod, shglnod)
c     ------------------------------------------------------------------
c        Computes shape functios at element boundary nodes
c        (NOT at the quadrature points)
c        Output:
c           shpnod  - Array of shape functions at boundary nodes
c           shglnod - Array of local shape function gradients at
c                     boundary nodes
c     ------------------------------------------------------------------

         include "common.h"

         real*8, intent(out) :: shpnod(nenbl, nshl),
     &                          shglnod(nenbl, nsd, nshl)

         ! Only 1st order shape functions implemented for now
         if (ipord .gt. 1) then
            call error('getNodalShapeFunctions',
     &                 'Order too high - not implemented',
     &                 ipord)
         end if

         ! Check element topology
         select case (lcsyst) ! <toppology cases>
            case (1) ! tet's
               ! 3 boundary vertices, bouondary nodes first
               call shpTet(ipord, (/1, 0, 0/), shpnod(1, 1:3),
     &                     shglnod(1, :, 1:4))
               call shpTet(ipord, (/0, 1, 0/), shpnod(2, 1:3),
     &                     shglnod(2, :, 1:4))
               call shpTet(ipord, (/0, 0, 1/), shpnod(3, 1:3),
     &                     shglnod(3, :, 1:4))

            case (2) ! hex's
               ! 4 boundary vertices, boundary nodes first
               ! TODO: check that order (parent space) is correct
               call shpHex(ipord, (/-1, -1, -1/), shpnod(1, 1:8),
     &                     shglnod(1, :, 1:8))
               call shpHex(ipord, (/+1, -1, -1/), shpnod(2, 1:8),
     &                     shglnod(2, :, 1:8))
               call shpHex(ipord, (/+1, +1, -1/), shpnod(3, 1:8),
     &                     shglnod(3, :, 1:8))
               call shpHex(ipord, (/-1, +1, -1/), shpnod(4, 1:8),
     &                     shglnod(4, :, 1:8))

c     TODO: instead of explicitly passing the local array of
c           nodal points, use  either a nodpt(MAXTOP,3,nenbl) array
c           or find the values using a for loop

            case default
               call error('getNodalShapeFunctions',
     &                    'Unknown element type - not implemented',
     &                    lcsyst)
         end select ! </toppology cases>

      end subroutine getNodalShapeFunctions

      subroutine slipCorrect(y)
c     ------------------------------------------------------------------
c        Modifies y variables to account for Dirichlet slip BC
c        Input:
c           y - Solution vector
c        Output:
c           y - Corrected solution vector
c     ------------------------------------------------------------------

         include "common.h"

         real*8, intent(in) :: y(nshg,nflow)
         integer :: vrt, nod
         real*8, allocatable :: shpnod(:,:), shpnodtmp(:,:),
     &                           shglnod(:,:,:), shglnodtmp(:,:,:)
         ! real*8, allocatable :: ycl(:,:,:)
         ! real*8, allocatable :: xl(:,:,:)
         ! real*8, allocatable :: yvl(:,:,:)

         ! Loop over boudary element blocks
         do iblk = 1, nelblb ! <loop blocks>
            iel    = lcblkb(1,iblk)
            ! lelCat = lcblkb(2,iblk)
            lcsyst = lcblkb(3,iblk)
            ! iorder = lcblkb(4,iblk)
            ! nenl   = lcblkb(5,iblk)  ! no. of vertices per element
            nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
            ! mattyp = lcblkb(7,iblk)
            ! ndofl  = lcblkb(8,iblk)
            nshl   = lcblkb(9,iblk)
            ! nshlb  = lcblkb(10,iblk)
            npro   = lcblkb(1,iblk+1) - iel
            ! if(lcsyst.eq.3) lcsyst=nenbl
            ! ngaussb = nintb(lcsyst)

            ! Get boundary nodes shape functions for block topology
            allocate(shpnodtmp(nenbl,nshl))
            allocate(shglnodtmp(nenbl,nsd,nshl))
            call getNodalShapeFunctions(shpnodtmp, shglnodtmp)

            ! Init local shape function arrays
            allocate(shpnod(npro,nshl))
            allocate(shglnod(npro,nshl,nsd))

            ! Loop over each boundary node
            do nod = 1, nenbl ! <loop nodes>
               ! Loop over all vertices or element nodes
               do vrt = 1, nshl
                  shpnod(:,vrt) = shpnodtmp(nod,vrt)
               end do

               ! Stuff here

            end do ! </loop nodes>
         enddo ! </loop blocks>

      endsubroutine slipCorrect
