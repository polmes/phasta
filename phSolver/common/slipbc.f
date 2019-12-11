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
            end if
         end do

      end subroutine setSlipiBC

      module slipGeometry
         real*8, allocatable :: xs(:,:)
      end module slipGeometry

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

         real*8, intent(out) :: shpnod(nenbl,nshl),
     &                          shglnod(nenbl,nsd,nshl)
         real*8, allocatable :: nodpt(:,:)
         integer :: nod

         ! Only 1st order shape functions implemented for now
         if (ipord .gt. 1) then
            call error('getNodalShapeFunctions',
     &                 'Order too high - not implemented',
     &                 ipord)
         end if

         ! Init array of nodes
         allocate(nodpt(nenbl,nsd))
         nodpt = zero

         ! Check element topology
         select case (lcsyst) ! <toppology cases>
            case (1) ! tet's
               ! 3 boundary vertices, bouondary nodes first
               ! TODO: check that order (in parent space) is correct
               nodpt(1,1) = 1
               nodpt(1,2) = 0
               nodpt(1,3) = 0
               nodpt(2,1) = 0
               nodpt(2,2) = 1
               nodpt(2,3) = 0
               nodpt(3,1) = 0
               nodpt(3,2) = 0
               nodpt(3,3) = 1

               do nod = 1, nenbl
                  call shpTet(ipord, nodpt(nod,:), shpnod(nod,:),
     &                        shglnod(nod,:,:))
               end do

            case (2) ! hex's
               ! 4 boundary vertices, boundary nodes first
               ! TODO: check that order (in parent space) is correct
               nodpt(1,1) = -1
               nodpt(1,2) = -1
               nodpt(1,3) = -1
               nodpt(2,1) = +1
               nodpt(2,2) = -1
               nodpt(2,3) = -1
               nodpt(3,1) = +1
               nodpt(3,2) = +1
               nodpt(3,3) = -1
               nodpt(4,1) = -1
               nodpt(4,2) = +1
               nodpt(4,3) = -1

               do nod = 1, nenbl
                  call shpHex(ipord, nodpt(nod,:), shpnod(nod,:),
     &                        shglnod(nod,:,:))
               end do

            case default
               call error('getNodalShapeFunctions',
     &                    'Unknown element type - not implemented',
     &                    lcsyst)
         end select ! </toppology cases>

      end subroutine getNodalShapeFunctions

      subroutine gradNodalShapeFunctions(shglnod, shgnod)
         use slipGeometry ! for global x(s)
         use pointer_data ! for ienb array

         include "common.h"

         real*8, intent(in) :: shglnod(npro,nshl,nsd)
         real*8, intent(out) :: shgnod(npro,nshl,nsd)
         real*8 :: dxdxib(npro,nsd,nsd), dxidxb(npro,nsd,nsd),
     &             temp(npro), xlb(npro,nshl,nsd)
         integer :: vrt, i, j

         ! Localize x node locations
         call localx(xs, xlb, mienb(iblk)%p, nsd, 'gather')

         ! Compute the deformation gradient
         dxdxib = zero
         do vrt = 1, nshl
            do i = 1, nsd
               do j = 1, nsd
                  dxdxib(:,i,j) = dxdxib(:,i,j) + xlb(:,vrt,i)
     &                            * shglnod(:,vrt,i)
               end do
            end do
         end do

         ! Compute the inverse deformation gradient
         dxidxb(:,1,1) = dxdxib(:,2,2) * dxdxib(:,3,3)
     &                   - dxdxib(:,3,2) * dxdxib(:,2,3)
         dxidxb(:,1,2) = dxdxib(:,3,2) * dxdxib(:,1,3)
     &                   - dxdxib(:,1,2) * dxdxib(:,3,3)
         dxidxb(:,1,3) = dxdxib(:,1,2) * dxdxib(:,2,3)
     &                   - dxdxib(:,1,3) * dxdxib(:,2,2)
         temp = one / (dxidxb(:,1,1) * dxdxib(:,1,1)
     &                 + dxidxb(:,1,2) * dxdxib(:,2,1)
     &                 + dxidxb(:,1,3) * dxdxib(:,3,1))
         dxidxb(:,1,1) = dxidxb(:,1,1) * temp
         dxidxb(:,1,2) = dxidxb(:,1,2) * temp
         dxidxb(:,1,3) = dxidxb(:,1,3) * temp
         dxidxb(:,2,1) = (dxdxib(:,2,3) * dxdxib(:,3,1)
     &                    - dxdxib(:,2,1) * dxdxib(:,3,3)) * temp
         dxidxb(:,2,2) = (dxdxib(:,1,1) * dxdxib(:,3,3)
     &                    - dxdxib(:,3,1) * dxdxib(:,1,3)) * temp
         dxidxb(:,2,3) = (dxdxib(:,2,1) * dxdxib(:,1,3)
     &                    - dxdxib(:,1,1) * dxdxib(:,2,3)) * temp
         dxidxb(:,3,1) = (dxdxib(:,2,1) * dxdxib(:,3,2)
     &                    - dxdxib(:,2,2) * dxdxib(:,3,1)) * temp
         dxidxb(:,3,2) = (dxdxib(:,3,1) * dxdxib(:,1,2)
     &                    - dxdxib(:,1,1) * dxdxib(:,3,2)) * temp
         dxidxb(:,3,3) = (dxdxib(:,1,1) * dxdxib(:,2,2)
     &                    - dxdxib(:,1,2) * dxdxib(:,2,1)) * temp

         ! Compute real space shape function gradients: N_{a,i}
         shgnod = zero
         do vrt = 1, nshl
            do i = 1, nsd
               do j = 1, nsd
               shgnod(:,vrt,i) = shgnod(:,vrt,i)
     &                           + shglnod(:,vrt,j)
     &                           * dxidxb(:,j,i)
               end do
            end do
         end do
      end subroutine gradNodalShapeFunctions

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
         integer :: vrt, nod, i
         real*8, allocatable :: shpnod(:,:), shpnodtmp(:,:),
     &                          shglnod(:,:,:), shglnodtmp(:,:,:),
     &                          shgnod(:,:,:)
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

            ! Allocate and init element shape function arrays
            allocate(shpnodtmp(nenbl,nshl))
            allocate(shglnodtmp(nenbl,nsd,nshl))
            shpnodtmp = zero
            shglnodtmp = zero

            ! Get boundary nodes shape functions for block topology
            call getNodalShapeFunctions(shpnodtmp, shglnodtmp)

            ! Allocate local shape function arrays
            allocate(shpnod(npro,nshl))
            allocate(shglnod(npro,nshl,nsd))
            allocate(shgnod(npro,nshl,nsd))

            ! Loop over each boundary node
            do nod = 1, nenbl ! <loop nodes>
               ! Init local shape function arrays
               shpnod = 0
               shglnod = 0
               shgnod = 0

               ! Loop over all vertices or element nodes
               ! to reshape shpnodtmp -> shpnod
               !            shglnodtmo -> shglnod
               do vrt = 1, nshl
                  shpnod(:,vrt) = shpnodtmp(nod,vrt)
                  do i = 1, nsd
                     shglnod(:,vrt,i) = shglnodtmp(nod,i,vrt)
                  end do
               end do

               ! Compute real space shape function gradients
               call gradNodalShapeFunctions(shglnodtmp, shgnod)

               ! Stuff


            end do ! </loop nodes>

            ! Deallocate all arrays for next iteration
            deallocate(shpnod, shglnod, shgnod, shpnodtmp, shglnodtmp)
         enddo ! </loop blocks>

      end subroutine slipCorrect
