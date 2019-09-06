      subroutine SSTStress(blk,rmu, g1yi, g2yi, g3yi, Sclr, TauSST)
c------------------------------------------------------------------------
c
c This routine computes the total viscous stresses (viscous plus Reynolds)
c according to the 2003 version of the Menter paper.
c The SST model is based on the Boussinesq approximation of the Reynolds
c stresses, so ReyStress = 2*mut*S_ij, making the total viscous stresses
c in the RANS momentum equation T_ij = 2*(nu+mut)*S_ij.
c The SST model adds anisotropy to the Reynolds stresses by saying that
c the Reynolds stresses are now
c     ReyStress_ij= 2*mut*S_ij - 0.66666*(d(u_k)/d(x_k))D_ij
c Where D_ij is the Kronecker Delta
c Inputs: rmu - nu+mut
c         g1yi - x gradient of the solution vector y
c         g2yi - y gradient of the solution vector y
c         g3yi - z gradient of the solution vector y
c         Sclr - scalar 1 (k) value for SST
c Outputs : TauSST - the 6 unique total viscous stresses (viscous+Reynolds)
c
c-------------------------------------------------------------------------
      use eblock
      include "common.h"
      type (LocalBlkData) blk


      dimension g1yi(blk%e,nflow), g2yi(blk%e,nflow), g3yi(blk%e,nflow),
     &          rmu(blk%e), TauSST(blk%e,6), Sclr(blk%e)
      real*8 Sij(blk%e,6), TauSSTA(blk%e,6), DivU(blk%e),
     &       mut(blk%e), rho(blk%e)

      
      rho = datmat(1,1,1)

c ... Compute the eddy viscosity
      do ii=1,blk%e
         mut(ii) = rmu(ii) - datmat(1,2,1)
      enddo

c ... Compute S_ij=1/2*(ui,j+uj,i)
      Sij(:,1) = g1yi(:,2) ! S_11
      Sij(:,2) = g2yi(:,3) ! S_22
      Sij(:,3) = g3yi(:,4) ! S_33
      Sij(:,4) = 0.50*(g2yi(:,2) + g1yi(:,3)) ! S_12
      Sij(:,5) = 0.50*(g3yi(:,3) + g2yi(:,4)) ! S_23
      Sij(:,6) = 0.50*(g3yi(:,2) + g1yi(:,4)) ! S_13

c ..  Compute the divergence of velocity

      if (matflg(1,1) .eq. 0) then !compressible
        DivU(:) = g1yi(:,2) + g2yi(:,3) + g3yi(:,4)
      else                         !incompressible
        DivU(:) = zero
      endif

c ... Compute the SST total (viscous+Reynolds) stresses

c ... Calculate the anisotropic addition   
c  | 1   4   6  |
c  |            |
c  | 4   2   5  |
c  |            |
c  | 6   5   3  |
      TauSSTA = zero
      TauSSTA(:,1) = mut * DivU(:) + rho(:) * Sclr(:)
      TauSSTA(:,2) = mut * DivU(:) + rho(:) * Sclr(:)
      TauSSTA(:,3) = mut * DivU(:) + rho(:) * Sclr(:)

      do ii=1,6
        TauSST(:,ii) = two*rmu(:)*Sij(:,ii) - two /three * TauSSTA(:,ii)
      enddo
      
      
      return
      end

