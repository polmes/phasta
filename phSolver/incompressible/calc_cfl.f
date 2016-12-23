      subroutine calc_cfl (blk, rho,          u1,       u2,
     &                     u3,           dxidx,    rmu,  
     &                     cfl_loc)  
c
c----------------------------------------------------------------------
c
c This routine computes the CFL for each element at an integration point. 
c
c input:
c  u1     (blk%e)           : x1-velocity component
c  u2     (blk%e)           : x2-velocity component
c  u3     (blk%e)           : x3-velocity component
c  dxidx  (blk%e,nsd,nsd)   : inverse of deformation gradient
c
c output:
c  cfl_loc(blk%e) 	   : CFL of the element
c
c----------------------------------------------------------------------
c
       include "common.h"
       include "eblock.h"
       type (LocalBlkData) blk
c
        dimension rho(blk%e),                 u1(blk%e),
     &            u2(blk%e),                  u3(blk%e),
     &            dxidx(blk%e,nsd,nsd), 
     &            rmu(blk%e), 
     &            dt2(blk%e), dt3(blk%e)
c
        dimension gijd(blk%e,6),  rnu(blk%e),  
     &            rhoinv(blk%e), cfl_loc(blk%e)
c
c.... get the metric tensor
c      
      call e3gijd(blk, dxidx, gijd )

      rhoinv=one/rho
      rnu=rmu*rhoinv

       dt2 = ( u1 * ( gijd(:,1) * u1
     4		    + gijd(:,4) * u2
     5		    + gijd(:,6) * u3 )
     6	     + u2 * ( gijd(:,4) * u1
     7		    + gijd(:,2) * u2
     8		    + gijd(:,5) * u3 )
     9	     + u3 * ( gijd(:,6) * u1
     a		    + gijd(:,5) * u2
     1		    + gijd(:,3) * u3 ) ) 
       dt3 = rnu ** 2
     3	          * ( gijd(:,1) ** 2
     4	            + gijd(:,2) ** 2
     5		    + gijd(:,3) ** 2
     6		    + 2.
     7		  * ( gijd(:,4) ** 2
     8		    + gijd(:,5) ** 2
     9		    + gijd(:,6) ** 2 ) )
c
         if (ires == 1) then
            cfl_loc= cfl_loc+sqrt(max(dt2,dt3/two))/(Dtgl*two)
         endif
c     
c.... return
c
        return
        end

