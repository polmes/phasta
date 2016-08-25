c-----------------------------------------------------------------------
c
c  compute the metrics of the mapping from global to local 
c  coordinates and the jacobian of the mapping (weighted by 
c  the quadrature weight
c
c-----------------------------------------------------------------------
      subroutine e3metric(blk,ith,  xl,      shgl,     dxidx,
     &                      shg,     WdetJ)

      include "common.h"
      include "eblock.h"
      type (LocalBlkData) blk

      
      real*8     xl(blk%e,blk%n,nsd),    shgl(blk%e,nsd,blk%s),
     &           dxidx(blk%e,nsd,nsd),  shg(blk%e,blk%s,nsd), 
     &           WdetJ(blk%e)

      real*8     dxdxi(blk%e,nsd,nsd),  tmp(blk%e)

c
c.... compute the deformation gradient
c
      dxdxi = zero
c
       do n = 1, blk%n
          dxdxi(1:blk%eblk%e,1,1) = dxdxi(1:blk%eblk%e,1,1) + xl(1:blk%eblk%e,n,1) * shgl(1:blk%eblk%e,1,n)
          dxdxi(1:blk%eblk%e,1,2) = dxdxi(1:blk%eblk%e,1,2) + xl(1:blk%eblk%e,n,1) * shgl(1:blk%eblk%e,2,n)
          dxdxi(1:blk%eblk%e,1,3) = dxdxi(1:blk%eblk%e,1,3) + xl(1:blk%eblk%e,n,1) * shgl(1:blk%eblk%e,3,n)
          dxdxi(1:blk%eblk%e,2,1) = dxdxi(1:blk%eblk%e,2,1) + xl(1:blk%eblk%e,n,2) * shgl(1:blk%eblk%e,1,n)
          dxdxi(1:blk%eblk%e,2,2) = dxdxi(1:blk%eblk%e,2,2) + xl(1:blk%eblk%e,n,2) * shgl(1:blk%eblk%e,2,n)
          dxdxi(1:blk%eblk%e,2,3) = dxdxi(1:blk%eblk%e,2,3) + xl(1:blk%eblk%e,n,2) * shgl(1:blk%eblk%e,3,n)
          dxdxi(1:blk%eblk%e,3,1) = dxdxi(1:blk%eblk%e,3,1) + xl(1:blk%eblk%e,n,3) * shgl(1:blk%eblk%e,1,n)
          dxdxi(1:blk%eblk%e,3,2) = dxdxi(1:blk%eblk%e,3,2) + xl(1:blk%eblk%e,n,3) * shgl(1:blk%eblk%e,2,n)
          dxdxi(1:blk%eblk%e,3,3) = dxdxi(1:blk%eblk%e,3,3) + xl(1:blk%eblk%e,n,3) * shgl(1:blk%eblk%e,3,n)
       enddo
c
c.... compute the inverse of deformation gradient
c
       dxidx(:,1,1) =   dxdxi(:,2,2) * dxdxi(:,3,3) 
     &                - dxdxi(:,3,2) * dxdxi(:,2,3)
       dxidx(:,1,2) =   dxdxi(:,3,2) * dxdxi(:,1,3) 
     &                - dxdxi(:,1,2) * dxdxi(:,3,3)
       dxidx(:,1,3) =  dxdxi(:,1,2) * dxdxi(:,2,3) 
     &                - dxdxi(:,1,3) * dxdxi(:,2,2)
       tmp          = one / ( dxidx(:,1,1) * dxdxi(:,1,1) 
     &                       + dxidx(:,1,2) * dxdxi(:,2,1)  
     &                       + dxidx(:,1,3) * dxdxi(:,3,1) )
       dxidx(:,1,1) = dxidx(:,1,1) * tmp
       dxidx(:,1,2) = dxidx(:,1,2) * tmp
       dxidx(:,1,3) = dxidx(:,1,3) * tmp
       dxidx(:,2,1) = (dxdxi(:,2,3) * dxdxi(:,3,1) 
     &                - dxdxi(:,2,1) * dxdxi(:,3,3)) * tmp
       dxidx(:,2,2) = (dxdxi(:,1,1) * dxdxi(:,3,3) 
     &                - dxdxi(:,3,1) * dxdxi(:,1,3)) * tmp
       dxidx(:,2,3) = (dxdxi(:,2,1) * dxdxi(:,1,3) 
     &                - dxdxi(:,1,1) * dxdxi(:,2,3)) * tmp
       dxidx(:,3,1) = (dxdxi(:,2,1) * dxdxi(:,3,2) 
     &                - dxdxi(:,2,2) * dxdxi(:,3,1)) * tmp
       dxidx(:,3,2) = (dxdxi(:,3,1) * dxdxi(:,1,2) 
     &                - dxdxi(:,1,1) * dxdxi(:,3,2)) * tmp
       dxidx(:,3,3) = (dxdxi(:,1,1) * dxdxi(:,2,2) 
     &                - dxdxi(:,1,2) * dxdxi(:,2,1)) * tmp
c
       WdetJ = Qwt(blk%l,ith) / tmp
c
c.... compute the global gradient of shape-functions
c
       do n = 1, blk%s
          shg(:,n,1) = shgl(:,1,n) * dxidx(:,1,1) + 
     &                 shgl(:,2,n) * dxidx(:,2,1) +
     &                 shgl(:,3,n) * dxidx(:,3,1)
          shg(:,n,2) = shgl(:,1,n) * dxidx(:,1,2) + 
     &                 shgl(:,2,n) * dxidx(:,2,2) +
     &                 shgl(:,3,n) * dxidx(:,3,2) 
          shg(:,n,3) = shgl(:,1,n) * dxidx(:,1,3) + 
     &                 shgl(:,2,n) * dxidx(:,2,3) +
     &                 shgl(:,3,n) * dxidx(:,3,3) 
       enddo

       return
       end



