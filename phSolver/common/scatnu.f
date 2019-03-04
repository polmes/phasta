      subroutine scatnu (blk, ien, strl, xmudmi, xnut, shp)
      
      use eblock

      include "common.h"
      type (LocalBlkData) blk

      dimension  ien(blk%e,blk%s),       strl(blk%e,blk%g),
     &           xmudmi(blk%e,blk%g),       shp(blk%s,blk%g)
      dimension  xnut(numnp)

      xmudmi=zero

      if(iLES.eq.5) return  ! Debugging with zero-ed model

      do ii = 1,blk%s
      do intp = 1,blk%g
        xmudmi(:,intp) = xmudmi(:,intp) + xnut(ien(:,ii)) * strl(:,intp)
     &        *shp(ii,intp)
      enddo  
      enddo
c
c  local clipping
c
      rmu=datmat(1,2,1)
      xmudmi=min(xmudmi,1000.0*rmu) !don't let it get larger than 1000 mu
      xmudmi=max(xmudmi, zero) ! don't let (xmudmi) < 0
c      xmudmi=max(xmudmi, -rmu) ! don't let (xmudmi + mu) < 0
      return
      end
