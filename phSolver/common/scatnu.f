      subroutine scatnu (blk, ien, strl, xmudmi, shp)
      
      use eblock
      use lesArrs

      include "common.h"
      type (LocalBlkData) blk

      dimension  ien(blk%e,blk%s),       strl(blk%e,blk%g),
     &           xmudmi(blk%e,blk%g),       shp(blk%s,blk%g)
      dimension  xnut(numnp)
      integer iel, inum
      real*8 nquadsInv, nquads, tmp(blk%e), whist, wcur

      xmudmi=zero
      xnut = cdelsq(:,1)
      iel = blk%i
      inum = iel+blk%e-1
      nquads = blk%g
      nquadsInv = one/nquads

c     Get time average weights
      if (irunTave.eq.0) then
        denom=max(1.0d0*(lstep),one)
        if(dtavei.lt.0) then
          wcur=one/denom
        else
          wcur=dtavei
        endif  
        whist=1.0-wcur
      else ! if irunTave=1 then doing running time average
        denom = max(1,lstep-irunTaveSt)
        wcur = one/denom
        whist = one-wcur  
      endif

      if(iLES.eq.5) return  ! Debugging with zero-ed model

      do ii = 1,blk%s
        do intp = 1,blk%g
          xmudmi(1:blk%e,intp) = xmudmi(1:blk%e,intp) + 
     &                         xnut(ien(1:blk%e,ii)) * strl(1:blk%e,intp)
     &                         * shp(ii,intp)
        enddo  
      enddo
c
c  local clipping
c
      rmu=datmat(1,2,1)
      xmudmi=min(xmudmi,1000.0*rmu) !don't let it get larger than 1000 mu
      xmudmi=max(xmudmi, zero) ! don't let (xmudmi) < 0
c      xmudmi=max(xmudmi, -rmu) ! don't let (xmudmi + mu) < 0

c     Compute average nu_t over element and time average to print to restart
      tmp = zero
      do intp=1,blk%g
         tmp = tmp+xmudmi(1:blk%e,intp)
      enddo
      lesnut(iel:inum,1) = tmp*nquadsInv
      lesnut(iel:inum,2) = whist*lesnut(iel:inum,2) + 
     &                     wcur*lesnut(iel:inum,1)

      return
      end
