c-----------------------------------------------------------------------
c     create the new statistics arrays
c-----------------------------------------------------------------------
      subroutine initStats(x,   iBC,    iper,   ilwork)
      
      use stats
      include "common.h"

      real*8  x(numnp,3)
      integer ilwork(nlwork), iper(nshg), iBC(nshg)
      
      write(*,*) 'Stats are not developed for compressible code'
      return
      end
