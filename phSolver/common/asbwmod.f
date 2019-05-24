      subroutine asbwmod (y,      ac,      x,      BC,   iBC,
     &                    iper,   ilwork,  ifath)
c
c----------------------------------------------------------------------
c
c This routine assembles traction BCs for a modeled wall
c
c----------------------------------------------------------------------
c
      use pointer_data
      use spanStats
      include "common.h"
c
      dimension y(nshg,ndof),         x(numnp, nsd),
     &          BC(nshg,ndofBC),      iBC(nshg),
     &          iper(nshg),           ilwork(nlwork),
     &          ifath(numnp),
     &          ac(nshg,ndof)
c
c.... compute and assemble the residuals corresponding to the 
c     boundary integral
c
              call settauw (y, x, BC, ifath)
c
c.... enforce the new BC for SA variable
c
           isclr = 1
           if (iRANS.eq.-1) then ! S-A RANS
              call itrBCSclr (y, ac, iBC,  BC, iper, ilwork)
           endif
c
           return
           end
