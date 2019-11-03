        subroutine AsBRes (blk, y, yc,      xlb,dwl,       shpb,    shglb,
     &                     ienb,    materb,  iBCB,    BCB,
     &                     rmes)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  boundary elements.
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        use eblock
        include "common.h"
        type (LocalBlkData) blk
c
        dimension y(nshg,nflow),           dwl(blk%e,blk%n),           
     &            yc(nshg,ndof),           shpb(nshl,ngaussb),
     &            shglb(nsd,nshl,ngaussb),
     &            ienb(npro,nshl),        materb(npro),
     &            iBCB(npro,ndiBCB),     BCB(npro,nshlb,ndBCB),
     &            rmes(nshg,nflow)
c
        dimension yl(npro,nshl,nflow),    xlb(npro,nenl,nsd),
     &            ycl(npro,nshl,ndof),    rml(npro,nshl,nflow)
c
        dimension sgn(npro,nshl)
c
c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
        if (ipord .gt. 1) then
           call getsgn(ienb,sgn)
        endif
c     
c.... gather the variables
c
        call localy(blk, y,      yl,     ienb,   nflow,  'gather  ')
        call localy(blk, yc,     ycl,    ienb,   ndof,  'gather  ')
c
c.... get the boundary element residuals
c
        rml = zero
        call e3b  (blk, yl,      ycl,     iBCB,    BCB,     shpb,    shglb,
     &             xlb,     rml,     rml,     sgn)

c
c.... assemble the residual and the modified residual
c
        if (iabres .eq. 1) rml = abs(rml)
c
        call local(blk, rmes,   rml,    ienb,   nflow,  'scatter ')
c
c.... end
c
        return
        end




