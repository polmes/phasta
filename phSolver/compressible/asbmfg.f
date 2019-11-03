        subroutine AsBMFG (blk, y,       xlb, dwl,       shpb,    shglb,
     &                     ienb,    materb,  iBCB,    BCB,
     &                     res,     rmes,    EGmass)
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
        dimension y(nshg,ndofl),           xlb(blk%e,blk%n,nsd),
     &            shpb(nshl,ngaussb),      dwl(blk%e,blk%n),
     &            shglb(nsd,nshl,ngaussb),        
     &            ienb(npro,nshl),          materb(npro),
     &            iBCB(npro,ndiBCB),        BCB(npro,nshlb,ndBCB),
     &            res(nshg,nflow),         rmes(nshg,nflow)
c
        dimension ycl(npro,nshl,ndofl),  
     &            rl(npro,nshl,nflow),
     &            rml(npro,nshl,nflow),
     &            EGmass(npro, nshl, nshl) 
c        
        dimension sgn(npro,nshl)
c
c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
        if (ipord .gt. 1) then
           call getsgn(blk, ienb,sgn)
        endif
c
c.... gather the variables
c

        call localy(blk,y,      ycl,     ienb,   ndofl,  'gather  ')
c

        !get the boundary element residuals

        rl  = zero
        rml = zero
c
!  pass the memory location of ycl to both yl and ycl in e3b.  This may
 !  seem dangerous since yl in e3b is :,nflow and ycl is :,ndof but they
 !  do not write to yl (out of bounds at least), only use the data there 
 !  so both will access data
 !  properly from this location.
c
        call e3b  (blk, ycl,     ycl,     iBCB,    BCB,     shpb,    shglb,
     &             xlb,     rl,      rml,     sgn,     EGmass)

        !assemble the residual and the modified residual
        call local(blk,res,    rl,     ienb,   nflow,  'scatter ')
        if (Navier .eq. 1 .and. ires.ne.1 )
     &    call local(blk,rmes,   rml,    ienb,   nflow,  'scatter ')
        
        !end
        return
        end
c
c
c
        subroutine AsBMFGSclr (blk, y,  xlb, dwl,       shpb,    shglb,
     &                         ienb,    materb,  iBCB, 
     &                         BCB,     rest,    rmest)
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
        dimension y(nshg,ndofl),           xl(blk%e,blk%n,nsd),
     &            shpb(nshl,maxsh),        dwl(blk%e,blk%n), 
     &            shglb(nsd,nshl,maxsh),         
     &            ienb(npro,nshl),       materb(npro),
     &            iBCB(npro,ndiBCB),   BCB(npro,nshlb,ndBCB),
     &            rest(nshg),         rmest(nshg)
c
        dimension ycl(npro,nshl,ndofl),   xlb(npro,nenl,nsd),
     &            rtl(npro,nshl),      
     &            rmtl(npro,nshl)
c
c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
        if (ipord .gt. 1) then
           call getsgn(blk, ienb,sgn)
        endif
c
c.... gather the variables
c
        call localy (blk,y,      ycl,     ienb,   ndofl,  'gather  ')
c
c.... get the boundary element residuals
c
        rtl  = zero
        rmtl = zero
c
c.... 3D
c
            call e3bSclr (blk,ycl,    iBCB,    BCB,     
     &                    shpb,  shglb,   sgn,
     &                    xlb,   rtl,     rmtl)
c
c.... assemble the residual and the modified residual
c

        call local (blk,rest,    rtl,     ienb,   1,  'scatter ')


c
        if (Navier .eq. 1)
     &  call local (blk,rmest,   rmtl,    ienb,   1,  'scatter ')
c
c.... end
c
        return
        end


