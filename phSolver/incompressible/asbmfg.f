      subroutine AsBMFG (blk,u,       y,       ac,      xlb,       
     &                   dwl,    shpb,    shglb,
     &                   ienb,    materb,  iBCB,    BCB,
     &                   res,     xKebe)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  boundary elements.
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c Alberto Figueroa, Winter 2004.  CMM-FSI
c Irene Vignon, Spring 2004.
c----------------------------------------------------------------------
c
      use turbSA                ! access to effvisc
      use eblock
      include "common.h"
      type (LocalBlkData) blk

c
        dimension y(nshg,ndofl),           x(numnp,nsd),
     &            ac(nshg,ndofl),          u(nshg,nsd),
     &            shpb(nshl,ngaussb),
     &            shglb(nsd,nshl,ngaussb),         
     &            ienb(npro,nshl),         materb(npro),
     &            iBCB(npro,ndiBCB),       BCB(npro,nshlb,ndBCB),
     &            res(nshg,nflow),        dwl(blk%e,nenl),
     &            evl(blk%e,blk%s)        
c
        dimension yl(blk%e,nshl,ndofl),     xlb(blk%e,nenl,nsd),
     &            rl(blk%e,nshl,nflow),     sgn(npro,nshl),
     &            ul(blk%e,nshl,nsd),       acl(blk%e,nshl,ndofl)
c
!disable        dimension xKebe(npro,9,nshl,nshl) 
     
c
c.... get the matrix of mode signs for the hierarchic basis functions
c
        if (ipord .gt. 1) then
           call getsgn(blk,ienb,sgn)
        endif
c
c.... gather the variables
c
        call localy(blk,y,      yl,     ienb,   ndofl,  'gather  ')
        call localy(blk,ac,     acl,    ienb,   ndofl,  'gather  ')
!        call localx(blk,x,      xlb,    ienb,   nsd,    'gather  ')
        call localx(blk,u,      ul,     ienb,   nsd,    'gather  ')
!        if(iRANS.eq.-2) then
!           call localx(blk,d2wall, dwl, ienb, 1, 'gather  ')
!        endif

        if ((iDNS.gt.0).and.(itwmod.eq.-2)) then
          call local(blk,effvisc, evl, ienb, 1, 'gather  ')
        endif
c
c.... zero the matrices if they are being recalculated
c
!disable       if (lhs. eq. 1)  then
!disable           xKebe = zero
!disable        endif   

c
c.... get the boundary element residuals
c
        rl  = zero
c
c.... 3D
c
        call e3b  (blk,ul,      yl,      acl,     iBCB,    BCB,     
     &             shpb,    shglb,
     &             xlb,     rl,      sgn,     dwl,     xKebe,
     &             evl)
c
c.... assemble the residual and the modified residual
c
        call local (blk,res,    rl,     ienb,   nflow,  'scatter ')

c     
c.... end
c
        return
        end
 

c
c----------------------------------------------------------------------
c
c     This routine computes and assembles the data corresponding to the
c     boundary elements for the temperature equation
c
c----------------------------------------------------------------------
c
      subroutine AsBSclr (blk,y,       xlb,  dwl,     shpb,    shglb,
     &                   ienb,    materb,  iBCB,    BCB,
     &                   res)
      use turbSA ! access to effvisc
      use eblock
      include "common.h"
      type (LocalBlkData) blk

c
        dimension y(nshg,ndofl),           x(numnp,nsd),
     &            shpb(nshl,*),
     &            shglb(nsd,nshl,*),         
     &            ienb(npro,nshl),         materb(npro),
     &            iBCB(npro,ndiBCB),       BCB(npro,nshlb,ndBCB),
     &            res(nshg)         
c
        dimension yl(blk%e,nshl,ndofl),     xlb(blk%e,nenl,nsd),
     &            rl(blk%e,nshl),     sgn(npro,nshl)
        real*8 dwl(blk%e,nshl),          evl(blk%e,blk%s)
c
c.... get the matrix of mode signs for the hierarchic basis functions
c
        if (ipord .gt. 1) then
           call getsgn(blk,ienb,sgn)
        endif
c
c.... gather the variables
c
        call localy(blk,y,      yl,     ienb,   ndofl,  'gather  ')
!        call localx(blk,x,      xlb,    ienb,   nsd,    'gather  ')
!        if(iRANS.eq.-2) then
!           call localx(blk,d2wall, dwl, ienb, 1, 'gather  ')
!        endif

        if ((iDNS.gt.0).and.(itwmod.eq.-2)) then
          call local(blk,effvisc, evl,    ienb,    1,      'gather  ')
        endif
c
c.... get the boundary element residuals
c
        rl  = zero

        call e3bSclr  (blk,yl,      iBCB,    BCB,     shpb,    shglb,
     &                 xlb,     rl,      sgn,     dwl,     evl)
c
c.... assemble the residual and the modified residual
c
        call local (blk,res,    rl,     ienb,   1,  'scatter ')
c     
c.... end
c
        return
        end


