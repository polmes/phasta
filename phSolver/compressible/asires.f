        subroutine AsIRes (blk, y,      yc,     xl, dwl,      xmudmi,    
     &                     shp,    shgl,   ien,    mater,
     &                     rmes,    ac)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  interior elements.
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        use eblock
        use rlssave     ! Use the resolved Leonard stresses at the nodes.
        include "common.h"
      type (LocalBlkData) blk
c
        dimension y(nshg,nflow),            yc(nshg,ndofl),
     &            ac(nshg,nflow),
     &            shp(nshl,ngauss),  
     &            shgl(nsd,nshl,ngauss),
     &            ien(npro,nshl),       mater(npro),
     &            rmes(nshg,nflow)
c
        dimension yl(npro,nshl,nflow),       ycl(npro,nshl,ndofl),
     &            xl(npro,nenl,nsd),         acl(npro,nshl,nflow),
     &           rml(npro,nshl,nflow), ql(npro,nshl,(nflow-1)*nsd),
     &            dwl(blk%e,blk%n)
c        
        dimension  xmudmi(npro,ngauss)        
        dimension sgn(npro,nshl)

        dimension rlsl(npro,nshl,6) 
c
        real*8 rerrl(npro,nshl,6)
c
c
c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
c     this is done in hierarchic.f
c$$$        do i=1,nshape
c$$$           where ( ien(:,i) < 0 )
c$$$              sgn(:,i) = -one
c$$$           elsewhere
c$$$              sgn(:,i) = one
c$$$           endwhere
c$$$        enddo
c
        if (ipord .gt. 1) then
           call getsgn(blk, ien,sgn)
        endif
c
c.... gather the variables
c
        call localy(blk,y,      yl,     ien,    nflow,  'gather  ')
c
        call localy(blk,yc,     ycl,    ien,    ndofl,  'gather  ')
        call localy(blk,ac,     acl,    ien,    nflow,  'gather  ')
       

        if( (iLES.gt.10).and.(iLES.lt.20)) then ! bardina 

           call local (blk,rls, rlsl,     ien,       6, 'gather  ')  
        else
           rlsl = zero
        endif

c
c.... get the element residual
c

        rml = zero
        
        EGmassd= one  ! just a dummy real since we don't have a LHS with MFI
        if(ierrcalc.eq.1) rerrl = zero        
        ttim(31) = ttim(31) - secs(0.0)
c	write(*,*) 'calling e3'

            call e3  (blk,yl,      ycl,     acl,     shp,
     &                shgl,    xl,      rml,     rml,
     &                xmudmi,  BDiagl,  ql,      sgn, rlsl, EGmassd,
     &                rerrl)

        ttim(31) = ttim(31) + secs(0.0)
c
c.... assemble the modified residual
c
        if (iabres .eq. 1) rml = abs(rml)
c
        call local (blk,rmes,   rml,    ien,    nflow,  'scatter ')
c
c.... end
c
        return
        end
