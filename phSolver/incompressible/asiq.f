        subroutine AsIq (blk, y,       xl, dwl,      shp,
     &                   shgl,    ien,     xmudmi,
     &                   qres,    rmass    )
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c interior elements for the global reconstruction of the diffusive
c flux vector.
c
c input:
c     y     (numnp,ndof)        : Y variables
c     x     (numnp,nsd)         : nodal coordinates
c     shp   (nen,nintg)         : element shape-functions
c     shgl  (nsd,nen,nintg)     : element local shape-function gradients
c     ien   (blk%e)              : nodal connectivity array
c
c output:
c     qres  (numnp,nsd,nsd)  : residual vector for diffusive flux
c     rmass  (numnp)            : lumped mass matrix
c
c----------------------------------------------------------------------
c
      use turbsa      ! access to d2wall & effvisc
      use eblock
      include "common.h"
      type (LocalBlkData) blk

c
        dimension y(nshg,ndof),     !          x(numnp,nsd),            
     &            shp(blk%s,blk%g),         shgl(nsd,blk%s,blk%g),
     &            ien(blk%e,blk%s),      dwl(blk%e,blk%n),
     &            qres(nshg,idflx),    rmass(nshg)
c
        dimension yl(blk%e,blk%s,ndof),          xl(blk%e,blk%n,nsd),
     &            ql(blk%e,blk%s,idflx),  rmassl(blk%e,blk%s),
     &            xmudmi(blk%e,blk%g)
c
        dimension sgn(blk%e,blk%s),       evl(blk%e,blk%s)
c
c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
        do i=1,blk%s
           where ( ien(:,i) < 0 )
              sgn(:,i) = -one
           elsewhere
              sgn(:,i) = one
           endwhere
        enddo

c
c.... gather the variables
c

        call localy(blk,y,      yl,     ien,    ndof,   'gather  ')
!        call localx (blk,x,      xl,     ien,    nsd,    'gather  ')
!        if (iRANS .eq. -2) then ! kay-epsilon
!           call localx (blk,d2wall,   dwl,     ien,    1,     'gather  ')
!        endif

        if ((iDNS.gt.0).and.(itwmod.eq.-2)) then
          call local(blk,effvisc, evl,    ien,    1,      'gather  ')
        endif

c
c.... get the element residuals 
c
        ql     = zero
        rmassl = zero

        call e3q  (blk,yl,         dwl,      shp,      shgl,    
     &             xl,         ql,       rmassl,
     &             xmudmi,     sgn,      evl  )

c
c.... assemble the diffusive flux residual 
c
        call local (blk,qres,   ql,  ien,  idflx,  'scatter ')
        call local (blk,rmass,  rmassl, ien,  1,          'scatter ')
c
c.... end
c
        return
        end


c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c interior elements for the global reconstruction of the diffusive
c flux vector.
c
c----------------------------------------------------------------------
       subroutine AsIqSclr (blk,y,       xl,   dwl,    shp,
     &                       shgl,    ien,     qres,    
     &                       rmass,   cfl,     icflhits   )
c
      use turbsa      ! access to d2wall (no but effvisc does still need this)
      use eblock
      include "common.h"
      type (LocalBlkData) blk

c
        dimension y(nshg,ndof),     !        x(numnp,nsd),            
     &            shp(blk%s,blk%g),         shgl(nsd,blk%s,blk%g),
     &            ien(blk%e,blk%s),      dwl(blk%e,blk%n),
     &            qres(nshg,nsd),           rmass(nshg),
     &            cfl(nshg),           icflhits(nshg)
c
        dimension yl(blk%e,blk%s,ndof),       xl(blk%e,blk%n,nsd),         
     &            ql(blk%e,blk%s,nsd),        rmassl(blk%e,blk%s),
     &            cfll(blk%e,blk%s)
c
        dimension sgn(blk%e,blk%s),       evl(blk%e,blk%s)

        if (blk%o .gt. 1) then
           call getsgn(blk,ien,sgn)
        endif
c
c.... gather the variables
c
        call localy(blk,y,      yl,     ien,    ndof,   'gather  ')
!        call localx (blk,x,      xl,     ien,    nsd,    'gather  ')
!        if (iRANS.eq.-2 .or. iRANS.eq.-5) then ! kay-epsilon and SST
!           call localx (blk,d2wall,   dwl,     ien,    1,     'gather  ')
!        endif

        if ((iDNS.gt.0).and.(itwmod.eq.-2)) then
          call local(blk,effvisc, evl,    ien,    1,      'gather  ')
        endif

c
c.... get the element residuals 
c
        ql     = zero
        rmassl = zero
        cfll = zero

        call e3qSclr  (blk,yl,      dwl,    shp,    shgl,    
     &                 xl,      ql,     rmassl, 
     &                 sgn,     evl,    cfll )

c
c.... assemble the temperature diffusive flux residual 
c
        call local (blk,qres,   ql,  ien,  nsd,  'scatter ')
        call local (blk,rmass,  rmassl, ien,  1, 'scatter ')
c
c.... assemble the CFL values.  cfl will contain the sum of
c     all contributing integration points.  Will divide by
c     the number of contributors to get the average CFL number.
        if (iLSet.eq.2) then
          call localSum (blk, cfl, cfll, ien, icflhits, 1)
        endif
c
c.... end
c
        return
        end

