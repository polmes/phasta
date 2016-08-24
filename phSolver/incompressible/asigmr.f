        subroutine AsIGMR (iblk, blk, y,       ac,      x,       xmudmi,
     &                     shp,     shgl,    ien,     
     &                     rl,     qres,
     &                     xKebe,   xGoC,    rerrl)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  interior elements.
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      use stats
      use rlssave  ! Use the resolved Leonard stresses at the nodes.
      use timedata    ! time series
      use turbsa                ! access to d2wall


      include "common.h"
c
        dimension y(nshg,ndofl),              ac(nshg,ndofl),
     &            x(numnp,nsd),              
     &            shp(blk%s,blk%g),            shgl(nsd,blk%s,blk%g),
     &            ien(ibksiz,blk%s),
     &            qres(nshg,idflx)

c
        dimension yl(ibksiz,blk%s,ndofl),         acl(ibksiz,blk%s,ndofl),
     &            xl(ibksiz,nenl,nsd),           dwl(ibksiz,nenl),      
     &            rl(ibksiz,blk%s,nflow), 
     &            ql(ibksiz,blk%s,idflx)
c        
        dimension xKebe(ibksiz,9,blk%s,blk%s), 
     &            xGoC(ibksiz,4,blk%s,blk%s)
c
        dimension rlsl(ibksiz,blk%s,6) 

c
        real*8    StsVecl(ibksiz,blk%s,nResDims)
        
        dimension xmudmi(ibksiz,blk%g)
        dimension sgn(ibksiz,blk%s)
        dimension CFLworst(ibksiz)
c
        real*8 rerrl(ibksiz,blk%s,6)
c
c.... gather the variables
c
c
c.... get the matrix of mode signs for the hierarchic basis functions. 
c
        if (ipord .gt. 1) then
           call getsgn(ien,sgn)
        endif
        
        call localy(y,      yl,     ien,    ndofl,  'gather  ')
        call localy(ac,    acl,     ien,    ndofl,  'gather  ')
        call localx(x,      xl,     ien,    nsd,    'gather  ')
        call local (qres,   ql,     ien,    idflx,  'gather  ')
        if (iRANS .eq. -2) then ! kay-epsilon
           call localx (d2wall,   dwl,     ien,    1,     'gather  ')
        endif
 
        if( (iLES.gt.10).and.(iLES.lt.20)) then  ! bardina 
           call local (rls, rlsl,     ien,       6, 'gather  ')  
        else
           rlsl = zero
        endif      

c
c.... zero the matrices if they are being recalculated
c
        if (lhs. eq. 1)  then
           xKebe = zero
           xGoC  = zero
        endif   
c
c.... get the element residuals, LHS matrix, and preconditioner
c
        rl     = zero

        if(ierrcalc.eq.1) rerrl = zero

        call e3  (yl,      acl,     dwl,     shp,
     &            shgl,    xl,      rl,      
     &            ql,      xKebe,   xGoC,    xmudmi, 
     &            sgn,     rerrl,  rlsl,     CFLworst)
c
c.... assemble the statistics residual
c
        if ( stsResFlg .eq. 1 ) then
           call e3StsRes ( xl, rl, StsVecl )
        endif
c
c.... end
c

        if (exts.and.ires.ne.2) then
           if ((iter.eq.1).and.(mod(lstep,freq).eq.0)) then
              call timeseries(yl,xl,ien,sgn)
           endif
        endif
        
        return
        end



C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c-----------------------------------------------------------------------
c=======================================================================


        subroutine AsIGMRSclr(y,       ac,      x,       
     &                     shp,     shgl,    ien,     
     &                     res,     qres,    xSebe, xmudmi )
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  interior elements.
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      use     turbSA  
      include "common.h"
c
        dimension y(nshg,ndofl),              ac(nshg,ndofl),
     &            x(numnp,nsd),              
     &            shp(blk%s,blk%g),            shgl(nsd,blk%s,blk%g),
     &            ien(ibksiz,blk%s),
     &            res(nshg),                  qres(nshg,nsd)

c
        real*8    yl(ibksiz,blk%s,ndofl),        acl(ibksiz,blk%s,ndofl),
     &            xl(ibksiz,nenl,nsd),         
     &            rl(ibksiz,blk%s),              ql(ibksiz,blk%s,nsd),
     &            dwl(ibksiz,nenl)            
c        
        real*8    xSebe(ibksiz,blk%s,blk%s),      xmudmi(ibksiz,blk%g) 
c
c.... gather the variables
c
        real*8 sgn(ibksiz,blk%s)
c
c.... get the matrix of mode signs for the hierarchic basis functions. 
c
        if (ipord .gt. 1) then
           call getsgn(ien,sgn)
        endif
        
        call localy(y,      yl,     ien,    ndofl,  'gather  ')
        call localy(ac,    acl,     ien,    ndofl,  'gather  ')
        call localx(x,      xl,     ien,    nsd,    'gather  ')
        if(iRANS.lt. 0) 
     &  call localx(d2wall, dwl,    ien,    1,      'gather  ')
        call local (qres,   ql,     ien,    nsd,    'gather  ')
c
c.... zero the matrices if they are being recalculated
c
        if (lhs. eq. 1)  then
           xSebe = zero
        endif   
c
c.... get the element residuals, LHS matrix, and preconditioner
c
      rl = zero
      call e3Sclr  (yl,      acl,     shp,
     &              shgl,    xl,      dwl,
     &              rl,      ql,      xSebe,   
     &              sgn, xmudmi)
c
c.... assemble the residual
c
        call local (res,    rl,     ien,    1,  'scatter ')
c
c.... end
c
        return
        end
