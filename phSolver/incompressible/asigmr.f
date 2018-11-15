        subroutine AsIGMR (blk,  y,       ac,      x,       xmudmi,
     &                     shp,     shgl,    ien,     
     &                     rl,     qres,
     &                     xlhs,   rerrl,   StsVecl,
     &                     cfl,     icflhits )
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
      use eblock
#ifdef HAVE_OMP
      use omp_lib
#endif


      include "common.h"
      type (LocalBlkData) blk

c
        dimension y(nshg,ndofl),              ac(nshg,ndofl),
     &            x(numnp,nsd),              
     &            shp(blk%s,blk%g),            shgl(nsd,blk%s,blk%g),
     &            ien(blk%e,blk%s),
     &            qres(nshg,idflx),           cfl(nshg),
     &            icflhits(nshg)

c
        dimension yl(bsz,blk%s,ndofl),         acl(bsz,blk%s,ndofl),
     &            xl(bsz,blk%n,nsd),           dwl(bsz,blk%n),      
     &            rl(bsz,blk%s,nflow), 
     &            ql(bsz,blk%s,idflx),
     &            evl(bsz,blk%s),
     &            cfll(bsz,blk%s)
c        
        dimension xlhs(bsz,16,blk%s,blk%s)
c
        dimension rlsl(bsz,blk%s,6) 

c
        real*8    StsVecl(bsz,blk%s,nResDims)
        
        dimension xmudmi(blk%e,blk%g)
        dimension sgn(blk%e,blk%s)
c
        real*8 rerrl(bsz,blk%s,6+isurf)
c
c.... gather the variables
c
c
c.... get the matrix of mode signs for the hierarchic basis functions. 
c
        if (blk%o .gt. 1) then
           write(*,*) 'not threaded'
           call getsgn(ien,sgn)
        endif
#ifdef HAVE_OMP_debug
!$OMP critical
        id = omp_get_thread_num ( )
        write (*,*) 'asigmr rank, th_num, iblk, ith ', myrank, id, blk%b, blk%t
!$OMP end critical
#endif

        
        call localy(blk,y,      yl,     ien,    ndofl,  'gather  ')
        call localy(blk,ac,    acl,     ien,    ndofl,  'gather  ')
        call localx(blk,x,      xl,     ien,    nsd,    'gather  ')
        call local (blk,qres,   ql,     ien,    idflx,  'gather  ')
        if (iRANS .eq. -2) then ! kay-epsilon
           call localx (blk,d2wall,   dwl,     ien,    1,     'gather  ')
        endif
 
        if( (iLES.gt.10).and.(iLES.lt.20)) then  ! bardina 
           call local (blk,rls, rlsl,     ien,       6, 'gather  ')  
        else
           rlsl = zero
        endif      
        if ((iDNS.gt.0).and.(itwmod.eq.-2)) then
          call local(blk,effvisc, evl,    ien,    1,      'gather  ')
        endif

c
c.... zero the matrices if they are being recalculated
c
        if (lhs. eq. 1)  then
           xlhs = zero
        endif   
c
c.... get the element residuals, LHS matrix, and preconditioner
c
        rl     = zero
        cfll   = zero

        if(ierrcalc.eq.1) rerrl = zero

        call e3  (blk,yl,      acl,     dwl,     shp,
     &            shgl,    xl,      rl,      
     &            ql,      xlhs, xmudmi, 
     &            sgn,     rerrl,  rlsl,
     &            cfll,  evl)
c
c.... assemble the statistics residual
c
        if ( stsResFlg .eq. 1 ) then
           call e3StsRes (blk, xl, rl, StsVecl )
        endif
c
c.... sum the CFL value from IPs.  These wil be divided by the number of
c     contributors in elmgmr to get average CFL value at node
c
        call localSum (blk,cfl, cfll, ien, icflhits, 1)
c
c.... end
c

        if (exts.and.ires.ne.2) then
           if ((iter.eq.1).and.(mod(lstep,freq).eq.0)) then
              call timeseries(blk,yl,xl,ien,sgn)
           endif
        endif
        
        return
        end



C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c-----------------------------------------------------------------------
c=======================================================================


        subroutine AsIGMRSclr(blk, y,       ac,      x,       
     &                     shp,     shgl,    ien,     
     &                     res,     qres,    xSebe, xmudmi,
     &                     cfl,     icflhits,  cflold )
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  interior elements.
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
      use     turbSA   ! access to d2wall and effvisc
      use eblock
      include "common.h"
      type (LocalBlkData) blk
c
        dimension y(nshg,ndofl),              ac(nshg,ndofl),
     &            x(numnp,nsd),              
     &            shp(blk%s,blk%g),            shgl(nsd,blk%s,blk%g),
     &            ien(blk%e,blk%s),
     &            res(nshg),                  qres(nshg,nsd),
     &            cfl(nshg),                  icflhits(nshg),
     &            cflold(nshg)

c
        real*8    yl(bsz,blk%s,ndofl),        acl(bsz,blk%s,ndofl),
     &            xl(bsz,blk%n,nsd),         
     &            rl(bsz,blk%s),              ql(bsz,blk%s,nsd),
     &            dwl(bsz,blk%n),             evl(bsz,blk%s),
     &            cfll(bsz,blk%s),            cfllold(bsz,blk%s)            
c        
        real*8    xSebe(bsz,blk%s,blk%s),      xmudmi(blk%e,blk%g) 
c
c.... gather the variables
c
        real*8 sgn(blk%e,blk%s)
c
c.... get the matrix of mode signs for the hierarchic basis functions. 
c
        if (blk%o .gt. 1) then
           call getsgn(blk, ien,sgn)
        endif
        
        call localy(blk,y,      yl,     ien,    ndofl,  'gather  ')
        call localy(blk,ac,    acl,     ien,    ndofl,  'gather  ')
        call localx(blk,x,      xl,     ien,    nsd,    'gather  ')
        if(iRANS.lt. 0) 
     &  call localx(blk,d2wall, dwl,    ien,    1,      'gather  ')
        call local (blk,qres,   ql,     ien,    nsd,    'gather  ')

        if ((iDNS.gt.0).and.(itwmod.eq.-2)) then
          call local(blk,effvisc, evl,    ien,    1,      'gather  ')
        endif
        if (iLSet.eq.2) then
          call local(blk,cflold, cfllold, ien,  1,  'gather  ')
        endif
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
      cfll = zero
      call e3Sclr  (blk,yl,      acl,     shp,
     &              shgl,    xl,      dwl,
     &              rl,      ql,      xSebe,   
     &              sgn, xmudmi,  cfll,
     &              cfllold, evl)
c
c.... assemble the residual
c
        call local (blk,res,    rl,     ien,    1,  'scatter ')
c
c.... assemble the CFL values.  cfl will contain the sum of
c     all contributing integration points.  Will divide by
c     the number of contributors to get the average CFL number.
        if (iLSet.eq.2) then
          call localSum (blk,cfl, cfll, ien, icflhits, 1)
        endif
c
c.... end
c
        return
        end
