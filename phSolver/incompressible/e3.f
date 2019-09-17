        subroutine e3 (blk, yl,      acl,     dwl,     shp,
     &                 shgl,    xl,      rl,      ql,
     &                 xlhs, xmudmi,  sgn, 
     &                 rerrl, rlsl,    cfll, evl)
c                                                                      
c----------------------------------------------------------------------
c
c     This routine calculates the residual and tangent matrix for the 
c     UBar formulation of the incompressible Navier Stokes equations.
c
c
c input:    e    a   1..5   when we think of U^e_a  and U is 5 variables
c  yl     (bsz,blk%s,ndof)      : Y variables (not U)
c  acl    (bsz,blk%s,ndof)      : Y acceleration (Y_{,t})
c  shp    (nen,blk%g)           : element shape-functions  N_a
c  shgl   (nsd,nen,blk%g)       : element local-shape-functions N_{a,xi}
c  wght   (blk%g)               : element weight (for quadrature)
c  xl     (bsz,blk%n,nsd)       : nodal coordinates at current step (x^e_a)
c  ql     (bsz,blk%s,nsd*nsd) : diffusive flux vector (don't worry)
c  rlsl   (bsz,blk%s,6)       : resolved Leonard stresses
c  evl    (bsz,blk%s)         : effective viscosity for turb. wall function
c
c output:
c  rl     (bsz,blk%s,nflow)      : element RHS residual    (G^e_a)
c  rml    (bsz,blk%s,nflow)      : element modified residual  (G^e_a tilde)
c  xlhs  (bsz,16,blk%s,blk%s)  : element LHS tangent mass matrix
c  cfll   (bsz,blk%s) 		: CFL of the element
c
c Note: This routine will calculate the element matrices for the
c        Hulbert's generalized alpha method integrator
c
c Mathematics done by:  Michel Mallet, Farzin Shakib (version 1)
c                       Farzin Shakib                (version 2)
c
c K. E. Jansen,   Winter 1999.   (advective form formulation)
c C. H. Whiting,  Winter 1999.   (advective form formulation)
c----------------------------------------------------------------------
c
! only needed if debugging      use omp_lib
      use spat_var_eps   ! use spatially-varying epl_ls
      use eblock
      include "common.h"
      type (LocalBlkData) blk


c
        dimension yl(bsz,blk%s,ndof),
     &            acl(bsz,blk%s,ndof),       
     &            shp(blk%s,blk%g),       shgl(nsd,blk%s,blk%g),
     &            xl(bsz,blk%n,nsd),      dwl(bsz,blk%n),
     &            rl(bsz,blk%s,nflow),     ql(bsz,blk%s,idflx),
     &            cfll(bsz,blk%s),        evl(bsz,blk%s)
c      
        real*8, allocatable, dimension(:,:,:,:,:) :: xK_qp
        real*8, allocatable, dimension(:,:,:,:) :: rl_qp,rerrl_qp

#ifdef SP_LHS
        real*4 xlhs(bsz,16,blk%s,blk%s)
#else
        real*8 xlhs(bsz,16,blk%s,blk%s)
#endif
c
c.... local declarations
c
        dimension g1yi(blk%e,ndof),        g2yi(blk%e,ndof),
     &            g3yi(blk%e,ndof),        shg(blk%e,blk%s,nsd),
     &            aci(blk%e,3),            dxidx(blk%e,nsd,nsd),       
     &            WdetJ(blk%e),            rho(blk%e),
     &            pres(blk%e),             u1(blk%e),
     &            u2(blk%e),               u3(blk%e),
     &            rLui(blk%e,nsd),         uBar(blk%e,nsd),
     &            xmudmi(blk%e,blk%g),     sgn(blk%e,blk%s), 
     &            shpfun(blk%e,blk%s),      shdrv(blk%e,nsd,blk%s),
     &            rmu(blk%e),              tauC(blk%e),
     &            tauM(blk%e),             tauBar(blk%e),
     &            src(blk%e,3),            gyti(blk%e,nsd),
     &            Sclr(blk%e),             mut(blk%e),
     &            ssq(blk%e)

        dimension rlsl(bsz,blk%s,6),      rlsli(blk%e,6)

        real*8    rerrl(bsz,blk%s,6+isurf)
        integer   aa

c
c... needed for CFL calculation
c
      real*8 cfll_loc(blk%e)

c
c     
c.... local reconstruction of diffusive flux vector for quadratics
c     or greater but NOT for bflux since local mass was not mapped
c
        if ( idiff==2 .and. ires .eq. 1 ) then
           call e3ql (blk, yl,        dwl,       shp,       shgl, 
     &                xl,        ql,        xmudmi, 
     &                sgn,       evl)
        endif
c
c.... loop through the integration points
c
! natural place for rdelta = TMRC() but moved lower to time just loop 
#ifdef HAVE_OMP_QP  
	allocate( rl_qp(bsz,blk%s,nflow,blk%g))
	allocate( rerrl_qp(bsz,blk%s,6+isurf,blk%g))
        allocate( xK_qp(bsz,16,blk%s,blk%s,blk%g))
        rl_qp=zero
        rerrl_qp=zero
        xK_qp=zero
#endif
! time just loop 
!       rdelta = TMRC() 
#ifdef HAVE_OMP_QP
!$OMP  parallel do 
!$OMP& private (ith,sgn,shpfun,shdrv,rmu,rho,aci,g1yi,g2yi,g3yi)
!$OMP& private (shg,dxidx,WdetJ,pres,u1,u2,u3,rLui,src,rlsi)
!$OMP& private (tauC,tauM,tauBar,uBar,id)
!$OMP& shared (shp,shgl,dwbszl,yl,xmudmi,xl,acl,rlsl)
!$OMP& shared (rl_qp,xK_qp)
#endif
        do ith = 1, blk%g
        if (Qwt(lcsyst,ith) .eq. zero) cycle          ! precaution
c
c.... get the hierarchic shape functions at this int point
c
        call getshp(blk,ith, shp,          shgl,      sgn, 
     &              shpfun,       shdrv)
c
c.... calculate the integration variables
c
        call e3ivar (blk, ith,yl,          acl,       shpfun,
     &               shdrv,       xl,
     &               aci,         g1yi,      g2yi,    
     &               g3yi,        shg,       dxidx,   
     &               WdetJ,       rho,       pres, 
     &               u1,          u2,        u3,              
     &               ql,          rLui,      src,
#ifdef HAVE_OMP_QP
     &               rerrl_qp(:,:,:,ith),      
#else
     &               rerrl,
#endif
     &               rlsl,      rlsli,
     &               dwl) 
c
c.... get necessary fluid properties (including eddy viscosity)
c
        if(iRANS.eq.-5) then ! solving K-W model
          ssq(:)    = sqrt( 2.0d0 * (g1yi(:,2) ** 2
     &                + g2yi(:,3) ** 2
     &                + g3yi(:,4) ** 2
     &                + 0.5d0
     &                * ( (g3yi(:,3) + g2yi(:,4)) ** 2
     &                + (g1yi(:,4) + g3yi(:,2)) ** 2
     &                + (g2yi(:,2) + g1yi(:,3)) ** 2 )))
          call AddEddyViscKomgq(blk,ith,yl, dwl, shdrv, shpfun, rho, 
     &                          rmu, ssq, xl)
        else
          call getdiff(blk,ith, dwl,  yl,     shpfun,     xmudmi, xl,   rmu, rho,
     &               elem_local_size(blk%i),
     &               evl)
        endif
c
c.... compute CFL number
c
        cfll_loc = zero
        call calc_cfl(blk, rho,          u1,       u2,
     &                u3,           dxidx,    rmu,    
     &                cfll_loc)

        do i=1,blk%s
          cfll(1:blk%e,i) = cfll(1:blk%e,i) + shpfun(:,i)*cfll_loc
        enddo
c
c.... compute the stabilization terms
c
        call e3stab (blk, rho,          u1,       u2,
     &               u3,           dxidx,    rLui,   
     &               rmu,          tauC,     tauM,   
     &               tauBar,       uBar )  
c
c.... compute the residual contribution at this integration point
c
        call e3Res ( blk, u1,        u2,         u3,
     &               uBar,      aci,        WdetJ,
     &               g1yi,      g2yi,       g3yi,
     &               rLui,      rmu,        rho,
     &               tauC,      tauM,       tauBar,
     &               shpfun,    shg,        src,
#ifdef HAVE_OMP_QP
     &               rl_qp(:,:,:,ith),      
#else
     &               rl,
#endif
     &               pres,
     &               rlsli)
c
c.... compute the tangent matrix contribution
c
        if (lhs .eq. 1) then
           call e3LHS ( blk,u1,        u2,         u3,
     &                  uBar,      WdetJ,      rho,
     &                  rLui,      rmu,
     &                  tauC,      tauM,       tauBar,
     &                  shpfun,    shg,        
#ifdef HAVE_OMP_QP
     &                  xK_qp(:,:,:,:,ith),
#else
     &                  xlhs)
#endif
        endif

c
c.... end of integration loop
c
      enddo
!just loop for now
!      rdelta = TMRC() - rdelta
!      rthreads = rthreads + rdelta
!just
#ifdef HAVE_OMP_QP
c
c here we accumulate the thread work from each quadrature point
c
      do i=1,blk%g
       rl(1:blk%e,1:blk%s,1:nflow)=rl(1:blk%e,1:blk%s,1:nflow)
     &                       +rl_qp(1:blk%e,1:blk%s,1:nflow,i)
       rerrl(1:blk%e,1:blk%s,1:6+isurf)=rerrl(1:blk%e,1:blk%s,1:isurf+1)
     &                       +rerrl_qp(1:blk%e,1:blk%s,1:isurf+1,i)
       xlhs(1:blk%e,1:12,1:blk%s,1:blk%s)=xlhs(1:blk%e,1:12,1:blk%s,1:blk%s)
     &                                +xK_qp(1:blk%e,1:12,1:blk%s,1:blk%s,i)
       xlhs(1:blk%e,16,1:blk%s,1:blk%s)=xlhs(1:blk%e,16,1:blk%s,1:blk%s)
     &                                +xK_qp(1:blk%e,16,1:blk%s,1:blk%s,i)
      enddo

      deallocate(xK_qp)
      deallocate(rl_qp)
#endif
!natural place for the rthreads timer when we want to capture all but moved to just loop for now
c
c.... symmetrize C and G from -G^T since only G was computed in e3lhs
      if (lhs .eq. 1) then
         do ib = 1, blk%s
            do iaa = 1, ib-1
               xlhs(1:blk%e,16,iaa,ib) = xlhs(1:blk%e,16,ib,iaa)
            enddo
            do iaa = 1,blk%s
               xlhs(1:blk%e,13:15,iaa,ib) = -xlhs(1:blk%e,4:12:4,ib,iaa)
            enddo
         enddo
      endif
c
c.... return
c
      return
      end


c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c###################################################################


      subroutine e3Sclr (blk,yl,      acl,     shp,
     &                     shgl,    xl,      dwl,
     &                     rl,      ql,      xSebe,   
     &                     sgn,     xmudmi,  cfll,
     &                   cfllold, evl, gradVl, IDDESfunl )
c                                                                      
c----------------------------------------------------------------------
c
c     This routine calculates the residual and tangent matrix for the 
c     advection - diffusion equation for scalar.
c
c K. E. Jansen,   Winter 1999.   (advective form formulation)
c C. H. Whiting,  Winter 1999.   (advective form formulation)
c----------------------------------------------------------------------
c
      use eblock
      include "common.h"
      type (LocalBlkData) blk

c
      real*8    yl(bsz,blk%s,ndof),     acl(bsz,blk%s,ndof),       
     &            shp(blk%s,blk%g),       shgl(nsd,blk%s,blk%g),
     &            xl(bsz,blk%n,nsd),      rl(bsz,blk%s),          
     &            ql(bsz,blk%s,nsd),     
     &            dwl(bsz,blk%n),         cfll(bsz,blk%s),
     &          cfllold(bsz,blk%s),     evl(bsz,blk%s),
     &          gradVl(bsz,blk%s,nsdsq)
#ifdef SP_LHS
      real*4    xSebe(bsz,blk%s,blk%s)
#else
      real*8    xSebe(bsz,blk%s,blk%s)
#endif
c
c.... local declarations
c
      real*8    gradS(blk%e,nsd),        shg(blk%e,blk%s,nsd),
     &            Sdot(blk%e),             Sclr(blk%e),
     &            dxidx(blk%e,nsd,nsd),    WdetJ(blk%e),      
     &            u1(blk%e),     u2(blk%e), u3(blk%e),
     &            sgn(blk%e,blk%s),         shpfun(blk%e,blk%s),       
     &            shdrv(blk%e,nsd,blk%s),   rLS(blk%e),
     &            tauS(blk%e),             diffus(blk%e),
     &            srcL(blk%e),             srcR(blk%e),
     &            gGradS(blk%e,nsd),       dcFct(blk%e),
     &            giju(blk%e,6),           mut(blk%e),
     &            F1(blk%e),               F2(blk%e)
c
c.... Source terms sometimes take the form (beta_i)*(phi,_i).  Since
c     the convective term has (u_i)*(phi,_i), it is useful to treat
c     beta_i as a "correction" to the velocity.  In calculating the
c     stabilization terms, the new "modified" velocity (u_i-beta_i) is 
c     then used in place of the pure velocity for stabilization terms,
c     and the source term sneaks into the RHS and LHS.
      real*8 uMod(blk%e,nsd), srcRat(blk%e), xmudmi(blk%e,blk%g)
      real*8 IDDESfun(blk%e,1), IDDEStmp(blk%e,1), IDDESfunl(bsz,blk%s,2)
      integer aa, b

c
c... needed for CFL calculation
c
      real*8 rmu_tmp(blk%e), rho_tmp(blk%e), cfll_loc(blk%e)
      rmu_tmp = zero
      rho_tmp = 1.0
      IDDEStmp = zero 
c     
c.... local reconstruction of diffusive flux vector
c
        if ( idiff==2 ) then
           call e3qlSclr (blk, yl, dwl, shp, shgl, xl, ql, sgn)
        endif
c
c.... loop through the integration points
c
        do intqp = 1, blk%g

        if (Qwt(lcsyst,intqp) .eq. zero) cycle          ! precaution
c
c.... get the hierarchic shape functions at this int point
c
        call getshp(blk, intqp,shp,          shgl,      sgn, 
     &              shpfun,        shdrv)
c
c.... get necessary fluid properties
c
        if(iRANS.eq.-5) then ! solving K-W model
          call getdiffsclrKW(blk,intqp,shpfun, shdrv, dwl, yl, xl, diffus, mut, 
     &                     F1, F2)
        else
          call getdiffsclr(blk, shpfun,dwl,yl,diffus,evl)
        endif
c
c.... calculate the integration variables
c
        call e3ivarSclr(blk, intqp, yl,          acl,       shpfun,
     &                  shdrv,       xl,        xmudmi,
     &                  Sclr,        Sdot,      gradS,
     &                  shg,         dxidx,     WdetJ,       
     &                  u1,          u2,        u3,              
     &                  ql,          rLS,       SrcR,
     &                  SrcL,        uMod,      dwl,
     &                  diffus,      srcRat,
     &                  cfllold, gradVl, IDDESfun )
c
c.... compute CFL number
c
        cfll_loc = zero
        call calc_cfl(blk, rho_tmp,      u1,       u2,
     &                u3,        dxidx,    rmu_tmp,
     &                cfll_loc)

        do i=1,blk%s
          cfll(1:blk%e,i) = cfll(1:blk%e,i) + shpfun(:,i)*cfll_loc
        enddo

c
c.... compute the stabilization terms
c
        call e3StabSclr (blk,uMod,    dxidx,   tauS, 
     &                   diffus,  srcR,    giju,
     &                   srcRat)
c
c... computing the DC factor for the discontinuity capturing
c
        if (idcsclr(1) .ne. 0) then
           if ((idcsclr(2).eq.1 .and. isclr.eq.1) .or. 
     &          (idcsclr(2).eq.2 .and. isclr.eq.2) .or.
     &          (idcsclr(2).eq.3 .and. isclr.eq.1) .or.
     &           (idcsclr(2).eq.3 .and. isclr.eq.2)) then ! scalar with dc
c
              call e3dcSclr (blk, gradS,    giju,     gGradS,
     &                        rLS,      tauS,     srcR,
     &                        dcFct)
           endif
        endif                   !end of idcsclr
c
c.... compute the residual contribution at this integration point
c
        call e3ResSclr ( blk, uMod,      gGradS,
     &                   Sclr,      Sdot,       gradS,  
     &                   WdetJ,     rLS,        tauS,
     &                   shpfun,    shg,        srcR,
     &                   diffus, 
     &                   rl )
c
c.... compute the tangent matrix contribution
c
        if (lhs .eq. 1) then
           call e3LHSSclr ( blk, uMod,      giju,       dcFct,
     &                      Sclr,      Sdot,       gradS,  
     &                      WdetJ,     rLS,        tauS,
     &                      shpfun,    shg,        srcL,
     &                      diffus,
     &                      xSebe )

        endif

c
c.... Acculumate IDDES functions for each integration point
        if (ispanIDDES.eq.1) then
          IDDEStmp(:,:) = IDDEStmp(:,:) + IDDESfun(:,:)
        endif 

c
c.... end of integration loop
c
      enddo
c
c.... Divide IDDEStmp by the number of integration points
      IDDEStmp = IDDEStmp/real(blk%g)
      do n=1,blk%s 
         IDDESfunl(1:blk%e,n,1) = one
         IDDESfunl(1:blk%e,n,2) = IDDEStmp(1:blk%e,1)
      enddo

c
c.... return
c
      return
      end
