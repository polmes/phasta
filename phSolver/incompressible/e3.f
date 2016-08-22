        subroutine e3 (yl,      acl,     dwl,     shp,
     &                 shgl,    xl,      rl,      ql,
     &                 xKebe,   xGoC,    xmudmi,  sgn, 
     &                 rerrl, rlsl)
c                                                                      
c----------------------------------------------------------------------
c
c     This routine calculates the residual and tangent matrix for the 
c     UBar formulation of the incompressible Navier Stokes equations.
c
c
c input:    e    a   1..5   when we think of U^e_a  and U is 5 variables
c  yl     (npro,nshl,ndof)      : Y variables (not U)
c  acl    (npro,nshl,ndof)      : Y acceleration (Y_{,t})
c  shp    (nen,ngauss)           : element shape-functions  N_a
c  shgl   (nsd,nen,ngauss)       : element local-shape-functions N_{a,xi}
c  wght   (ngauss)               : element weight (for quadrature)
c  xl     (npro,nenl,nsd)       : nodal coordinates at current step (x^e_a)
c  ql     (npro,nshl,nsd*nsd) : diffusive flux vector (don't worry)
c  rlsl   (npro,nshl,6)       : resolved Leonard stresses
c
c output:
c  rl     (npro,nshl,nflow)      : element RHS residual    (G^e_a)
c  rml    (npro,nshl,nflow)      : element modified residual  (G^e_a tilde)
c  xKebe  (npro,9,nshl,nshl)  : element LHS tangent mass matrix
c  xGoC   (npro,4,nshl,nshl)    : element LHS tangent mass matrix
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
      use omp_lib
      include "common.h"

c
        dimension yl(npro,nshl,ndof),
     &            acl(npro,nshl,ndof),       
     &            shp(nshl,ngauss),       shgl(nsd,nshl,ngauss),
     &            xl(npro,nenl,nsd),      dwl(npro,nenl),
     &            rl(npro,nshl,nflow),     ql(npro,nshl,idflx)
c      
        real*8, allocatable, dimension(:,:,:,:,:) :: xK_qp, xG_qp
        real*8, allocatable, dimension(:,:,:,:) :: rl_qp

        dimension xKebe(npro,9,nshl,nshl), xGoC(npro,4,nshl,nshl)
c
c.... local declarations
c
        dimension g1yi(npro,ndof),        g2yi(npro,ndof),
     &            g3yi(npro,ndof),        shg(npro,nshl,nsd),
     &            aci(npro,3),            dxidx(npro,nsd,nsd),       
     &            WdetJ(npro),            rho(npro),
     &            pres(npro),             u1(npro),
     &            u2(npro),               u3(npro),
     &            rLui(npro,nsd),         uBar(npro,nsd),
     &            xmudmi(npro,ngauss),     sgn(npro,nshl), 
     &            shpfun(npro,nshl),      shdrv(npro,nsd,nshl),
     &            rmu(npro),              tauC(npro),
     &            tauM(npro),             tauBar(npro),
     &            src(npro,3)

        dimension rlsl(npro,nshl,6),      rlsli(npro,6)

        real*8    rerrl(npro,nshl,6)
        integer   aa

c
c     
c.... local reconstruction of diffusive flux vector for quadratics
c     or greater but NOT for bflux since local mass was not mapped
c
        if ( idiff==2 .and. ires .eq. 1 ) then
           call e3ql (yl,        dwl,       shp,       shgl, 
     &                xl,        ql,        xmudmi, 
     &                sgn)
        endif
c
c.... loop through the integration points
c
! for now if we thread it is ngauss        nthread_qp=1
	allocate( rl_qp(npro,nshl,nflow,ngauss))
        allocate( xK_qp(npro,9,nshl,nshl,ngauss))
        allocate( xG_qp(npro,4,nshl,nshl,ngauss))
        rl_qp=zero
        xK_qp=zero
        xG_qp=zero
!$OMP  parallel do 
!$OMP& private (ith,sgn,shpfun,shdrv,rmu,rho,aci,g1yi,g2yi,g3yi)
!$OMP& private (shg,dxidx,WdetJ,pres,u1,u2,u3,rLui,src,rlsi)
!$OMP& private (tauC,tauM,tauBar,uBar,id)
!$OMP& shared (shp,shgl,dwl,yl,xmudmi,xl,acl,rlsl)
!$OMP& shared (rl_qp,xK_qp,xG_qp)
        do ith = 1, ngauss
        if(npro.eq.-1) then
          id = omp_get_thread_num ( )
          write (*,*) 'HELLO from rank, process,ith ', myrank,id,ith
        endif

        if (Qwt(lcsyst,ith) .eq. zero) cycle          ! precaution
c
c.... get the hierarchic shape functions at this int point
c
        call getshp(ith, shp,          shgl,      sgn, 
     &              shpfun,       shdrv)
c
c.... get necessary fluid properties (including eddy viscosity)
c
        call getdiff(ith, dwl,  yl,     shpfun,     xmudmi, xl,   rmu, rho)
c
c.... calculate the integration variables
c
        call e3ivar (ith,yl,          acl,       shpfun,
     &               shdrv,       xl,
     &               aci,         g1yi,      g2yi,    
     &               g3yi,        shg,       dxidx,   
     &               WdetJ,       rho,       pres, 
     &               u1,          u2,        u3,              
     &               ql,          rLui,      src,
     &               rerrl,       rlsl,      rlsli,
     &               dwl) 
c
c.... compute the stabilization terms
c
        call e3stab (rho,          u1,       u2,
     &               u3,           dxidx,    rLui,   
     &               rmu,          tauC,     tauM,   
     &               tauBar,       uBar )  
c
c.... compute the residual contribution at this integration point
c
        call e3Res ( u1,        u2,         u3,
     &               uBar,      aci,        WdetJ,
     &               g1yi,      g2yi,       g3yi,
     &               rLui,      rmu,        rho,
     &               tauC,      tauM,       tauBar,
     &               shpfun,    shg,        src,
     &               rl_qp(:,:,:,ith),        pres,       acl,
     &               rlsli)
c
c.... compute the tangent matrix contribution
c
        if (lhs .eq. 1) then
           call e3LHS ( u1,        u2,         u3,
     &                  uBar,      WdetJ,      rho,
     &                  rLui,      rmu,
     &                  tauC,      tauM,       tauBar,
     &                  shpfun,    shg,        xK_qp(:,:,:,:,ith),
     &                  xG_qp(:,:,:,:,ith))
        endif

c
c.... end of integration loop
c
      enddo
c
c here we accumulate the thread work
c
      do i=1,ngauss
       rl(1:npro,1:nshl,1:nflow)=rl(1:npro,1:nshl,1:nflow)
     &                       +rl_qp(1:npro,1:nshl,1:nflow,i)
       xKebe(1:npro,1:9,1:nshl,1:nshl)=xKebe(1:npro,1:9,1:nshl,1:nshl)
     &                                +xK_qp(1:npro,1:9,1:nshl,1:nshl,i)
       xGoC(1:npro,1:4,1:nshl,1:nshl)=xGoC(1:npro,1:4,1:nshl,1:nshl)
     &                              +xG_qp(1:npro,1:4,1:nshl,1:nshl,i)
      enddo

      deallocate(xG_qp)
      deallocate(xK_qp)
      deallocate(rl_qp)
c
c.... symmetrize C
c
      if (lhs .eq. 1) then
         do ib = 1, nshl
            do iaa = 1, ib-1
               xGoC(1:npro,4,iaa,ib) = xGoC(1:npro,4,ib,iaa)
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


      subroutine e3Sclr (yl,      acl,     shp,
     &                     shgl,    xl,      dwl,
     &                     rl,      ql,      xSebe,   
     &                     sgn,     xmudmi)
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
      include "common.h"
c
      real*8    yl(npro,nshl,ndof),     acl(npro,nshl,ndof),       
     &            shp(nshl,ngauss),       shgl(nsd,nshl,ngauss),
     &            xl(npro,nenl,nsd),      rl(npro,nshl),          
     &            ql(npro,nshl,nsd),      xSebe(npro,nshl,nshl),
     &            dwl(npro,nshl)
c
c.... local declarations
c
      real*8    gradS(npro,nsd),        shg(npro,nshl,nsd),
     &            Sdot(npro),             Sclr(npro),
     &            dxidx(npro,nsd,nsd),    WdetJ(npro),      
     &            u1(npro),     u2(npro), u3(npro),
     &            sgn(npro,nshl),         shpfun(npro,nshl),       
     &            shdrv(npro,nsd,nshl),   rLS(npro),
     &            tauS(npro),             diffus(npro),
     &            srcL(npro),             srcR(npro),
     &            gGradS(npro,nsd),       dcFct(npro),
     &            giju(npro,6)
c
c.... Source terms sometimes take the form (beta_i)*(phi,_i).  Since
c     the convective term has (u_i)*(phi,_i), it is useful to treat
c     beta_i as a "correction" to the velocity.  In calculating the
c     stabilization terms, the new "modified" velocity (u_i-beta_i) is 
c     then used in place of the pure velocity for stabilization terms,
c     and the source term sneaks into the RHS and LHS.
      real*8    uMod(npro,nsd), srcRat(npro), xmudmi(npro,ngauss)
c
      integer   aa, b
c     
c.... local reconstruction of diffusive flux vector
c
        if ( idiff==2 ) then
           call e3qlSclr (yl, dwl, shp, shgl, xl, ql, sgn)
        endif
c
c.... loop through the integration points
c
        do intp = 1, ngauss

        if (Qwt(lcsyst,intp) .eq. zero) cycle          ! precaution
c
c.... get the hierarchic shape functions at this int point
c
        call getshp(intp,shp,          shgl,      sgn, 
     &              shpfun,        shdrv)
c
c.... get necessary fluid properties
c
        call getdiffsclr(shpfun,dwl,yl,diffus)
c
c.... calculate the integration variables
c
        call e3ivarSclr(yl,          acl,       shpfun,
     &                  shdrv,       xl,        xmudmi,
     &                  Sclr,        Sdot,      gradS,
     &                  shg,         dxidx,     WdetJ,       
     &                  u1,          u2,        u3,              
     &                  ql,          rLS,       SrcR,
     &                  SrcL,        uMod,      dwl,
     &                  diffus,      srcRat)


c
c.... compute the stabilization terms
c
        call e3StabSclr (uMod,    dxidx,   tauS, 
     &                   diffus,  srcR,    giju,
     &                   srcRat)
c
c... computing the DC factor for the discontinuity capturing
c
        if (idcsclr(1) .ne. 0) then
           if ((idcsclr(2).eq.1 .and. isclr.eq.1) .or. 
     &          (idcsclr(2).eq.2 .and. isclr.eq.2)) then ! scalar with dc
c
              call e3dcSclr ( gradS,    giju,     gGradS,
     &                        rLS,      tauS,     srcR,
     &                        dcFct)
           endif
        endif                   !end of idcsclr
c
c.... compute the residual contribution at this integration point
c
        call e3ResSclr ( uMod,      gGradS,
     &                   Sclr,      Sdot,       gradS,  
     &                   WdetJ,     rLS,        tauS,
     &                   shpfun,    shg,        srcR,
     &                   diffus, 
     &                   rl )
c
c.... compute the tangent matrix contribution
c
        if (lhs .eq. 1) then
           call e3LHSSclr ( uMod,      giju,       dcFct,
     &                      Sclr,      Sdot,       gradS,  
     &                      WdetJ,     rLS,        tauS,
     &                      shpfun,    shg,        srcL,
     &                      diffus,
     &                      xSebe )

        endif

c
c.... end of integration loop
c
      enddo

c
c.... return
c
      return
      end
