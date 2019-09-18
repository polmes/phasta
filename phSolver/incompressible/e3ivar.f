      subroutine e3ivar (blk, ith,yl,          acl,       shpfun,
     &                   shgl,        xl,       
     &                   aci,         g1yi,      g2yi,    
     &                   g3yi,        shg,       dxidx,   
     &                   WdetJ,       rho,       pres, 
     &                   u1,          u2,        u3,              
     &                   ql,          rLui,      src,
     &                   rerrl,       rlsl,      rlsli,
     &                   dwl)
c
c----------------------------------------------------------------------
c
c  This routine computes the variables at integration point.
c
c input:
c  yl     (bsz,blk%s,ndof)      : primitive variables
c  acl    (bsz,blk%s,ndof)      : prim.var. accel. 
c  shp    (nen)                 : element shape-functions
c  shgl   (nsd,nen)             : element local-grad-shape-functions
c  xl     (bsz,blk%n,nsd)       : nodal coordinates at current step
c  ql     (bsz,blk%s,nsd*nsd) : diffusive flux vector
c  rlsl   (bsz,blk%s,6)       : resolved Leonard stresses
c
c output:
c  aci    (blk%e,3)              : primvar accel. variables 
c  g1yi   (blk%e,ndof)           : grad-y in direction 1
c  g2yi   (blk%e,ndof)           : grad-y in direction 2
c  g3yi   (blk%e,ndof)           : grad-y in direction 3
c  shg    (blk%e,blk%s,nsd)       : element global grad-shape-functions
c  dxidx  (blk%e,nsd,nsd)        : inverse of deformation gradient
c  WdetJ  (blk%e)                : weighted Jacobian
c  rho    (blk%e)                : density
c  pres   (blk%e)                : pressure
c  u1     (blk%e)                : x1-velocity component
c  u2     (blk%e)                : x2-velocity component
c  u3     (blk%e)                : x3-velocity component
c  rLui   (blk%e,nsd)            : xi-momentum residual
c  src    (blk%e,nsd)            : body force term (not density weighted)
c  rlsli  (blk%e,6)              : resolved Leonard stresses at quad pt
c
c locally calculated and used
c  divqi  (blk%e,nsd+isurf)      : divergence of reconstructed quantity
c
c Zdenek Johan, Summer 1990. (Modified from e2ivar.f)
c Zdenek Johan, Winter 1991. (Fortran 90)
c Kenneth Jansen, Winter 1997. Primitive Variables
c Christian Whiting, Winter 1999. (uBar formulation)
c
c----------------------------------------------------------------------
c
      use  spat_var_eps ! for spatially varying epsilon_ls
      use eblock
      include "common.h"
      type (LocalBlkData) blk

c
c  passed arrays
c
      dimension yl(bsz,blk%s,ndof),        dwl(bsz,blk%n),       
     &            acl(bsz,blk%s,ndof),       shpfun(blk%e,blk%s),
     &            shgl(blk%e,nsd,blk%s),       xl(bsz,blk%n,nsd),
     &            aci(blk%e,nsd),             g1yi(blk%e,ndof),
     &            g2yi(blk%e,ndof),           g3yi(blk%e,ndof),
     &            shg(blk%e,blk%s,nsd),        dxidx(blk%e,nsd,nsd),
     &            WdetJ(blk%e),               
     &            rho(blk%e),                 pres(blk%e),
     &            u1(blk%e),                  u2(blk%e),
     &            u3(blk%e),                  divqi(blk%e,nflow-1+isurf),
     &            ql(bsz,blk%s,idflx),       rLui(blk%e,nsd),
     &            src(blk%e,nsd), Temp(blk%e),xx(blk%e,nsd)
c
        dimension tmp(blk%e), tmp1(blk%e), dkei(blk%e), dist2w(blk%e)
c
        dimension rlsl(bsz,blk%s,6),         rlsli(blk%e,6)
c
        real*8    rerrl(bsz,blk%s,6+isurf), omega(3), divu(blk%e)
        dimension gyti(blk%e,nsd),            gradh(blk%e,nsd),
     &            sforce(blk%e,3),            weber(blk%e),
     &            Sclr(blk%e)
      real*8 epsilon_ls_tmp
c
c.... ------------->  Primitive variables at int. point  <--------------
c
c.... compute primitive variables
c
       rho  = datmat(1,1,1)
       pres = zero
       u1   = zero
       u2   = zero
       u3   = zero
c
       do n = 1, blk%s 
          pres = pres + shpfun(1:blk%e,n) * yl(1:blk%e,n,1)
          u1   = u1   + shpfun(1:blk%e,n) * yl(1:blk%e,n,2)
          u2   = u2   + shpfun(1:blk%e,n) * yl(1:blk%e,n,3)
          u3   = u3   + shpfun(1:blk%e,n) * yl(1:blk%e,n,4)
       enddo
       if(matflg(5,1).eq.2) then ! boussinesq body force
          Temp = zero
          do n = 1, blk%s
             Temp = Temp + shpfun(1:blk%e,n) * yl(1:blk%e,n,5)
          enddo
       endif
       if(matflg(5,1).eq.3.or.matflg(6,1).eq.1) then
c         user-specified body force or coriolis force specified
                 xx = zero
          do n  = 1,blk%n
             xx(1:blk%e,1) = xx(1:blk%e,1)  + shpfun(1:blk%e,n) * xl(1:blk%e,n,1)
             xx(1:blk%e,2) = xx(1:blk%e,2)  + shpfun(1:blk%e,n) * xl(1:blk%e,n,2)
             xx(1:blk%e,3) = xx(1:blk%e,3)  + shpfun(1:blk%e,n) * xl(1:blk%e,n,3)
          enddo
       endif
c
       if(iRANS.eq.-2.or.iRANS.eq.-5) then ! kay-epsilon and SST
          dist2w = zero
          do n = 1, blk%n
             dist2w = dist2w + shpfun(1:blk%e,n) * dwl(1:blk%e,n)
          enddo
       endif
c
 
       if( (iLES.gt.10).and.(iLES.lt.20))  then  ! bardina
       rlsli = zero
       do n = 1, blk%s 

          rlsli(1:blk%e,1) = rlsli(1:blk%e,1) + shpfun(1:blk%e,n) * rlsl(1:blk%e,n,1)
          rlsli(1:blk%e,2) = rlsli(1:blk%e,2) + shpfun(1:blk%e,n) * rlsl(1:blk%e,n,2)
          rlsli(1:blk%e,3) = rlsli(1:blk%e,3) + shpfun(1:blk%e,n) * rlsl(1:blk%e,n,3)
          rlsli(1:blk%e,4) = rlsli(1:blk%e,4) + shpfun(1:blk%e,n) * rlsl(1:blk%e,n,4)
          rlsli(1:blk%e,5) = rlsli(1:blk%e,5) + shpfun(1:blk%e,n) * rlsl(1:blk%e,n,5)
          rlsli(1:blk%e,6) = rlsli(1:blk%e,6) + shpfun(1:blk%e,n) * rlsl(1:blk%e,n,6)

       enddo
       else
          rlsli = zero
       endif
c
c.... ----------------------->  accel. at int. point  <----------------------
c
       aci = zero
       do n = 1, blk%s
          aci(1:blk%e,1) = aci(1:blk%e,1) + shpfun(1:blk%e,n) * acl(1:blk%e,n,2)
          aci(1:blk%e,2) = aci(1:blk%e,2) + shpfun(1:blk%e,n) * acl(1:blk%e,n,3)
          aci(1:blk%e,3) = aci(1:blk%e,3) + shpfun(1:blk%e,n) * acl(1:blk%e,n,4)
       enddo
c
c.... --------------------->  Element Metrics  <-----------------------
c
       call e3metric(blk,ith, xl,         shgl,       dxidx,  
     &                shg,        WdetJ)
c
c.... compute the global gradient of u and P
c
c
       g1yi = zero
       g2yi = zero
       g3yi = zero
       do n = 1, blk%s
          g1yi(1:blk%e,1) = g1yi(1:blk%e,1) + shg(1:blk%e,n,1) * yl(1:blk%e,n,1)
          g1yi(1:blk%e,2) = g1yi(1:blk%e,2) + shg(1:blk%e,n,1) * yl(1:blk%e,n,2)
          g1yi(1:blk%e,3) = g1yi(1:blk%e,3) + shg(1:blk%e,n,1) * yl(1:blk%e,n,3)
          g1yi(1:blk%e,4) = g1yi(1:blk%e,4) + shg(1:blk%e,n,1) * yl(1:blk%e,n,4)
c
          g2yi(1:blk%e,1) = g2yi(1:blk%e,1) + shg(1:blk%e,n,2) * yl(1:blk%e,n,1)
          g2yi(1:blk%e,2) = g2yi(1:blk%e,2) + shg(1:blk%e,n,2) * yl(1:blk%e,n,2)
          g2yi(1:blk%e,3) = g2yi(1:blk%e,3) + shg(1:blk%e,n,2) * yl(1:blk%e,n,3)
          g2yi(1:blk%e,4) = g2yi(1:blk%e,4) + shg(1:blk%e,n,2) * yl(1:blk%e,n,4)
c
          g3yi(1:blk%e,1) = g3yi(1:blk%e,1) + shg(1:blk%e,n,3) * yl(1:blk%e,n,1)
          g3yi(1:blk%e,2) = g3yi(1:blk%e,2) + shg(1:blk%e,n,3) * yl(1:blk%e,n,2)
          g3yi(1:blk%e,3) = g3yi(1:blk%e,3) + shg(1:blk%e,n,3) * yl(1:blk%e,n,3)
          g3yi(1:blk%e,4) = g3yi(1:blk%e,4) + shg(1:blk%e,n,3) * yl(1:blk%e,n,4)
       enddo

       divqi = zero
       idflow = 3
       if ( idiff >= 1 .or. isurf==1 ) then
c     
c.... compute divergence of diffusive flux vector, qi,i
c     
          if(idiff >= 1) then
             do n=1, blk%s
                divqi(1:blk%e,1) = divqi(1:blk%e,1) + shg(1:blk%e,n,1)*ql(1:blk%e,n,1 ) 
     &                                  + shg(1:blk%e,n,2)*ql(1:blk%e,n,4 )
     &                                  + shg(1:blk%e,n,3)*ql(1:blk%e,n,7 )

                divqi(1:blk%e,2) = divqi(1:blk%e,2) + shg(1:blk%e,n,1)*ql(1:blk%e,n,2 ) 
     &                                  + shg(1:blk%e,n,2)*ql(1:blk%e,n,5 )
     &                                  + shg(1:blk%e,n,3)*ql(1:blk%e,n,8)

                divqi(1:blk%e,3) = divqi(1:blk%e,3) + shg(1:blk%e,n,1)*ql(1:blk%e,n,3 ) 
     &                                  + shg(1:blk%e,n,2)*ql(1:blk%e,n,6 )
     &                                  + shg(1:blk%e,n,3)*ql(1:blk%e,n,9 )

          enddo

          endif                 !end of idiff
c     
          if (isurf .eq. 1) then   
c     .... divergence of normal calculation (curvature)
c  note ql(:,:,idflx-2 through idflx) contains the normal vector
c  e3q computed grad(phi) and qpbc divided by the magnitude of grad(phi).
             do n=1, blk%s
                divqi(1:blk%e,idflow+1) = divqi(1:blk%e,idflow+1) 
     &               + shg(1:blk%e,n,1)*ql(1:blk%e,n,idflx-2)
     &               + shg(1:blk%e,n,2)*ql(1:blk%e,n,idflx-1)
     &               + shg(1:blk%e,n,3)*ql(1:blk%e,n,idflx)
             enddo 
c     .... initialization of some variables
             Sclr = zero
             gradh= zero
             gyti = zero
             sforce=zero
             do i = 1, blk%e
                do n = 1, blk%s      
                   Sclr(i) = Sclr(i) + shpfun(i,n) * yl(i,n,6) !scalar
c     
c     .... compute the global gradient of Scalar variable
c     
                   gyti(i,1) = gyti(i,1) + shg(i,n,1) * yl(i,n,6) 
                   gyti(i,2) = gyti(i,2) + shg(i,n,2) * yl(i,n,6)
                   gyti(i,3) = gyti(i,3) + shg(i,n,3) * yl(i,n,6)
c     
                enddo
c
c..Now we have to use spatially varied epsilon_ls
c
                epsilon_ls_tmp = epsilon_ls * 
     &               elem_local_size(blk%i +i-1)
  
                if (abs (sclr(i)) .le. epsilon_ls_tmp) then
                   gradh(i,1) = 0.5/epsilon_ls_tmp * (1.0 
     &                  + cos(pi*Sclr(i)/epsilon_ls_tmp)) * gyti(i,1)
                   gradh(i,2) = 0.5/epsilon_ls_tmp * (1.0 
     &                  + cos(pi*Sclr(i)/epsilon_ls_tmp)) * gyti(i,2) 
                   gradh(i,3) = 0.5/epsilon_ls_tmp * (1.0 
     &                  + cos(pi*Sclr(i)/epsilon_ls_tmp)) * gyti(i,3)
                endif
             enddo              !end of the loop over blk%e
c     
c .. surface tension force calculation
c .. divide by density now as it gets multiplied in e3res.f, as surface
c    tension force is already in the form of force per unti volume
c     
             weber(:) = Bo
             sforce(:,1) = -(1.0/weber(:)) * divqi(:,idflow+1) !x-direction
     &            *gradh(:,1) /rho(:)
             sforce(:,2) = -(1.0/weber(:)) * divqi(:,idflow+1) !y-direction
     &            *gradh(:,2) /rho(:)
             sforce(:,3) = -(1.0/weber(:)) * divqi(:,idflow+1) !z-direction
     &            *gradh(:,3) /rho(:)          
c
          endif        ! end of the surface tension force calculation
       endif           ! diffusive flux computation
c
c Calculate strong form of pde as well as the source term
c      
       call e3resStrongPDE(blk,
     &      aci,  u1,   u2,   u3,   Temp, rho,  xx,
     &            g1yi, g2yi, g3yi,
     &      rLui, src, divqi, sforce)      
c
c.... -------------------> error calculation  <-----------------
c     
        if((ierrcalc.eq.1).and.(nitr.eq.iter)) then
          do ia=1,blk%s
             tmp=shpfun(:,ia)*WdetJ(:)
             tmp1=shpfun(:,ia)*Qwt(lcsyst,ith) 
             rerrl(1:blk%e,ia,1) = rerrl(1:blk%e,ia,1) +
     &                       tmp1(1:blk%e)*rLui(1:blk%e,1)*rLui(1:blk%e,1)
             rerrl(1:blk%e,ia,2) = rerrl(1:blk%e,ia,2) +
     &                       tmp1(1:blk%e)*rLui(1:blk%e,2)*rLui(1:blk%e,2)
             rerrl(1:blk%e,ia,3) = rerrl(1:blk%e,ia,3) +
     &                       tmp1(1:blk%e)*rLui(1:blk%e,3)*rLui(1:blk%e,3)

             rerrl(1:blk%e,ia,4) = rerrl(1:blk%e,ia,4) +
     &                       tmp(1:blk%e)*divqi(1:blk%e,1)*divqi(1:blk%e,1)
             rerrl(1:blk%e,ia,5) = rerrl(1:blk%e,ia,5) +
     &                       tmp(1:blk%e)*divqi(1:blk%e,2)*divqi(1:blk%e,2)
             rerrl(1:blk%e,ia,6) = rerrl(1:blk%e,ia,6) +
     &                       tmp(1:blk%e)*divqi(1:blk%e,3)*divqi(1:blk%e,3)
          if (isurf .eq. 1) rerrl(1:blk%e,ia,7) = rerrl(1:blk%e,ia,7) +
     &                       tmp(:)*divqi(:,idflow+1)*divqi(:,idflow+1)
          enddo
       endif
       distcalc=0  ! return to 1 if you want to compute T-S instability
       if(distcalc.eq.1) then
c
c.... ----------------------->  dist. kin energy at int. point  <--------------
c
       
       if (ires .ne. 2 .and. iter.eq.1)  then  !only do at beginning of step
c
c calc exact velocity for a channel at quadrature points.
c
       dkei=0.0
c
       do n = 1, blk%n 
          dkei = dkei + shpfun(1:blk%e,n) * (1.0-xl(1:blk%e,n,2)**2) !u_ex^~ (in FEM space)
       enddo
          dkei = (u1(:)-dkei)**2 +u2(:)**2  ! u'^2+v'^2
          dkei = dkei*WdetJ  ! mult function*W*det of jacobian to
c                              get this quadrature point contribution
          dke  = dke+sum(dkei) ! we move the sum over elements inside of the
c                              sum over quadrature to save memory (we want
c                              a scalar only)
       endif
       endif
c     
c.... return
c
       return
       end

c-----------------------------------------------------------------------
c 
c     Calculate the variables for the scalar advection-diffusion
c     equation.
c
c-----------------------------------------------------------------------
      subroutine e3ivarSclr (blk, ith, yl,          acl,       shpfun,
     &                      shgl,        xl,        xmudmi,
     &                      Sclr,        Sdot,      gradS,  
     &                      shg,         dxidx,     WdetJ,
     &                      u1,          u2,        u3,              
     &                      ql,          rLS ,       SrcR,
     &                      SrcL,        uMod,      dwl,
     &                      diffus,      srcRat,
     &                      cfll, gradVl, IDDESfun )
c
      use spat_var_eps   ! use spatially-varying epl_ls
      use eblock
      use spanStats
      include "common.h"
      type (LocalBlkData) blk

c
c  passed arrays
c
      dimension yl(bsz,blk%s,ndof),        acl(bsz,blk%s,ndof), 
     &          Sclr(blk%e),                Sdot(blk%e),
     &          gradS(blk%e,nsd),           shpfun(blk%e,blk%s),
     &          shgl(blk%e,nsd,blk%s),       xl(bsz,blk%n,nsd),
     &          shg(blk%e,blk%s,nsd),        dxidx(blk%e,nsd,nsd),
     &          WdetJ(blk%e),              
     &          u1(blk%e),                  u2(blk%e),
     &          u3(blk%e),                  divS(blk%e),
     &          ql(bsz,blk%s,nsd),         rLS(blk%e),
     &          SrcR(blk%e),                 SrcL(blk%e),
     &          dwl(bsz,blk%n),            diffus(blk%e),
     &          umod(blk%e,nsd), Temp(blk%e),xx(blk%e,nsd),
     &          divqi(blk%e) ,
     &          cfll(bsz,blk%s),   gradVl(bsz,blk%s,nsdsq) 
c
      dimension tmp(blk%e), srcRat(blk%e)
      real*8 rLui(blk%e,nsd),     aci(blk%e,nsd),
     &       g1yi(blk%e,nflow),   g2yi(blk%e,nflow),
     &       g3yi(blk%e,nflow),
     &       src(blk%e,nsd),      rho(blk%e),
     &       rmu(blk%e),          ssq(blk%e)
      real*8 uBar(blk%e,nsd), xmudmi(blk%e,ngauss)
      real*8 denominv(blk%e)
      real*8 IDDESfun(blk%e,nfun)
      integer comp

c
c.... ------------->  Primitive variables at int. point  <--------------
c
c.... compute primitive variables
c
      u1   = zero
      u2   = zero
      u3   = zero
      Sclr = zero
c
      id=isclr+5
      do n = 1, blk%s 
         u1   = u1   + shpfun(1:blk%e,n) * yl(1:blk%e,n,2)
         u2   = u2   + shpfun(1:blk%e,n) * yl(1:blk%e,n,3)
         u3   = u3   + shpfun(1:blk%e,n) * yl(1:blk%e,n,4)
         Sclr = Sclr + shpfun(1:blk%e,n) * yl(1:blk%e,n,id)
      enddo
c
c
c.... ----------------------->  dS/dt at int. point  <----------------------
c
      Sdot = zero
      do n = 1, blk%s
         Sdot = Sdot + shpfun(1:blk%e,n) * acl(1:blk%e,n,id)
      enddo
c
c.... --------------------->  Element Metrics  <-----------------------
c

      call e3metric(blk,ith, xl,         shgl,        dxidx,  
     &               shg,        WdetJ)

c
c.... compute the global gradient of u and P
c
c
       gradS = zero
       do n = 1, blk%s
          gradS(1:blk%e,1) = gradS(1:blk%e,1) + shg(1:blk%e,n,1) * yl(1:blk%e,n,id)
          gradS(1:blk%e,2) = gradS(1:blk%e,2) + shg(1:blk%e,n,2) * yl(1:blk%e,n,id)
          gradS(1:blk%e,3) = gradS(1:blk%e,3) + shg(1:blk%e,n,3) * yl(1:blk%e,n,id)
       enddo

       divS = zero
       if ( idiff >= 1 ) then
c
c.... compute divergence of diffusive flux vector, qi,i
c
          do n=1, blk%s
             divS(1:blk%e) = divS(1:blk%e) + shg(1:blk%e,n,1)*ql(1:blk%e,n,1 ) 
     &                         + shg(1:blk%e,n,2)*ql(1:blk%e,n,2 ) 
     &                         + shg(1:blk%e,n,3)*ql(1:blk%e,n,3 ) 
          enddo
       endif                    ! diffusive flux computation

       if(consrv_sclr_conv_vel.eq.1) then
c         Calculate uBar = u - TauM*L, where TauM is the momentum
c         stabilization factor and L is the momentum residual

          if(matflg(5,1).eq.2) then ! boussinesq body force
             Temp = zero
             do n = 1, blk%s
                Temp = Temp + shpfun(1:blk%e,n) * yl(1:blk%e,n,5)
             enddo
          endif
          if(matflg(5,1).eq.3.or.matflg(6,1).eq.1) then
c     user-specified body force or coriolis force specified
             xx = zero
             do n  = 1,blk%n
                xx(1:blk%e,1) = xx(1:blk%e,1)  + shpfun(1:blk%e,n) * xl(1:blk%e,n,1)
                xx(1:blk%e,2) = xx(1:blk%e,2)  + shpfun(1:blk%e,n) * xl(1:blk%e,n,2)
                xx(1:blk%e,3) = xx(1:blk%e,3)  + shpfun(1:blk%e,n) * xl(1:blk%e,n,3)
             enddo
          endif
          aci = zero
          do n = 1, blk%s
             aci(1:blk%e,1) = aci(1:blk%e,1) + shpfun(1:blk%e,n) * acl(1:blk%e,n,2)
             aci(1:blk%e,2) = aci(1:blk%e,2) + shpfun(1:blk%e,n) * acl(1:blk%e,n,3)
             aci(1:blk%e,3) = aci(1:blk%e,3) + shpfun(1:blk%e,n) * acl(1:blk%e,n,4)
          enddo
          g1yi = zero
          g2yi = zero
          g3yi = zero
          do n = 1, blk%s
             g1yi(1:blk%e,1) = g1yi(1:blk%e,1) + shg(1:blk%e,n,1) * yl(1:blk%e,n,1)
             g1yi(1:blk%e,2) = g1yi(1:blk%e,2) + shg(1:blk%e,n,1) * yl(1:blk%e,n,2)
             g1yi(1:blk%e,3) = g1yi(1:blk%e,3) + shg(1:blk%e,n,1) * yl(1:blk%e,n,3)
             g1yi(1:blk%e,4) = g1yi(1:blk%e,4) + shg(1:blk%e,n,1) * yl(1:blk%e,n,4)
c     
             g2yi(1:blk%e,1) = g2yi(1:blk%e,1) + shg(1:blk%e,n,2) * yl(1:blk%e,n,1)
             g2yi(1:blk%e,2) = g2yi(1:blk%e,2) + shg(1:blk%e,n,2) * yl(1:blk%e,n,2)
             g2yi(1:blk%e,3) = g2yi(1:blk%e,3) + shg(1:blk%e,n,2) * yl(1:blk%e,n,3)
             g2yi(1:blk%e,4) = g2yi(1:blk%e,4) + shg(1:blk%e,n,2) * yl(1:blk%e,n,4)
c     
             g3yi(1:blk%e,1) = g3yi(1:blk%e,1) + shg(1:blk%e,n,3) * yl(1:blk%e,n,1)
             g3yi(1:blk%e,2) = g3yi(1:blk%e,2) + shg(1:blk%e,n,3) * yl(1:blk%e,n,2)
             g3yi(1:blk%e,3) = g3yi(1:blk%e,3) + shg(1:blk%e,n,3) * yl(1:blk%e,n,3)
             g3yi(1:blk%e,4) = g3yi(1:blk%e,4) + shg(1:blk%e,n,3) * yl(1:blk%e,n,4)
          enddo
c          
          if (iLSet .eq. 0)then
             rho  = datmat(1,1,1)
             rmu = datmat(1,2,1)
          else
             write(*,*) 'Not sure if we can handle level set with K-E'
             write(*,*) '(different uMods? correct value of rho?)'
          endif
          sforce=zero ! until we bring sforce computation to scalar solve
          divqi=zero  ! until we reconstruct q_flow for scalar solve
          call e3resStrongPDE( blk,
     &         aci,  u1,   u2,   u3,   Temp, rho,  x,
     &               g1yi, g2yi, g3yi,
     &         rLui, src, divqi,sforce)
c     this call to getdiff is for the flow diffusive properties

      if(iRANS.eq.-5) then ! solving K-W model
            ssq(:)    = sqrt( 2.0d0 * (g1yi(:,2) ** 2
     &         + g2yi(:,3) ** 2
     &         + g3yi(:,4) ** 2
     &         + 0.5d0
     &         * ( (g3yi(:,3) + g2yi(:,4)) ** 2
     &           + (g1yi(:,4) + g3yi(:,2)) ** 2
     &           + (g2yi(:,2) + g1yi(:,3)) ** 2 )))
            call AddEddyViscKomgq(yl, dwl, shgl, shpfun, rho, rmu, ssq, xl)
       else
         call getdiff(blk,ith,dwl, yl, shpfun, xmudmi, xl, rmu, rho,
     &               elem_local_size(blk%i),
     &               evl)
       endif
! bad to have a separate routine for something like Tau that may vary          call e3uBar(blk,rho, src, dxidx, rLui, rmu, uBar)
! changing to call e3stab that now has conditionals to skip computation of tauc and taubar so should be o.k. that they are note
! dimensioned in this case since memory not accessed
          call e3stab (blk, rho,          u1,       u2,
     &                 u3,           dxidx,    rLui,
     &                 rmu,          tauC,     tauM,
     &                 tauBar,       uBar )

          u1=ubar(:,1)          ! the entire scalar residual
          u2=ubar(:,2)          ! is based on the modified
          u3=ubar(:,3)          ! velocity for conservation
       endif
c
c.... Initialize uMod, the modified velocity uMod
c      We initialize it to u_i and then calculate
c      the correction in e3sourcesclr
c

       umod(:,1) = u1
       umod(:,2) = u2
       umod(:,3) = u3
c     
c.... compute  source terms
c
!
!    if we are solving the redistancing equation, the umod(:,:) are 
!    modified in e3sourceSclr.  
!
!  if we are redistancing levelset variable we want to use a use the  
!  convective term from the equation.  


       if(nosource.ne.1) then
        srcRat = zero
        call e3sourceSclr ( blk, Sclr,         Sdot,      gradS,  dwl,
     &                      shpfun,       shg,       yl,     dxidx,
     &                      diffus,       u1,        u2,     u3,
     &                      xl,           srcR,      srcL,   uMod,
     &                      srcRat,  cfll , gradVl, IDDESfun )
       else
        srcRat = zero
        srcR   = zero
        srcL   = zero
       endif
c
c.... -------------------> Scalar residual  <-----------------
c

         rLS(:) = ( Sdot(:) +  (u1*gradS(:,1) + 
     &                              u2*gradS(:,2) +
     &                              u3*gradS(:,3)) )
     &        - divS(:)           

c
c.... return
c
       return
       end

