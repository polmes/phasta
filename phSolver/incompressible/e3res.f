      subroutine e3Res ( blk, u1,        u2,         u3,
     &                   uBar,      aci,        WdetJ,
     &                   g1yi,      g2yi,       g3yi,
     &                   rLui,      rmu,        rho,
     &                   tauC,      tauM,       tauBar,
     &                   shpfun,    shg,        src,
     &                   rl,        pres,       
     &                   rlsli)
c------------------------------------------------------------------------
c 
c  This routine computes the residual vector at an
c  integration point.
c
c  input:
c     u1(blk%e)                  : x1-velocity
c     u2(blk%e)                  : x2-velocity
c     u3(blk%e)                  : x3-velocity
c     uBar(blk%e,3)              : u - tauM * Li
c     aci(blk%e,3)               : acceleration
c     rlsli(blk%e,6)             : resolved Leonard stresses
c     WdetJ(blk%e)               : weighted jacobian determinant
c     g1yi(blk%e,ndof)              : x1-gradient of variables
c     g2yi(blk%e,ndof)              : x2-gradient of variables
c     g3yi(blk%e,ndof)              : x3-gradient of variables
c     rLui(blk%e,3)              : total residual of NS equations
c     rmu(blk%e)                 : fluid viscosity
c     rho(blk%e)                 : density
c     tauC(blk%e)                : continuity tau
c     tauM(blk%e)                : momentum tau
c     tauBar(blk%e)              : additional tau
c     shpfun(blk%e,blk%s)         : element shape functions
c     shg(blk%e,blk%s,nsd)        : global grad of element shape functions
c     src(blk%e,nsd)             : body force term
c
c  output:
c     rl(bsz,blk%s,nflow)
c
c------------------------------------------------------------------------
      use eblock
      include "common.h"
      type (LocalBlkData) blk

      dimension u1(blk%e),         u2(blk%e),         u3(blk%e),
     &          uBar(blk%e,nsd),   aci(blk%e,nsd),    WdetJ(blk%e),
     &          g1yi(blk%e,nflow), g2yi(blk%e,nflow), g3yi(blk%e,nflow),
     &          rLui(blk%e,nsd),   rmu(blk%e),        rho(blk%e),
     &          tauC(blk%e),       tauM(blk%e),       tauBar(blk%e),
     &          shpfun(blk%e,blk%s),shg(blk%e,blk%s,nsd), src(blk%e,nsd),
     &          pres(blk%e)
      
      dimension rl(bsz,blk%s,nflow),
     &          rlsli(blk%e,6)
c
c.... local declarations
c
      real*8    tmp1(blk%e),   tmp2(blk%e),          tmp3(blk%e), 
     &          tmp(blk%e),    rGNa(blk%e,nsd,nsd),  rNa(blk%e,nsd),
     &          locmass(blk%e,blk%s),omega(3)

      integer aa
c     
c.... initialize multipliers for Na and Na_{,i}
c
      rNa  = zero
      rGNa = zero
c
c.... compute the Na multiplier
c
      tmps     = one-flmpr  ! consistant mass factor

c
c no density yet...it comes later
c
      rNa(:,1) = aci(:,1)  * tmps
     &         - src(:,1)
      rNa(:,2) = aci(:,2)  * tmps
     &         - src(:,2)
      rNa(:,3) = aci(:,3)  * tmps
     &         - src(:,3)

c
c.... rotating frame terms if needed
c

      if(matflg(6,1).eq.1) then ! rotation

         omega(1)=datmat(1,6,1)
         omega(2)=datmat(2,6,1)
         omega(3)=datmat(3,6,1)
c
c no density yet...it comes later
c
         
         rNa(:,1) = rNa(:,1) + (omega(2)-omega(3)) * tauM * rLui(:,1)
         rNa(:,2) = rNa(:,2) + (omega(3)-omega(1)) * tauM * rLui(:,2)
         rNa(:,3) = rNa(:,3) + (omega(1)-omega(2)) * tauM * rLui(:,3)
      endif


c
c.... compute the Na,i multiplier
c
      tmp  = -pres + tauC * (g1yi(:,2) + g2yi(:,3) + g3yi(:,4))
      tmp1 =  rmu * ( g2yi(:,2) + g1yi(:,3) )
      tmp2 =  rmu * ( g3yi(:,3) + g2yi(:,4) )
      tmp3 =  rmu * ( g1yi(:,4) + g3yi(:,2) )


      if(iconvflow.eq.2) then  ! advective form (NO IBP either)
c
c no density yet...it comes later
c
         rNa(:,1) = rNa(:,1) 
     &            + ubar(:,1) * g1yi(:,2)
     &            + ubar(:,2) * g2yi(:,2)
     &            + ubar(:,3) * g3yi(:,2)
         rNa(:,2) = rNa(:,2)
     &            + ubar(:,1) * g1yi(:,3)
     &            + ubar(:,2) * g2yi(:,3)
     &            + ubar(:,3) * g3yi(:,3)
         rNa(:,3) = rNa(:,3)
     &            + ubar(:,1) * g1yi(:,4)
     &            + ubar(:,2) * g2yi(:,4)
     &            + ubar(:,3) * g3yi(:,4)

         rGNa(:,1,1) = two * rmu * g1yi(:,2) + tmp
         rGNa(:,1,2) = tmp1
         rGNa(:,1,3) = tmp3
         rGNa(:,2,1) = tmp1
         rGNa(:,2,2) = two * rmu * g2yi(:,3) + tmp
         rGNa(:,2,3) = tmp2
         rGNa(:,3,1) = tmp3
         rGNa(:,3,2) = tmp2
         rGNa(:,3,3) = two * rmu * g3yi(:,4) + tmp
      else   ! conservative form (with IBP)

c                                            IBP conservative convection
c                                                      |||||
c                                                      vvvvv
         rGNa(:,1,1) = two * rmu * g1yi(:,2) + tmp - u1(:)*u1(:)*rho(:)
         rGNa(:,1,2) = tmp1                        - u1(:)*u2(:)*rho(:)
         rGNa(:,1,3) = tmp3                        - u1(:)*u3(:)*rho(:)
         rGNa(:,2,1) = tmp1                        - u1(:)*u2(:)*rho(:)
         rGNa(:,2,2) = two * rmu * g2yi(:,3) + tmp - u2(:)*u2(:)*rho(:)
         rGNa(:,2,3) = tmp2                        - u3(:)*u2(:)*rho(:)
         rGNa(:,3,1) = tmp3                        - u1(:)*u3(:)*rho(:)
         rGNa(:,3,2) = tmp2                        - u3(:)*u2(:)*rho(:)
         rGNa(:,3,3) = two * rmu * g3yi(:,4) + tmp - u3(:)*u3(:)*rho(:)
      endif
      if((iLES.gt.10).and.(iLES.lt.20)) then    ! bard
         rGNa(:,1,1) = rGNa(:,1,1) - rlsli(:,1)*rho(:)
         rGNa(:,1,2) = rGNa(:,1,2) - rlsli(:,4)*rho(:)
         rGNa(:,1,3) = rGNa(:,1,3) - rlsli(:,5)*rho(:)
         rGNa(:,2,1) = rGNa(:,2,1) - rlsli(:,4)*rho(:)
         rGNa(:,2,2) = rGNa(:,2,2) - rlsli(:,2)*rho(:)
         rGNa(:,2,3) = rGNa(:,2,3) - rlsli(:,6)*rho(:)
         rGNa(:,3,1) = rGNa(:,3,1) - rlsli(:,5)*rho(:)
         rGNa(:,3,2) = rGNa(:,3,2) - rlsli(:,6)*rho(:)
         rGNa(:,3,3) = rGNa(:,3,3) - rlsli(:,3)*rho(:)
      endif  
   
      tmp1        = tauM * rLui(:,1) 
      tmp2        = tauM * rLui(:,2) 
      tmp3        = tauM * rLui(:,3)
      
      rGNa(:,1,1) = rGNa(:,1,1) + tmp1 * u1 
      rGNa(:,1,2) = rGNa(:,1,2) + tmp1 * u2
      rGNa(:,1,3) = rGNa(:,1,3) + tmp1 * u3
      rGNa(:,2,1) = rGNa(:,2,1) + tmp2 * u1
      rGNa(:,2,2) = rGNa(:,2,2) + tmp2 * u2
      rGNa(:,2,3) = rGNa(:,2,3) + tmp2 * u3
      rGNa(:,3,1) = rGNa(:,3,1) + tmp3 * u1
      rGNa(:,3,2) = rGNa(:,3,2) + tmp3 * u2
      rGNa(:,3,3) = rGNa(:,3,3) + tmp3 * u3

      if(iconvflow.eq.1) then  
c
c... get the u_j w_{i,i} term in there to match A_j^T w_{i,j} tau L_i
c    to match the SUPG of incompressible limit
c 
         rGNa(:,1,1) = rGNa(:,1,1) + tmp1 * u1
         rGNa(:,1,2) = rGNa(:,1,2) + tmp2 * u1
         rGNa(:,1,3) = rGNa(:,1,3) + tmp3 * u1
         rGNa(:,2,1) = rGNa(:,2,1) + tmp1 * u2
         rGNa(:,2,2) = rGNa(:,2,2) + tmp2 * u2
         rGNa(:,2,3) = rGNa(:,2,3) + tmp3 * u2
         rGNa(:,3,1) = rGNa(:,3,1) + tmp1 * u3
         rGNa(:,3,2) = rGNa(:,3,2) + tmp2 * u3
         rGNa(:,3,3) = rGNa(:,3,3) + tmp3 * u3
      endif

      if(iconvflow.eq.2) then  ! advective form has a taubar term to restore con
         tmp1 = tauBar
     &     * ( rLui(:,1) * g1yi(:,2)
     &       + rLui(:,2) * g2yi(:,2)
     &       + rLui(:,3) * g3yi(:,2) )
         tmp2 = tauBar
     &     * ( rLui(:,1) * g1yi(:,3)
     &       + rLui(:,2) * g2yi(:,3)
     &       + rLui(:,3) * g3yi(:,3) )
         tmp3 = tauBar
     &     * ( rLui(:,1) * g1yi(:,4)
     &       + rLui(:,2) * g2yi(:,4)
     &       + rLui(:,3) * g3yi(:,4) )

         rGNa(:,1,1) = rGNa(:,1,1) + tmp1 * rLui(:,1)
         rGNa(:,1,2) = rGNa(:,1,2) + tmp1 * rLui(:,2)
         rGNa(:,1,3) = rGNa(:,1,3) + tmp1 * rLui(:,3)
         rGNa(:,2,1) = rGNa(:,2,1) + tmp2 * rLui(:,1)
         rGNa(:,2,2) = rGNa(:,2,2) + tmp2 * rLui(:,2)
         rGNa(:,2,3) = rGNa(:,2,3) + tmp2 * rLui(:,3)
         rGNa(:,3,1) = rGNa(:,3,1) + tmp3 * rLui(:,1)
         rGNa(:,3,2) = rGNa(:,3,2) + tmp3 * rLui(:,2)
         rGNa(:,3,3) = rGNa(:,3,3) + tmp3 * rLui(:,3)
      endif   ! end of advective form
c
c.... everything that gets multiplied by rNa was supposed
c     to have density multiplying it.  Do it now.

      rNa(:,1) = rNa(:,1) * rho
      rNa(:,2) = rNa(:,2) * rho
      rNa(:,3) = rNa(:,3) * rho

c
c
c.... multiply the residual pieces by the weight space
c
      do aa = 1,blk%s
c
c.... continuity
c
         rl(1:blk%e,aa,4) = rl(1:blk%e,aa,4) + WdetJ
     &              * ( shg(1:blk%e,aa,1) * uBar(1:blk%e,1) 
     &                + shg(1:blk%e,aa,2) * uBar(1:blk%e,2) 
     &                + shg(1:blk%e,aa,3) * uBar(1:blk%e,3) )
c
c.... momentum
c
         rl(1:blk%e,aa,1) = rl(1:blk%e,aa,1) - WdetJ
     &              * ( shpfun(1:blk%e,aa) * rNa(1:blk%e,1)
     &                + shg(1:blk%e,aa,1) * rGNa(1:blk%e,1,1)
     &                + shg(1:blk%e,aa,2) * rGNa(1:blk%e,1,2)
     &                + shg(1:blk%e,aa,3) * rGNa(1:blk%e,1,3) )
         rl(1:blk%e,aa,2) = rl(1:blk%e,aa,2) - WdetJ
     &              * ( shpfun(1:blk%e,aa) * rNa(1:blk%e,2)
     &                + shg(1:blk%e,aa,1) * rGNa(1:blk%e,2,1)
     &                + shg(1:blk%e,aa,2) * rGNa(1:blk%e,2,2)
     &                + shg(1:blk%e,aa,3) * rGNa(1:blk%e,2,3) )
         rl(1:blk%e,aa,3) = rl(1:blk%e,aa,3) - WdetJ
     &              * ( shpfun(1:blk%e,aa) * rNa(1:blk%e,3)
     &                + shg(1:blk%e,aa,1) * rGNa(1:blk%e,3,1)
     &                + shg(1:blk%e,aa,2) * rGNa(1:blk%e,3,2)
     &                + shg(1:blk%e,aa,3) * rGNa(1:blk%e,3,3) )
      
      enddo                 
c
c.... return
c
      return
      end



c------------------------------------------------------------------------
c
c     calculate the residual for the advection-diffusion equation
c
c------------------------------------------------------------------------
      subroutine e3ResSclr ( blk, uMod,              gGradS,
     &                       Sclr,		Sdot,	gradS,  
     &                       WdetJ,		rLS,	tauS,
     &                       shpfun,            shg,    src,
     &                       diffus,
     &                       rl )
c
      use eblock
      include "common.h"
      type (LocalBlkData) blk

      real*8    uMod(blk%e,nsd),   gGradS(blk%e, nsd),
     &          Sclr(blk%e),       Sdot(blk%e),	gradS(blk%e,nsd),
     &          WdetJ(blk%e),      rLS(blk%e),	rho(blk%e),
     &          tauS(blk%e),       shpfun(blk%e,blk%s), src(blk%e), 
     &          shg(blk%e,blk%s,3), rl(bsz,blk%s)
      
      real*8    diffus(blk%e)
c
c.... local declarations
c
      real*8    rGNa(blk%e,nsd),   rNa(blk%e),  rcp(blk%e), tmp(blk%e)

      integer   aa
c     
c.... initialize multipliers for Na and Na_{,i}
c
      rNa  = zero
      rGNa = zero
c
c.... Na multiplier
c
      tmps     = one-flmpr  ! consistant mass factor
      rcp = one ! rho * cp
      

         rNa = rcp*(tmps*Sdot + uMod(:,1) * gradS(:,1)
     &                         + uMod(:,2) * gradS(:,2)
     &                         + uMod(:,3) * gradS(:,3) ) 
     &        - src



      tmp = rcp * tauS * (rLS -src)
c
c.... Na,i multiplier
c
      rGNa(:,1) = diffus * gradS(:,1) + uMod(:,1) * tmp
      rGNa(:,2) = diffus * gradS(:,2) + uMod(:,2) * tmp
      rGNa(:,3) = diffus * gradS(:,3) + uMod(:,3) * tmp
c
      if (idcsclr(1) .ne. 0) then
         if ((idcsclr(2).eq.1 .and. isclr.eq.1) .or. 
     &        (idcsclr(2).eq.2 .and. isclr.eq.2)) then ! scalar with dc
c
c.... add the contribution of DC to residual
c
            rGNa(:,1) = rGNa(:,1) + gGradS(:,1) ! gGradS is 
            rGNa(:,2) = rGNa(:,2) + gGradS(:,2) ! g^{ij}*Y_{j}*dcFct
            rGNa(:,3) = rGNa(:,3) + gGradS(:,3) ! calculated in e3dc.f
c
         endif
      endif                     ! end of idcsclr
c
c.... multiply the residual pieces by the weight space
c
      do aa = 1,blk%s
c
         rl(1:blk%e,aa) = rl(1:blk%e,aa)	- WdetJ
     &                        * ( shpfun(:,aa) * rNa(:)
     &                        + shg(:,aa,1) * rGNa(:,1)
     &                        + shg(:,aa,2) * rGNa(:,2)
     &                        + shg(:,aa,3) * rGNa(:,3) )

      enddo
c
c.... return
c
      return
      end



c----------------------------------------------------------------------
c     Calculate the strong PDE residual.
c----------------------------------------------------------------------
      subroutine e3resStrongPDE( blk,
     &     aci,  u1,   u2,   u3,   Temp, rho,  xx,
     &           g1yi, g2yi, g3yi,
     &     rLui, src, divqi)
      use eblock
      include "common.h"
      type (LocalBlkData) blk

c     INPUTS
      double precision, intent(in), dimension(blk%e,nsd) :: 
     &     aci, xx
      double precision, intent(in), dimension(blk%e,nflow) :: 
     &     g1yi, g2yi, g3yi
      double precision, intent(in), dimension(blk%e) ::
     &     u1, u2, u3, Temp, rho
c     OUTPUTS
      double precision, intent(out), dimension(blk%e,nsd) ::
     &     rLui, src
c     LOCALS
      double precision, dimension(blk%e) ::
     &     divu
      double precision, dimension(blk%e,nsd) ::
     &     divqi
      double precision, dimension(nsd) ::
     &     omega
c.... compute source term
      src = zero
      if(matflg(5,1) .ge. 1) then
c        body force contribution to src
         bfx      = datmat(1,5,1) ! Boussinesq, g*alfap
         bfy      = datmat(2,5,1)
         bfz      = datmat(3,5,1)
         select case ( matflg(5,1) )
         case ( 1 )             ! standard linear body force
            src(:,1) = bfx
            src(:,2) = bfy
            src(:,3) = bfz
         case ( 2 )             ! boussinesq body force
            Tref = datmat(2,2,1)
            src(:,1) = bfx * (Temp(:)-Tref)
            src(:,2) = bfy * (Temp(:)-Tref)
            src(:,3) = bfz * (Temp(:)-Tref)
         case ( 3 )             ! user specified f(x,y,z)
            call e3source(xx, src)
         end select
      endif
c     
      if(matflg(6,1).eq.1) then
c        coriolis force contribution to src
         omega(1)=datmat(1,6,1)
         omega(2)=datmat(2,6,1)
         omega(3)=datmat(3,6,1)
c  note that we calculate f as if it contains the usual source
c  plus the Coriolis and the centrifugal forces taken to the rhs (sign change)
c  as long as we are doing SUPG with no accounting for these terms in the
c  LHS this is the only change (which will find its way to the RHS momentum
c  equation (both Galerkin and SUPG parts)).
c
c  uncomment later if you want rotation always about z axis
c                 orig_src - om x om x r       - two om x u
c
c$$$          src(:,1)=src(:,1)+omega(3)*omega(3)*xx(:,1)+two*omega(3)*u2
c$$$          src(:,2)=src(:,2)+omega(3)*omega(3)*xx(:,2)-two*omega(3)*u1
c
c more general for testing
c
         src(:,1)=src(:,1)
     &        -omega(2)*(omega(1)*xx(:,2)-omega(2)*xx(:,1))
     &        -omega(3)*(omega(1)*xx(:,3)-omega(3)*xx(:,1))
     &        -two*(omega(2)*u3-omega(3)*u2)
         src(:,2)=src(:,2)
     &        -omega(1)*(omega(2)*xx(:,1)-omega(1)*xx(:,2))
     &        -omega(3)*(omega(2)*xx(:,3)-omega(3)*xx(:,2))
     &        -two*(omega(3)*u1-omega(1)*u3)
         src(:,3)=src(:,3)
     &        -omega(1)*(omega(3)*xx(:,1)-omega(1)*xx(:,3))
     &        -omega(2)*(omega(3)*xx(:,2)-omega(2)*xx(:,3))
     &        -two*(omega(1)*u2-omega(2)*u1)
      endif
c     calculate momentum residual
      rLui(:,1) =(aci(:,1) + u1 * g1yi(:,2)
     &     + u2 * g2yi(:,2)
     &     + u3 * g3yi(:,2) - src(:,1) ) * rho
     &     + g1yi(:,1)
     &        - divqi(:,1)
      rLui(:,2) =(aci(:,2) + u1 * g1yi(:,3)
     &     + u2 * g2yi(:,3)
     &     + u3 * g3yi(:,3) - src(:,2) ) * rho
     &     + g2yi(:,1)
     &        - divqi(:,2)
      rLui(:,3) =(aci(:,3) + u1 * g1yi(:,4)
     &     + u2 * g2yi(:,4)
     &     + u3 * g3yi(:,4) - src(:,3) ) * rho
     &     + g3yi(:,1)
     &        - divqi(:,3)
      if(iconvflow.eq.1) then
         divu(:)  = (g1yi(:,2) + g2yi(:,3) + g3yi(:,4))*rho
         rLui(:,1)=rlui(:,1)+u1(:)*divu(:)
         rLui(:,2)=rlui(:,2)+u2(:)*divu(:)
         rLui(:,3)=rlui(:,3)+u3(:)*divu(:)
      endif
c
      return
      end subroutine e3resStrongPDE
