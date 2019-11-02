      subroutine e3LHS ( blk, u1,        u2,         u3,
     &                   uBar,      WdetJ,      rho,
     &                   rLui,      rmu,       
     &                   tauC,      tauM,       tauBar,
     &                   shpfun,    shg,        xlhs)
c------------------------------------------------------------------------
c 
c  This routine computes the left hand side tangent matrix at an 
c  integration point.
c
c  input:
c     u1(blk%e)                  : x1-velocity
c     u2(blk%e)                  : x2-velocity
c     u3(blk%e)                  : x3-velocity
c     uBar(blk%e,3)              : u - tauM * Li
c     WdetJ(blk%e)               : weighted jacobian determinant
c     rLui(blk%e,3)              : total residual of NS equations
c     rmu(blk%e)                 : fluid viscosity
c     rho(blk%e)				  : fluid density
c     tauC(blk%e)                : continuity tau
c     tauM(blk%e)                : momentum tau
c     tauBar(blk%e)              : additional tau
c     shpfun(blk%e,blk%s)         : element shpfun functions
c     shg(blk%e,blk%s,3)          : global grad of element shape functions
c
c  output:
c     xlhs(blk%e,16,blk%s,blk%s) : left hand side
c
c
c------------------------------------------------------------------------
      use eblock
      include "common.h"
      type (LocalBlkData) blk


      dimension u1(blk%e),         u2(blk%e),       u3(blk%e),
     &          uBar(blk%e,3),     WdetJ(blk%e),    rho(blk%e),
     &          rLui(blk%e,3),     rmu(blk%e),   
     &          tauC(blk%e),       tauM(blk%e),     tauBar(blk%e),
     &          shpfun(blk%e,blk%s),shg(blk%e,blk%s,3)
      
#ifdef SP_LHS
      real*4 xlhs(bsz,16,blk%s,blk%s)
#else
      real*8 xlhs(bsz,16,blk%s,blk%s)
#endif
c
c.... local declarations
c
      dimension t1(blk%e,3),       t2(blk%e,3),      t3(blk%e,3),
     &          tmp1(blk%e),       tmp2(blk%e),    
     &          tmp(blk%e),        tlW(blk%e)

      integer   aa, b
      
      real*8    lhmFct, lhsFct,           tsFct(blk%e)
      
      lhsFct = alfi * gami * Delt(itseq)
      lhmFct = almi * (one - flmpl) 
c
c.... scale variables for efficiency
c
      tlW      = lhsFct * WdetJ     
      tmp1      = tlW * rho
      tauM      = tlW * tauM 
      tauC      = tlW * tauC 
      rmu       = tlW * rmu 
      tsFct     = lhmFct * WdetJ * rho
      if(iconvflow.eq.2) then  ! 2 is ubar form 3 is cons form but ubar tang. 
         tauBar    = lhsFct * WdetJ * tauBar 
         uBar(:,1) = tmp1 * uBar(:,1)
         uBar(:,2) = tmp1 * uBar(:,2)
         uBar(:,3) = tmp1 * uBar(:,3)
      else
         tauBar = zero  !lazy tangent...if effective code it
         uBar(:,1) = tmp1 * u1(:)
         uBar(:,2) = tmp1 * u2(:)
         uBar(:,3) = tmp1 * u3(:)
      endif

c
c.... compute mass and convection terms
c
      do b = 1, blk%s
         t1(:,1) = uBar(:,1) * shg(:,b,1)
     &           + uBar(:,2) * shg(:,b,2)
     &           + uBar(:,3) * shg(:,b,3)
c
c t1=ubar_k N^b,k*rho*alpha_f*gamma*deltat*WdetJ  
c
         do aa = 1, blk%s
            tmp1 = tsFct * shpfun(:,aa) * shpfun(:,b)
            tmp2 = tmp1 + t1(:,1) * shpfun(:,aa)
c
c tmp1=alpha_m*(1-lmp)*WdetJ*N^aN^b*rho   the time term CORRECT
c tmp2=tmp1+N^a*ubar_k N^b,k*rho*alpha_f*gamma*deltat*WdetJ   the 
c    second term is convective term CORRECT
c            
            xlhs(1:blk%e,1,aa,b) = xlhs(1:blk%e,1,aa,b) + tmp2
            xlhs(1:blk%e,6,aa,b) = xlhs(1:blk%e,6,aa,b) + tmp2
            xlhs(1:blk%e,11,aa,b) = xlhs(1:blk%e,11,aa,b) + tmp2
         enddo
      enddo
c
c.... compute the rest of K (symmetric terms)
c      
      do b = 1, blk%s
         
         t1(:,1) = tauC * shg(:,b,1)
         t1(:,2) = tauC * shg(:,b,2)
         t1(:,3) = tauC * shg(:,b,3)

c t1 is tauC*N^b_i,j*alpha_f*gamma*deltat*WdetJ
         
         t2(:,1) = rmu  * shg(:,b,1)
         t2(:,2) = rmu  * shg(:,b,2)
         t2(:,3) = rmu  * shg(:,b,3)
c t2 is mu*N^b_j,k*alpha_f*gamma*deltat*WdetJ
      
         tmp1 = tauM   * ( u1 * shg(:,b,1)  
     &                   + u2 * shg(:,b,2) 
     &                   + u3 * shg(:,b,3) )*rho
c tmp1 is tauM*(rho u_m N^b_j,m)*alpha_f*gamma*deltat*WdetJ

         tmp2 = tauBar * ( rLui(:,1) * shg(:,b,1)
     &                   + rLui(:,2) * shg(:,b,2)
     &                   + rLui(:,3) * shg(:,b,3) )
c tmp2 is taubar*(L_m N^b_j,m)*alpha_f*gamma*deltat*WdetJ
         t3(:,1) = t2(:,1) + tmp1 * u1 + tmp2 * rLui(:,1)
         t3(:,2) = t2(:,2) + tmp1 * u2 + tmp2 * rLui(:,2)
         t3(:,3) = t2(:,3) + tmp1 * u3 + tmp2 * rLui(:,3)

c t3 is   (mu*N^b_j,k + u_k tauM*(rho u_m N^b_j,m)+ L_k*taubar*(L_mN^b_j,m ) 
c   *alpha_f*gamma*deltat*WdetJ     which isline 2 page 40 of whiting
c   ALMOST (waiting to get hit with N^a_{i,k}
c mu correct NOW (wrong before) and rho weight on tauM term
c
c.... first do the (nodal) diagonal blocks         
c
         aa  = b
         
         tmp = t3(:,1) * shg(:,aa,1)
     &       + t3(:,2) * shg(:,aa,2)
     &       + t3(:,3) * shg(:,aa,3)
c previous command is the N^a_{i,k} dot product with t3 defined above

         xlhs(1:blk%e,1,aa,b) = xlhs(1:blk%e,1,aa,b) + tmp
     &                      + t1(1:blk%e,1) * shg(1:blk%e,aa,1)
     &                      + t2(1:blk%e,1) * shg(1:blk%e,aa,1)
         xlhs(1:blk%e,6,aa,b) = xlhs(1:blk%e,6,aa,b) + tmp
     &                      + t1(1:blk%e,2) * shg(1:blk%e,aa,2)
     &                      + t2(1:blk%e,2) * shg(1:blk%e,aa,2)
         xlhs(1:blk%e,11,aa,b) = xlhs(1:blk%e,11,aa,b) + tmp
     &                      + t1(1:blk%e,3) * shg(1:blk%e,aa,3)
     &                      + t2(1:blk%e,3) * shg(1:blk%e,aa,3)
c
         tmp1               = t1(:,1) * shg(:,aa,2)
     &                      + t2(:,2) * shg(:,aa,1)
         xlhs(1:blk%e,2,aa,b) = xlhs(1:blk%e,2,aa,b) + tmp1 
         xlhs(1:blk%e,5,b,aa) = xlhs(1:blk%e,5,b,aa) + tmp1 
c
         tmp1               = t1(:,1) * shg(:,aa,3)
     &                      + t2(:,3) * shg(:,aa,1)
         xlhs(1:blk%e,3,aa,b) = xlhs(1:blk%e,3,aa,b) + tmp1 
         xlhs(1:blk%e,9,b,aa) = xlhs(1:blk%e,9,b,aa) + tmp1 
c
         tmp1               = t1(:,2) * shg(:,aa,3)
     &                      + t2(:,3) * shg(:,aa,2)
         xlhs(1:blk%e,7,aa,b) = xlhs(1:blk%e,7,aa,b) + tmp1 
         xlhs(1:blk%e,10,b,aa) = xlhs(1:blk%e,10,b,aa) + tmp1 
c
c.... now the off-diagonal (nodal) blocks
c
         do aa = b+1, blk%s
            tmp             = t3(:,1) * shg(:,aa,1)
     &                      + t3(:,2) * shg(:,aa,2)
     &                      + t3(:,3) * shg(:,aa,3)
c
            tmp1            = tmp
     &                      + t1(:,1) * shg(:,aa,1)
     &                      + t2(:,1) * shg(:,aa,1)
            xlhs(1:blk%e,1,aa,b) = xlhs(1:blk%e,1,aa,b) + tmp1
            xlhs(1:blk%e,1,b,aa) = xlhs(1:blk%e,1,b,aa) + tmp1
c
            tmp1            = tmp
     &                      + t1(1:blk%e,2) * shg(:,aa,2)
     &                      + t2(:,2) * shg(:,aa,2)
            xlhs(1:blk%e,6,aa,b) = xlhs(1:blk%e,6,aa,b) + tmp1
            xlhs(1:blk%e,6,b,aa) = xlhs(1:blk%e,6,b,aa) + tmp1
c
            tmp1            = tmp
     &                      + t1(:,3) * shg(:,aa,3)
     &                      + t2(:,3) * shg(:,aa,3)
            xlhs(1:blk%e,11,aa,b) = xlhs(1:blk%e,11,aa,b) + tmp1
            xlhs(1:blk%e,11,b,aa) = xlhs(1:blk%e,11,b,aa) + tmp1
c
c.... ( i != j )
c
            tmp1               = t1(:,1) * shg(:,aa,2)
     &                         + t2(:,2) * shg(:,aa,1)
            xlhs(1:blk%e,2,aa,b) = xlhs(1:blk%e,2,aa,b) + tmp1
            xlhs(1:blk%e,5,b,aa) = xlhs(1:blk%e,5,b,aa) + tmp1
c
            tmp1               = t1(:,1) * shg(:,aa,3)
     &                         + t2(:,3) * shg(:,aa,1)
            xlhs(1:blk%e,3,aa,b) = xlhs(1:blk%e,3,aa,b) + tmp1
            xlhs(1:blk%e,9,b,aa) = xlhs(1:blk%e,9,b,aa) + tmp1
c
            tmp1               = t1(:,2) * shg(:,aa,1)
     &                         + t2(:,1) * shg(:,aa,2)
            xlhs(1:blk%e,5,aa,b) = xlhs(1:blk%e,5,aa,b) + tmp1
            xlhs(1:blk%e,2,b,aa) = xlhs(1:blk%e,2,b,aa) + tmp1
c
            tmp1               = t1(:,2) * shg(:,aa,3)
     &                         + t2(:,3) * shg(:,aa,2)
            xlhs(1:blk%e,7,aa,b) = xlhs(1:blk%e,7,aa,b) + tmp1
            xlhs(1:blk%e,10,b,aa) = xlhs(1:blk%e,10,b,aa) + tmp1
c
            tmp1               = t1(:,3) * shg(:,aa,1)
     &                         + t2(:,1) * shg(:,aa,3)
            xlhs(1:blk%e,9,aa,b) = xlhs(1:blk%e,9,aa,b) + tmp1
            xlhs(1:blk%e,3,b,aa) = xlhs(1:blk%e,3,b,aa) + tmp1
c
            tmp1               = t1(:,3) * shg(:,aa,2)
     &                         + t2(:,2) * shg(:,aa,3)
            xlhs(1:blk%e,10,aa,b) = xlhs(1:blk%e,10,aa,b) + tmp1
            xlhs(1:blk%e,7,b,aa) = xlhs(1:blk%e,7,b,aa) + tmp1
c
         enddo
      enddo
c
c.... compute G   Nai Nbp,j
c
      
      do b = 1, blk%s
         t1(:,1) = tlW * shg(:,b,1)
         t1(:,2) = tlW * shg(:,b,2)
         t1(:,3) = tlW * shg(:,b,3)
         do aa = 1, blk%s
            xlhs(1:blk%e,4,aa,b) = xlhs(1:blk%e,4,aa,b) 
     &      + t1(1:blk%e,1) * shpfun(1:blk%e,aa)  
            xlhs(1:blk%e,8,aa,b) = xlhs(1:blk%e,8,aa,b) 
     &      + t1(1:blk%e,2) * shpfun(1:blk%e,aa)  
            xlhs(1:blk%e,12,aa,b) = xlhs(1:blk%e,12,aa,b) 
     &      + t1(1:blk%e,3) * shpfun(1:blk%e,aa)  
         enddo
      enddo
c
c.... compute C
c we divide by rho because the L on the weight space is density divided
c      form
c
      tauM=tauM/rho
      do b = 1, blk%s
         t1(:,1) = tauM * shg(:,b,1)
         t1(:,2) = tauM * shg(:,b,2)
         t1(:,3) = tauM * shg(:,b,3)
         do aa = b, blk%s ! note loops starts at b;1:b-1 after quadrature loop in e3 
            xlhs(1:blk%e,16,aa,b) = xlhs(1:blk%e,16,aa,b) 
     &                      + t1(1:blk%e,1) * shg(1:blk%e,aa,1)
     &                      + t1(1:blk%e,2) * shg(1:blk%e,aa,2)
     &                      + t1(1:blk%e,3) * shg(1:blk%e,aa,3)
         enddo
      enddo

c
c.... return
c
      return
      end


c------------------------------------------------------------------------
c
c     calculate the tangent matrix for the advection-diffusion equation
c
c------------------------------------------------------------------------
      subroutine e3LHSSclr ( blk,uMod,      giju,       dcFct,
     &                       Sclr,      Sdot,       gradS,  
     &                       WdetJ,     rLS,        tauS,
     &                       shpfun,    shg,        srcL,
     &                       diffus,
     &                       xSebe )

c
      use eblock
      include "common.h"
      type (LocalBlkData) blk


      real*8    uMod(blk%e,nsd),
     &          Sclr(blk%e),       Sdot(blk%e),   gradS(blk%e,nsd),
     &          WdetJ(blk%e),      rLS(blk%e),        rho(blk%e), 
     &          tauS(blk%e),       shpfun(blk%e,blk%s),  
     &          srcL(blk%e),        shg(blk%e,blk%s,3)
#ifdef SP_LHS
      real*4    xSebe(blk%e,blk%s,blk%s)
#else
      real*8    xSebe(blk%e,blk%s,blk%s)
#endif

      
      real*8    diffus(blk%e),  cp,  kptmp(blk%e),tauSo(blk%e)

c
c.... local declarations
c
      dimension t1(blk%e,3),       tmp1(blk%e),       tmp2(blk%e),
     &          tmp(blk%e),        dcFct(blk%e),      giju(blk%e,6)

      integer   aa, b
      
      real*8    lhsFct,           tsFct(blk%e)
      
      lhsFct = alfi * gami * Delt(itseq)
c
c.... scale variables for efficiency
c     
      tauSo     = tauS
      tauS      = lhsFct * WdetJ * tauS 
      kptmp     = lhsFct * WdetJ * diffus
      tsFct     = almi   * WdetJ * (one - flmpl)
      srcL       = srcL    * WdetJ * lhsFct
c
c.... compute mass and convection terms
c
      do b = 1, blk%s
         t1(:,1) = WdetJ * ( uMod(:,1) * shg(:,b,1)
     &                     + uMod(:,2) * shg(:,b,2)
     &                     + uMod(:,3) * shg(:,b,3) )
         t1(:,2) = t1(:,1) * tauSo
         do aa = 1, blk%s
            tmp1 = shpfun(:,aa) * shpfun(:,b)
            tmp2 = shpfun(:,aa) * lhsFct
            xSebe(1:blk%e,aa,b) = xSebe(1:blk%e,aa,b) + tmp1 * (tsFct + srcL)
     &                                    + tmp2 * t1(1:blk%e,1)
c
c.... compute mass term for stab u_j N_{a,j} tau N_b (note that a and b
c            flipped on both sides below)
c
            xSebe(1:blk%e,b,aa) = xSebe(1:blk%e,b,aa) + t1(:,2)*shpfun(:,aa)*
     &                                      almi*(one-flmpl)
         enddo
      enddo
c
c.... compute the rest of S (symmetric terms)
c      
      do b = 1, blk%s
         tmp     = tauS(:) 
     &             * ( uMod(:,1) * shg(:,b,1)
     &               + uMod(:,2) * shg(:,b,2)
     &               + uMod(:,3) * shg(:,b,3) )

         t1(:,1) = kptmp * shg(:,b,1) + uMod(:,1) * tmp
         t1(:,2) = kptmp * shg(:,b,2) + uMod(:,2) * tmp
         t1(:,3) = kptmp * shg(:,b,3) + uMod(:,3) * tmp
         if (idcsclr(1) .ne. 0) then
            if ((idcsclr(2).eq.1 .and. isclr.eq.1) .or. 
     &          (idcsclr(2).eq.2 .and. isclr.eq.2) .or.
     &          (idcsclr(2).eq.3 .and. isclr.eq.1) .or.
     &           (idcsclr(2).eq.3 .and. isclr.eq.2)) then ! scalar with dc
c
               tmp = WdetJ * dcFct * lhsFct
c
               giju(:,1)	= tmp * giju(:,1)
               giju(:,2)	= tmp * giju(:,2)
               giju(:,3)	= tmp * giju(:,3)
               giju(:,4)	= tmp * giju(:,4)
               giju(:,5)	= tmp * giju(:,5)
               giju(:,6)	= tmp * giju(:,6)
c       
               t1(:,1) = t1(:,1) + giju(:,1) * shg(:,b,1) 
     2                           + giju(:,4) * shg(:,b,2) 
     3			         + giju(:,6) * shg(:,b,3)
               t1(:,2) = t1(:,2) + giju(:,4) * shg(:,b,1) 
     2                           + giju(:,2) * shg(:,b,2) 
     3			         + giju(:,5) * shg(:,b,3)
               t1(:,3) = t1(:,3) + giju(:,6) * shg(:,b,1) 
     2                           + giju(:,5) * shg(:,b,2) 
     3			         + giju(:,3) * shg(:,b,3)
            endif
         endif                  !end of idcsclr
c
c.... first do the (nodal) diagonal blocks         
c
         aa  = b
         
         xSebe(1:blk%e,aa,b) = xSebe(1:blk%e,aa,b) + t1(1:blk%e,1) * shg(1:blk%e,aa,1)
     &                                 + t1(1:blk%e,2) * shg(1:blk%e,aa,2)
     &                                 + t1(1:blk%e,3) * shg(1:blk%e,aa,3)

c
c.... now the off-diagonal (nodal) blocks
c
         do aa = b+1, blk%s
            tmp             = t1(:,1) * shg(:,aa,1)
     &                      + t1(:,2) * shg(:,aa,2)
     &                      + t1(:,3) * shg(:,aa,3)
c
            xSebe(1:blk%e,aa,b) = xSebe(1:blk%e,aa,b) + tmp
            xSebe(1:blk%e,b,aa) = xSebe(1:blk%e,b,aa) + tmp
c
         enddo
      enddo
      
c
c.... return
c
      return
      end

