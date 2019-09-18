      subroutine elm3SST (blk, yl, xl, shg, kay , omega ,
     &     dwall, dwl, gradV,
     &     srcRat1 , src1 , srcJac, add, IDDESfun )
     
      use turbKW
      use eblock
      use spanStats
      include "common.h"
      type (LocalBlkData) blk

      real*8 kay(blk%e),    omega(blk%e),
     &       dwall(blk%e),  gradV(blk%e,nsd,nsd)
      real*8 srcRat1(blk%e) ,  src1(blk%e) , srcJac(blk%e,4)
      real*8 yl(bsz,blk%s,ndof),    shp(blk%e,blk%s),
     &       shg(blk%e,blk%s,nsd),   xl(bsz,blk%n,nsd),
     &       dxidx(blk%e,nsd,nsd),   WdetJ(blk%e),
     &       gradK(blk%e,nsd),    gradW(blk%e,nsd), !!check v
     &       add(blk%e,nsd),  dwl(bsz,blk%n)
      real*8 rmu(blk%e), rho(blk%e)
      real*8 srcRat(blk%e,2), src(blk%e,2)
      integer advdiff, e, n ! n is possibly redundant. clean up required
      real*8 k, kq, omg, y, mut, mut_k,
     &       mut_w, nu, F1, F2, P1, F2Inv, Rou, mu
      real*8 k_Prod(blk%e), k_Dest(blk%e), k_Rat(blk%e)
      real*8 w_Prod(blk%e), w_Dest(blk%e), w_CD(blk%e), w_Rat(blk%e)
      real*8 F1array(blk%e), F2array(blk%e)
      real*8 alpha, theta, chi, psi, ss, ssq, ssqInv
      real*8 alp, beta, sigk, sigw, omicron
      real*8 Sk, Sk_k, Sw, Sw_w, delKdelW
      real*8 kInv,omgInv, mutInv, den
      real*8 Pk1, Pk2, Pk, Pk_max, Pkrat, Pkrat_k, Pk_k, Pk1_k, Pk2_k
      real*8 Dk, Dkrat, Dkrat_k, Dk_k, tmp
      real*8 Pw, Pwrat, Pw_w, Pwrat_w
      real*8 Dw, Dwrat, Dw_w, Dwrat_w
      real*8 CDw, CDw_w
      real*8 Sk_sust, Sw_sust, k_amb, Omega_amb
      real*8 yloc, ysum
c     DDES and IDDES variables
      real*8 dx, dy, dz, hmax, CDES, lLES, lRANS, qfac, dwallsqqfact,
     &       rd, fd, lDDES, hwn, delta, rdt, rdl, fl, ft, fe2, alphaI,
     &       fe1, fe, fbarg, fb, fdt, lIDDES, kP3, lIDDESInv, xcenter,
     &       cf, utau, retau
      real*8 IDDESfun(blk%e,nfun) 
c     Density and molecular viscosity
      rho(:)=datmat(1,1,1)
      rmu(:)=datmat(1,2,1)

c     Debug flag to do only advection diffusion problem
      advdiff = 0

c     NOT doing advection-diffucion only, compute source terms
      if(advdiff . eq . 0) then
        call gradKWgen(blk, yl, shg, gradK, gradW)
c       Loop over elements in this block
        do e = 1, blk%e
           Rou = rho(e)
           nu = rmu(e)/Rou
           mu = rmu(e)
           k = max(zero,kay(e)) ! prevent negative values
!           k = abs(kay(e))
           kInv = one / k
           if (k.lt.1.0d-12) kInv=1.0d12
           omg = max(zero,omega(e)) ! prevent negative values
           omgInv = one / omg
           if (omg.lt.1.d-12) omgInv=1.0d12
           y = dwall(e)

c          alpha = dui/dxj*dui/dxj
           alpha = gradV(e,1,1) ** 2 + gradV(e,1,2) ** 2
     &       + gradV(e,1,3) ** 2 + gradV(e,2,1) ** 2
     &       + gradV(e,2,2) ** 2 + gradV(e,2,3) ** 2
     &       + gradV(e,3,1) ** 2 + gradV(e,3,2) ** 2
     &       + gradV(e,3,3) ** 2

c          theta = dui/dxj*duj/dxi
           theta = gradV(e,1,1) ** 2 + gradV(e,2,2) ** 2
     &       + gradV(e,3,3) ** 2 + 2*gradV(e,2,1)*gradV(e,1,2)
     &       + 2*gradV(e,3,1)*gradV(e,1,3)
     &       + 2*gradV(e,3,2)*gradV(e,2,3)

c          chi = duk/dxk*dum/dxm = div(u)^2
           if (matflg(1,1).eq.0) then !compressible
             chi = (gradV(e,1,1) + gradV(e,2,2)
     &               + gradV(e,3,3)) ** 2
           else                         !incompressible
             chi = zero
           endif

c          psi = duk/dxk = div(u) = sqrt(chi)
           if (matflg(1,1).eq.0) then
             psi = (gradV(e,1,1) + gradV(e,2,2)
     &              + gradV(e,3,3))
           else
             psi = zero
           endif

c          Shear stress squared ss=S_ij*S_ij
           ss = 0.5 * (alpha + theta)
c          Shear stress invariant ssq=sqrt(2*S_ij*S_ij)
           ssq = sqrt(two * ss)
           ssqInv = one / ssq 
           if (ssq.lt.1d-12) ssqInv = 1d12

c          Cross diffusion term dk/dx_j * dw/dx_j
           delKdelW = gradK(e,1)*gradW(e,1)  +  gradK(e,2)*gradW(e,2)
     &           +  gradK(e,3)*gradW(e,3)

c          Compute the blending functions F1 and F2
           call getblendfunc (delKdelW, k, omg, y, Rou, nu, F1, F2 )
           F1array(e) = F1
           F2array(e) = F2

c          Get the blended value of the model constants
           alp = F1 * alp1 + (one - F1) * alp2
           beta = F1 * beta1 + (one - F1) * beta2
           sigk = F1 * sigk1 + (one - F1) * sigk2
           sigw = F1 * sigw1 + (one - F1) * sigw2
           omicron = F1 * gam1 + (one - F1) * gam2 ! this is gamma in standard definition

c          Calculate eddy viscosity and its variants
           if (a1 * omg . gt . F2 * ssq) then
              mut = Rou * k * omgInv ! simple k/w
              mut = max(mut, 1.0d-4*mu) ! limit mut to 0.01% of mu from below
              mut_k = Rou *  omgInv  ! dmut/dk
              mut_w = - k * Rou * (omgInv ** 2) ! dmut/dw
           else
              den = max(F2 * ssq, 1.0d-12)
              mut = a1 * Rou * k / den
              mut = max(mut, 1.0d-4*mu) ! limit mut to 0.01% of mu from below
              mut_k = a1 * Rou / den
              mut_w = zero
           endif
           mutInv = one / mut

c          Do some computations for DDES and IDDES models if desired
           if (iLES.lt.0) then
               dx = maxval(xl(e,:,1))-minval(xl(e,:,1))
               dy = maxval(xl(e,:,2))-minval(xl(e,:,2))
               dz = maxval(xl(e,:,3))-minval(xl(e,:,3))
               if (iles.eq.-1.or.iles.eq.-2) then ! regular DES and DDES
                  hmax = max(dx,max(dy,dz))
                  CDES = 0.780d0*F1+0.610d0*(one-F1)
                  lLES = CDES*hmax
                  lRANS = sqrt(k)*omgInv/CmuKW
                  qfac = sqrt(
     &                gradV(e,1,1)**2+gradV(e,1,2)**2+gradV(e,1,3)**2
     &              + gradV(e,2,1)**2+gradV(e,2,2)**2+gradV(e,2,3)**2
     &              + gradV(e,3,1)**2+gradV(e,3,2)**2+gradV(e,3,3)**2 )
                  dwallsqqfact = max(kappa**2*dwall(e)**2*qfac,1.0d-12)
                  rd = (mut/Rou+nu)/dwallsqqfact
                  fd = one-tanh((20.0d0*rd)**3)
                  lDDES = lRANS-fd*max(zero,lRANS-lLES)
               elseif (iles.eq.-3) then ! IDDES
                  hmax = max(dx,max(dy,dz)) ! definition in IDDES paper
c                 hwn is the local step in the wall normal direction
                  hwn = maxval(dwl(e,:)) - minval(dwl(e,:))
                  delta=min(hmax,max(0.15*dwall(e),max(0.15*hmax,hwn))) !Cw=0.15
                  CDES = 0.780d0*F1+0.610d0*(one-F1)
                  lLES = CDES*delta
                  lRANS = sqrt(k)*omgInv/CmuKW
                  qfac = sqrt(
     &                gradV(e,1,1)**2+gradV(e,1,2)**2+gradV(e,1,3)**2
     &              + gradV(e,2,1)**2+gradV(e,2,2)**2+gradV(e,2,3)**2
     &              + gradV(e,3,1)**2+gradV(e,3,2)**2+gradV(e,3,3)**2 )
                  dwallsqqfact = max(kappa**2*dwall(e)**2*qfac,1.0d-12)
                  rdt = (mut/Rou)/dwallsqqfact
                  rdl = (nu)/dwallsqqfact
                  fdt = one-tanh((20.0d0*rdt)**3.0d0)
                  fl = tanh((5.0d0**2*rdl)**10.0d0)
                  ft = tanh((1.870d0**2*rdt)**3.0d0)
                  fe2 = one-max(ft,fl)
                  if (iLocInterface.eq.0) then
                    alphaI = 0.250d0-dwall(e)/hmax ! y=dwall(e) above
                    fbarg = two*exp(-9.0d0*alphaI**2)
                    fb = min(fbarg,one)
                    fd = max((one-fdt),fb)
                  elseif (iLocInterface.eq.1) then
                    xcenter = (minval(xl(e,:,1))+maxval(xl(e,:,1)))*0.50
                    tmp = 0.01665*xcenter+0.16158 ! delta_995(x) for CRS flat plate with IDDES-SA
                    tmp = tmp/10
                    alphaI = 0.250d0-dwall(e)/tmp
                    fbarg = two*exp(-9.0d0*alphaI**2)
                    fb = min(fbarg,one)
                    fd = fb
                  elseif (iLocInterface.eq.2) then
                    xcenter = (minval(xl(e,:,1))+maxval(xl(e,:,1)))*0.50
                    cf =  0.000016954*xcenter**2 ! Cf(x) for CRS flat plate with IDDES-SA
     &                   - 0.000090025*xcenter+0.0038619
                    utau = sqrt(cf*1.0**2/two) ! uTau=sqrt(Cf*Uinf^2/2) and Uinf=1
                    retau = (0.01665*xcenter+0.16158)*utau/1.5d-5 ! ReTau=delta_995*uTau/nu
                    tmp = 3.9*sqrt(retau) ! h+ = 3.9*sqrt(ReTau)
                    tmp = tmp*1.5e-5/utau
                    alphaI = 0.250d0-dwall(e)/tmp
                    fbarg = two*exp(-9.0d0*alphaI**2)
                    fb = min(fbarg,one)
                    fd = fb
                  endif
                  if (alphaI.ge.zero) then
                     fe1 = two*exp(-11.09d0*alphaI**2)
                  else
                     fe1 = two*exp(-9.0d0*alphaI**2)
                  endif
                  fe = fe2*max((fe1-one),zero)
                  lIDDES = fd*(one+fe)*lRANS+(one-fd)*lLES
                  IDDESfun(e,1) = fd
                  IDDESfun(e,2) = fe1
                  IDDESfun(e,3) = fe2
                  IDDESfun(e,4) = fe
                  IDDESfun(e,5) = lIDDES
               endif
           endif

c          Start calculating source terms in k and omega eqns
c          k production
           Pk1 = mut * ( alpha + theta - 2.0d0 * chi/3.0d0 ) ! chi=0 for incompressible
           Pk2 = 2.0/3.0 * Rou * k * psi ! psi=0 for incompressible
           Pk = max(Pk1 - Pk2, 0d0) 
           Pk_max = max(10.d0 * CmuKW * Rou * omg * k, 0d0) ! production limiter
           if (Pk.lt.Pk_max) then ! Pk smaller than limiter
             Pkrat = Pk * kInv
             Pkrat_k = zero
             Pk1_k = (alpha + theta - 2.0/3.0 * chi) * mut_k !dPk1/dk
             Pk2_k = 2.0/3.0 * Rou * psi !dPk2/dk
             Pk_k = Pk1_k-Pk2_k
           else ! P1 equals production limiter
!             write(*,*) 'Production limited to 10 times dissipation'
             Pk = Pk_max
             Pkrat = Pk * kInv
             Pkrat_k = zero
             Pk_k = max(10.d0 * CmuKW * Rou * omg, zero)
           endif
c          k destruction
           if (iLES.eq.0) then 
             Dk = max(CmuKW * Rou * omg * k, zero) !Dissipation term for RANS
             Dkrat = max(CmuKW * Rou * omg, zero)
             Dkrat_k = zero
             Dk_k = CmuKW * Rou * omg
           else if (iLES.eq.-2) then
             kP3 = k*k*k
             lDDESInv = max(lDDES,1.0d-10)
             Dk = Rou * sqrt(kP3) / lDDESInv !Dissipation term for DDES
           else if (iLES.eq.-3) then
             kP3 = k*k*k
             if (lIDDES.lt.1.0d-12) then
                !write(*,*) 'Small l_IDDES'
             endif
             lIDDESInv = max(lIDDES,1.0d-12)
             Dk = Rou * sqrt(kP3) / lIDDESInv !Dissipation term for IDDES
             Dkrat = Dk * kInv
             Dk_k = Rou * CmuKW * omg / (fd * (1+fe))
             Dkrat_k = zero
           endif
c          k full source and LHS
           Sk = Pk - Dk ! Sk=prod-diss full source for k
           tmp = zero
           if (Pkrat.lt.zero) tmp = Pkrat
           if (Pkrat_k.lt.zero) tmp = tmp+Pkrat_k*k
           if (Dkrat.gt.zero) tmp = tmp-Dkrat
           if (Dkrat_k.gt.zero) tmp = tmp-Dkrat_k*k
c           if ((Dkrat-Pkrat).lt.zero) tmp = Dkrat-Pkrat
c           if ((Dkrat_k-Pkrat_k).lt.zero) tmp = tmp-
!           Sk_k = -one*tmp
           Sk_k = CmuKW * Rou * omg ! Menter's suggestion in 1992 paper
           k_Prod(e) = Pk
           k_Dest(e) = Dk
c          omega production
           Pw = Pk*omicron*Rou*mutInv
           Pwrat = Pw * omgInv
           Pw_w = zero
           Pwrat_w = zero
c          omega destruction
           Dw = max(beta * Rou * ( omg ** 2  ), 0d0) ! destruction
           Dwrat = max(beta * Rou * omg, 0d0)
           Dw_w = two * beta * Rou * omg
           Dwrat_w = zero
c          omega cross diffusion
           CDw = two * (one-F1) * Rou * sigw2 * delKdelW * omgInv ! cross diffusion
           CDw_w = -CDw*omgInv
           add(e,1) = zero
           add(e,2) = zero
           add(e,3) = zero
c          add is then added to the velocity to make uMod
           !add(e,1) = (one-F1)*two*Rou*sigw2*gradK(e,1)*omgInv
           !add(e,2) = (one-F1)*two*Rou*sigw2*gradK(e,2)*omgInv
           !add(e,3) = (one-F1)*two*Rou*sigw2*gradK(e,3)*omgInv
           !CDw = zero
           !CDw_w = zero
c          omega full source and LHS
           Sw = Pw - Dw + CDw ! the full source term for w equation
           w_Prod(e) = Pw
           w_Dest(e) = Dw
           w_CD(e) = CDw ! cross diffusion contribution to source
           tmp = zero
           if (Pwrat.lt.zero) tmp = Pwrat
           if (Pwrat_w.lt.zero) tmp = tmp+Pwrat*omg
           if (Dwrat.gt.zero) tmp = tmp-Dwat
           if (Dwrat_w.gt.zero) tmp = tmp-Dwrat_w*omg
           if (CDw_w.lt.zero) tmp = tmp+CDw_w
!           Sw_w = -one*tmp
           Sw_w = (abs(CDw)+two*Dw)*omgInv ! Menter's suggestion in 1992 paper

c          Rate of change of source term (srcRat)
           k_Rat(e) = Pkrat
           w_Rat(e) = Pwrat

c          SST-sust (refer to NASA turbulence modelling website)
           k_amb = zero
           Omega_amb = zero
           Sk_sust = CmuKW * Rou * k_amb * Omega_amb
           Sw_sust = beta * Rou * Omega_amb ** 2
           Sk = Sk + Sk_sust 
           Sw = Sw + Sw_sust 

c          Add the k and w equation source terms to the src srray
           src(e,1) = Sk
           src(e,2) = Sw

c          Fill srcRat
           srcRat(e,1) = Sk * kInv
           srcRat(e,2) = Sw * omgInv

c          Fill srcJac for the LHS of the source terms
           srcJac(e,1) = Sk_k ! d(Sk)/d(k)
           srcJac(e,3) = zero ! d(Sk)/d(omega)
           srcJac(e,2) = zero ! d(Sw)/d(k)
           srcJac(e,4) = Sw_w ! d(Sw)/d(omega)
           if (srcJac(e,1).le.zero.or.srcJac(e,4).le.zero) then
           !   write(*,*) 'Jacobian is zero or negative'
           endif
!           srcJac(e,1) = max( srcJac(e,1), zero )
!           srcJac(e,4) = max( srcJac(e,4), zero )

        enddo ! End loop over elements in this block
      
c       Pick which scalar to pass out
        if(isclr.eq.1)then        ! kay
           srcRat1 = srcRat(:,1)
           src1 = src(:,1)
        else if (isclr.eq.2) then ! omega
           srcRat1 = srcRat(:,2)
           src1 = src(:,2)
        endif

      else ! Doing Advection-diffusion
         srcrat1 = zero
         src1 = zero
         srcjac = zero
      endif


      return
      end

      
      subroutine getblendfunc (delKdelW, kay, omg, dwl, rho, nu, F1, F2)
        use turbKW ! access to KW model constants
        include "common.h"
        real*8 kay, omg, nu, dwl, rho, delKdelW
        real*8 F1, F2
        real*8 kq, omgInv
        real*8 arg11, arg11den, arg12, arg12num, arg12den
        real*8 arg13, arg13num, arg13den, argCD_kw, CD_kw
        real*8 arg1, arg21, arg22, arg2


        omgInv = one / omg
        if (omg . lt . 1.0d-12) omgInv=1.0d12
        kq = sqrt(kay)
        arg11den = CmuKW*omg*dwl
        arg11    = kq / arg11den
        arg12num = 500.0d0*nu
        arg12den = dwl * dwl * omg
        arg12    = arg12num / arg12den
        arg13num = 4.0d0*rho*sigw2*kay
        argCD_kw = 2.0d0*(rho*sigw2*omgInv)*delKdelW
        CD_kw = max(argCD_kw,1.0d-10)
        arg13den = CD_kw*dwl*dwl
        arg13  = arg13num/arg13den
        arg1 = min(max(arg11 , arg12) , arg13 )
        F1 = tanh(arg1 ** 4)
        arg2 = max(2.0d0*arg11 , arg12)
        F2 = tanh(arg2 * arg2)


      return
      end


!calculates the value of gradK and gradW for use in delkdelw. Cannot use 
!e3qvarKWsclr as we have access to shg, not shgl in the elm3komega
!subroutine.


      subroutine gradKWgen (blk, yl, shg, gradK, gradW )
 
      use eblock
      include "common.h"
      type (LocalBlkData) blk
      real*8   yl(bsz,blk%s,ndof),    shp(blk%e,blk%s),
     &         shg(blk%e,blk%s,nsd),  
     &         gradK(blk%e,nsd),    gradW(blk%e,nsd)
      integer n
 
      gradK = zero
      gradW = zero
      do n = 1, blk%s
           gradK(:,1) = gradK(:,1) + shg(:,n,1) * yl(1:blk%e,n,6)
           gradK(:,2) = gradK(:,2) + shg(:,n,2) * yl(1:blk%e,n,6)
           gradK(:,3) = gradK(:,3) + shg(:,n,3) * yl(1:blk%e,n,6)
      enddo
      do n = 1, blk%s
           gradW(:,1) = gradW(:,1) + shg(:,n,1) * yl(1:blk%e,n,7)
           gradW(:,2) = gradW(:,2) + shg(:,n,2) * yl(1:blk%e,n,7)
           gradW(:,3) = gradW(:,3) + shg(:,n,3) * yl(1:blk%e,n,7)
      enddo
 
      return
      end

