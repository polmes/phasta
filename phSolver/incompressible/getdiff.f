      subroutine getDiff(blk, ith, dwl,yl, shape, xmudmi, xl, rmu,  rho,
     &                    elem_size)
c-----------------------------------------------------------------------
c  compute and add the contribution of the turbulent
c  eddy viscosity to the molecular viscosity.
c-----------------------------------------------------------------------
      use     turbSA
      use  spat_var_eps !* for spatially varying epsilon_ls  
      use eblock
      include "common.h"
      type (LocalBlkData) blk


      real*8  yl(bsz,blk%s,ndof), rmu(blk%e), xmudmi(blk%e,blk%g),
     &        shape(blk%e,blk%s),   rho(blk%e),
     &        dwl(bsz,blk%n),     sclr(blk%e),
     &        xl(bsz,blk%n,nsd)
      real*8 elem_size(blk%e)
      integer n, e

      real*8  epsilon_ls, kay, epsilon
     &        h_param, prop_blend(blk%e),test_it(blk%e)
c
c    
c.... get the material properties (2 fluid models will need to determine
c     the "interpolated in phase space" properties....constant for now.
c     two options exist in the interpolation 1) smooth (recommended) 
c     interpolation of nodal data, 2) discontinuous "sampling" phase at 
c     quadrature points.
c
!
!    prop_blend is a smoothing function to avoid possible large density 
!   gradients, e.g., water and air flows where density ratios can approach
!   1000.
!
!    epsilon_ls is an adjustment to control the width of the band over which
!    the properties are blended. 



      if (iLSet .eq. 0)then

         rho  = datmat(1,1,1)	! single fluid model, i.e., only 1 density
         rmu = datmat(1,2,1)

      else     !  two fluid properties used in this model

!        Smooth the tranistion of properties for a "distance" of epsilon_ls
!        around the interface.  Here "distance" is define as the value of the 
!        levelset function.  If the levelset function is properly defined, 
!        this is the true distance normal from the front.  Of course, the 
!        distance is in a driection normal to the front.

         Sclr = zero
         isc=abs(iRANS)+6
         do n = 1, blk%s
            Sclr = Sclr + shape(1:blk%e,n) * yl(1:blk%e,n,isc)
         enddo
         if (icode .ge. 20) then !scalar 2
         do i= 1, blk%e
	     epsilon_lsd_tmp = epsilon_lsd*elem_size(i)
c
             if (sclr(i) .lt. - epsilon_lsd_tmp) then
                 prop_blend(i) = zero
              elseif  (abs(sclr(i)) .le. epsilon_lsd_tmp) then
                 prop_blend(i) = 0.5*(one + Sclr(i) / epsilon_lsd_tmp +
     &                           (sin(pi*Sclr(i)/epsilon_lsd_tmp))/pi)
              elseif (sclr(i) .gt. epsilon_lsd_tmp) then
                 prop_blend(i) = one
              endif
           enddo
        else !scalar 1 or flow
            do i= 1, blk%e
              epsilon_ls_tmp = epsilon_ls*elem_size(i)
              if (sclr(i) .lt. - epsilon_ls_tmp) then
                 prop_blend(i) = zero	
              elseif  (abs(sclr(i)) .le. epsilon_ls_tmp) then
                 prop_blend(i) = 0.5*(one + Sclr(i)/epsilon_ls_tmp +
     &                           (sin(pi*Sclr(i)/epsilon_ls_tmp))/pi)
              elseif (sclr(i) .gt. epsilon_ls_tmp) then
                 prop_blend(i) = one
              endif
           enddo
        endif
c
        rho = datmat(1,1,2) + (datmat(1,1,1)-datmat(1,1,2))*prop_blend
        rmu = datmat(1,2,2) + (datmat(1,2,1)-datmat(1,2,2))*prop_blend

      endif

!	At this point we have a rho that is bounded by the two values for
! 	density 1, datmat(1,1,1), the fluid,  and density 2, datmat(1,1,2)
!     the gas

c
c  The above approach evaluates all intermediate quantities at the 
c  quadrature point, then combines them to form the needed quantities there.
c  1 alternative is calculating all quanties (only rho here) at the nodes and 
c  then interpolating the result to the quadrature points.  If this is done,
c  do not forget to do the same thing for rou in e3b!!!
c  ^^^^^^^^^^
c  ||||||||||
c  WARNING
c
c.... dynamic model
c      
      if (iLES .gt. 0 .and. iRANS.eq.0) then   ! simple LES
         rmu = rmu + xmudmi(:,ith)
      else if (iRANS.lt.0) then 
         if (iRANS .eq. -1) then ! RANS (Spalart-Allmaras)
            call AddEddyViscSA(blk,yl, shape, rmu)
         else if(iRANS.eq.-2) then ! RANS-KE
            sigmaInv=1.0        ! use full eddy viscosity for flow equations
            call AddEddyViscKE(blk,yl, dwl, shape, rho, sigmaInv, rmu)
         endif
         if (iLES.gt.0) then    ! this is DES so we have to blend in
                                ! xmudmi based on max edge length of
                                ! element
            call EviscDESIC (blk,xl,rmu,xmudmi)
         endif
      endif                     ! check for LES or RANS
c
      return
      end

      subroutine EviscDESIC(blk,xl,xmut,xmudmi)
     
      use eblock
      include "common.h"
      type (LocalBlkData) blk
      
      real*8 xmut(blk%e),xl(bsz,blk%n,nsd),xmudmi(blk%e,blk%g)


      do i=1,blk%e
         dx=maxval(xl(i,:,1))-minval(xl(i,:,1))
         dy=maxval(xl(i,:,2))-minval(xl(i,:,2))
         dz=maxval(xl(i,:,3))-minval(xl(i,:,3))
         emax=max(dx,max(dy,dz))
         if(emax.lt.eles) then  ! pure les
            xmut(i)=xmudmi(i,ith)
         else if(emax.lt.two*eles) then ! blend
            xi=(emax-eles)/(eles)
            xmut(i)=xi*xmut(i)+(one-xi)*(xmudmi(1,ith)+datmat(1,2,2))
         endif                  ! leave at RANS value as edge is twice pure les
      enddo
 !this was made messy by the addEddyVisc routines  Must clean up later.


      
      return
      end
       
      subroutine getdiffsclr(blk,shape, dwl, yl, diffus)

      use turbSA
      use turbKE ! access to KE model constants
      use eblock
      include "common.h"
      type (LocalBlkData) blk
      
      real*8   diffus(blk%e), rho(blk%e)
      real*8   yl(bsz,blk%s,ndof), dwl(bsz,blk%n), shape(blk%e,blk%s)
      integer n, e
      rho(:)  = datmat(1,1,1)	! single fluid model, i.e., only 1 density
      if(isclr.eq.0) then  ! solving the temperature equation
         diffus(:) = datmat(1,4,1)
      else if(iRANS.eq.-1) then ! solving SA model
         diffus(:) = datmat(1,2,1)
         call AddSAVar(blk,yl, shape, diffus)
      else if(iRANS.eq.-2)then ! solving KE model
         diffus(:) = datmat(1,2,1)
         if(isclr.eq.2) then
            sigmaInv=1.0/ke_sigma ! different eddy viscosity for epsilon
         else
            sigmaInv=1.0 ! full eddy viscosity for solving kay equation
         endif
         call AddEddyViscKE(blk,yl, dwl, shape, rho, sigmaInv, diffus)
      else                      ! solving scalar advection diffusion equations
         diffus = scdiff(isclr)
      endif
c
      return
      end

      function ev2sa(xmut,rm,cv1)
      implicit none
      real*8 err,ev2sa,rm,cv1,f,dfds,rat,efac
      real*8 pt5,kappa,B,xmut,chi3,denom,cv1_3
      integer iter
      pt5=0.5
      err=1.0d-6
      ev2sa=rm*cv1*1.2599       ! inflection point chi=cv1*cuberoot(2)
      kappa=0.4
c$$$        B=5.5
      efac=0.1108               ! exp(-kappa*B)
      do iter=1,50
         chi3=ev2sa/rm
         chi3=chi3*chi3*chi3
         cv1_3=cv1**3
         denom=chi3+cv1_3
         
         f=ev2sa*chi3/denom - xmut
         dfds=chi3*(chi3+4.0*cv1_3)/(denom**2)
         rat=-f/dfds
         ev2sa=ev2sa+rat
         if(abs(rat).le.err) goto 20
      enddo
      write(*,*)'ev2sa failed to converge'
      write(*,*) 'dfds,        rat,        ev2sa,        mu'
      write(*,*) dfds,rat,ev2sa,rm
 20   continue
      return
      end  
c     
      

      subroutine AddEddyViscSA(blk, yl,shape,rmu)
      use turbSA
      use eblock
      include "common.h"
      type (LocalBlkData) blk
c     INPUTS
      double precision, intent(in), dimension(bsz,blk%s,ndof) ::
     &     yl
      double precision, intent(in), dimension(blk%e,blk%s) ::
     &     shape
c     INPUT-OUTPUTS
      double precision, intent(inout), dimension(blk%e) ::
     &     rmu
c     LOCALS
      logical, dimension(blk%s) ::
     &     wallnode
      integer e, n
      double precision xki, xki3, fv1, evisc
  
      if(itwmod.eq.-2) then  ! effective viscosity
        do e = 1, blk%e
c         assume no wall nodes on this element
          wallnode(:) = .false.
c         mark the wall nodes for this element, if there are any
          do n = 1, blk%s
            u1=yl(e,n,2)
            u2=yl(e,n,3)
            u3=yl(e,n,4)
            if((u1.eq.zero).and.(u2.eq.zero).and.(u3.eq.zero)) then
              wallnode(n)=.true.
            endif
          enddo
c
          if( any(wallnode(:)) ) then
c if there are wall nodes for this elt, then we are using effective-
c viscosity near-wall modeling, and eddy viscosity has been stored
c at the wall nodes in place of the spalart-allmaras variable; the
c eddy viscosity for the whole element is taken to be the avg of the
c wall values
            evisc = zero
            nwnode=0
            do n = 1, blk%s
              if(wallnode(n)) then
                evisc = evisc + yl(e,n,6)
                nwnode = nwnode + 1
              endif
            enddo
            evisc = evisc/nwnode
            rmu(e) = rmu(e) + abs(evisc)
c this is what we would use instead of the above if we were allowing
c the eddy viscosity to vary through the element based on non-wall nodes
c$$$               evisc = zero
c$$$               Turb = zero
c$$$               do n = 1, blk%s
c$$$                  if(wallmask(n).eq.1) then
c$$$                     evisc = evisc + shape(e,n) * yl(e,n,6)
c$$$                  else
c$$$                     Turb = Turb + shape(e,n) * yl(e,n,6)
c$$$                  endif
c$$$               enddo
c$$$               xki    = abs(Turb)/rmu(e)
c$$$               xki3   = xki * xki * xki
c$$$               fv1    = xki3 / (xki3 + saCv1P3)
c$$$               rmu(e) = rmu(e) + fv1*abs(Turb)               
c$$$               rmu(e) = rmu(e) + abs(evisc)
          else !no wall nodes
            Turb = zero
            do n = 1, blk%s
              Turb = Turb + shape(e,n) * yl(e,n,6)
            enddo
            xki    = abs(Turb)/rmu(e)
            xki3   = xki * xki * xki
            fv1    = xki3 / (xki3 + saCv1P3)
            rmu(e) = rmu(e) + fv1*abs(Turb)
          endif       
        enddo !loop over elements  
      else
c else one of the following is the case:
c   using slip-velocity
c   using no wall model; walls are resolved 
c in all of these cases, eddy viscosity is calculated normally
c
c     Loop over elements in this block
        if(blk%s.eq.4) then !hard code tets
          do e = 1, blk%e
            Turb=shape(e,1)*yl(e,1,6)
     &          +shape(e,2)*yl(e,2,6)
     &          +shape(e,3)*yl(e,3,6)
     &          +shape(e,4)*yl(e,4,6)
            Turb   = abs(Turb)
            xki    = Turb/rmu(e)
            xki3   = xki * xki * xki
            fv1    = xki3 / (xki3 + saCv1P3)
            rmu(e) = rmu(e) + fv1*Turb
          enddo
            
        else
          do e = 1, blk%e
            Turb = zero
            do n = 1, blk%s
               Turb = Turb + shape(e,n) * yl(e,n,6)
            enddo
            Turb   = abs(Turb)
            xki    = Turb/rmu(e)
            xki3   = xki * xki * xki
            fv1    = xki3 / (xki3 + saCv1P3)
            rmu(e) = rmu(e) + fv1*Turb
          enddo                     ! end loop over elts
        endif
      endif
      return
      end subroutine AddEddyViscSA



      subroutine AddSAVar(blk,yl,shape,rmu)
      use turbSA
      use eblock
      include "common.h"
      type (LocalBlkData) blk
c     INPUTS
      double precision, intent(in), dimension(bsz,blk%s,ndof) ::
     &     yl
      double precision, intent(in), dimension(blk%e,blk%s) ::
     &     shape
c     INPUT-OUTPUTS
      double precision, intent(inout), dimension(blk%e) ::
     &     rmu
c     LOCALS
      logical, dimension(blk%s) ::
     &     wallnode
      integer e, n
      double precision savar, savarw
c     Loop over elements in this block
      do e = 1, blk%e
c        assume no wall nodes on this element
         wallnode(:) = .false.
         if(itwmod.eq.-2) then  ! effective viscosity
c           mark the wall nodes for this element, if there are any
            do n = 1, blk%s
               u1=yl(e,n,2)
               u2=yl(e,n,3)
               u3=yl(e,n,4)
               if((u1.eq.zero).and.(u2.eq.zero).and.(u3.eq.zero))
     &              then
                  wallnode(n)=.true.
               endif
            enddo
         endif
c
         savar=zero
         do n = 1, blk%s
            if( wallnode(n) ) then
c if wallmask was set, we're using effective-viscosity wall-model and
c this node is on a wall.  Eddy viscosity has been stored at the wall 
c nodes in place of the S-A variable, so we must convert it
               savarw = ev2sa(yl(e,n,6),datmat(1,2,1),saCv1)
               savar  = savar + shape(e,n) * savarw
            else
c if wallmask wasn't set, then one of the following is the case:
c   using effective-viscosity, but this isn't a wall node
c   using slip-velocity
c   using no wall model; wall is resolved
c in all these cases, the S-A variable is calculated normally
               savar  = savar + shape(e,n) * yl(e,n,6)
            endif   
         enddo
         rmu(e)=datmat(1,2,1)
         rmu(e) = (rmu(e) + abs(savar)) * saSigmaInv
      enddo                     ! end loop over elts
      return
      end subroutine AddSAVar



      subroutine AddEddyViscKE(blk, yl, dwl, shape, rho, sigmaInv, rmu)
      use turbKE ! access to KE model constants
      use eblock
      include "common.h"
      type (LocalBlkData) blk
c     INPUTS
      double precision, intent(in), dimension(bsz,blk%s,ndof) ::
     &     yl
      double precision, intent(in), dimension(bsz,blk%n) ::
     &     shape, dwl
      double precision, intent(in), dimension(blk%e) ::
     &     rho
      double precision sigmaInv
c     INPUT-OUTPUTS
      double precision, intent(inout), dimension(blk%e) ::
     &     rmu
c     LOCALS
      double precision eviscKE, kay, epsilon, dw, CmuKE
      double precision epsInv, Rey, Ret, RetInv, tmp1, fmuKE
      integer e,n
c
      do e = 1, blk%e
         kay = 0.0
         epsilon = 0.0
         dw = 0.0
         do n = 1, blk%s
            kay = kay + shape(e,n)*yl(e,n,6)
            epsilon = epsilon + shape(e,n)*yl(e,n,7)
         enddo
         do n = 1, blk%n
            dw = dw + shape(e,n)*dwl(e,n)
         enddo
         kay = abs(kay)
         if(kay.lt.1.0e-32) kay=0.0
         epsInv	    = 0
         if ( abs(epsilon) .gt.1.e-32) then
            epsInv        = 1. / abs(epsilon)
         endif
         
         Rey                 = sqrt(kay) *    dw * rho(e) / rmu(e)       
         Ret                 = kay*kay   * epsInv * rho(e) / rmu(e)
         RetInv              = 0
         if(Ret.lt.1.d100.AND.Ret.gt.zero) RetInv=1./Ret
         tmp1     = exp(-0.0165*Rey)          
         fmuKE    = (1. -tmp1) ** 2 * (1.+20.5*RetInv) ! fmu of Lam-Bremhorst

         eviscKE=rho(e)*ke_C_mu*fmuKE*kay*kay*epsInv
         
         rmu(e) = rmu(e) + eviscKE*sigmaInv
      enddo
      return
      end subroutine AddEddyViscKE


