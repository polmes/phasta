      subroutine asithf (blk, y, x, strnrm, ien, fres, shgl, shp, Qwtf)

      use eblock

      include "common.h"
      type (LocalBlkData) blk

      dimension y(nshg,ndofl),            fres(nshg,24)
      dimension x(numnp,nsd),            xl(bsz,blk%n,nsd)
      dimension ien(blk%e,blk%s),        ycl(bsz,blk%s,ndofl),
     &          fresl(bsz,24),        WdetJ(blk%e),
     &          u1(blk%e),              u2(blk%e),
     &          u3(blk%e),              dxdxi(blk%e,nsd,nsd),
     &          strnrm(blk%e,blk%g),    dxidx(blk%e,nsd,nsd),
     &          shgl(nsd,blk%s,blk%g),       shg(blk%e,blk%s,nsd),
     &          shp(blk%s,blk%g),
     &          fresli(blk%e,24),       Qwtf(ngaussf)

      dimension tmp(blk%e)

      call localy (blk, y,      ycl,     ien,    ndofl,  'gather  ')
      call localx (blk, x,      xl,     ien,    nsd,  'gather  ')
c

      if(matflg(1,1).eq.0) then ! compressible
         ycl (:,:,1) = ycl(:,:,1) / (Rgas * ycl(:,:,5)) !get density
      else
         ycl(:,:,1) = one ! Even if density non unity, it would cancel out
      endif

      fresl = zero

      do intp = 1, ngaussf


c  calculate the metrics
c
c
c.... --------------------->  Element Metrics  <-----------------------
c
c.... compute the deformation gradient
c
        dxdxi = zero
c
          do n = 1, blk%n
            dxdxi(1:blk%e,1,1) = dxdxi(1:blk%e,1,1) + xl(1:blk%e,n,1) * shgl(1,n,intp)
            dxdxi(1:blk%e,1,2) = dxdxi(1:blk%e,1,2) + xl(1:blk%e,n,1) * shgl(2,n,intp)
            dxdxi(1:blk%e,1,3) = dxdxi(1:blk%e,1,3) + xl(1:blk%e,n,1) * shgl(3,n,intp)
            dxdxi(1:blk%e,2,1) = dxdxi(1:blk%e,2,1) + xl(1:blk%e,n,2) * shgl(1,n,intp)
            dxdxi(1:blk%e,2,2) = dxdxi(1:blk%e,2,2) + xl(1:blk%e,n,2) * shgl(2,n,intp)
            dxdxi(1:blk%e,2,3) = dxdxi(1:blk%e,2,3) + xl(1:blk%e,n,2) * shgl(3,n,intp)
            dxdxi(1:blk%e,3,1) = dxdxi(1:blk%e,3,1) + xl(1:blk%e,n,3) * shgl(1,n,intp)
            dxdxi(1:blk%e,3,2) = dxdxi(1:blk%e,3,2) + xl(1:blk%e,n,3) * shgl(2,n,intp)
            dxdxi(1:blk%e,3,3) = dxdxi(1:blk%e,3,3) + xl(1:blk%e,n,3) * shgl(3,n,intp)
          enddo
c
c.... compute the inverse of deformation gradient
c
        dxidx(1:blk%e,1,1) =   dxdxi(1:blk%e,2,2) * dxdxi(1:blk%e,3,3)
     &                 - dxdxi(1:blk%e,3,2) * dxdxi(1:blk%e,2,3)
        dxidx(1:blk%e,1,2) =   dxdxi(1:blk%e,3,2) * dxdxi(1:blk%e,1,3)
     &                 - dxdxi(1:blk%e,1,2) * dxdxi(1:blk%e,3,3)
        dxidx(1:blk%e,1,3) =   dxdxi(1:blk%e,1,2) * dxdxi(1:blk%e,2,3)
     &                 - dxdxi(1:blk%e,1,3) * dxdxi(1:blk%e,2,2)
        tmp          = one / ( dxidx(1:blk%e,1,1) * dxdxi(1:blk%e,1,1)
     &                       + dxidx(1:blk%e,1,2) * dxdxi(1:blk%e,2,1)
     &                       + dxidx(1:blk%e,1,3) * dxdxi(1:blk%e,3,1) )
        dxidx(1:blk%e,1,1) = dxidx(1:blk%e,1,1) * tmp
        dxidx(1:blk%e,1,2) = dxidx(1:blk%e,1,2) * tmp
        dxidx(1:blk%e,1,3) = dxidx(1:blk%e,1,3) * tmp
        dxidx(1:blk%e,2,1) = (dxdxi(1:blk%e,2,3) * dxdxi(1:blk%e,3,1)
     &                - dxdxi(1:blk%e,2,1) * dxdxi(1:blk%e,3,3)) * tmp
        dxidx(1:blk%e,2,2) = (dxdxi(1:blk%e,1,1) * dxdxi(1:blk%e,3,3)
     &                - dxdxi(1:blk%e,3,1) * dxdxi(1:blk%e,1,3)) * tmp
        dxidx(1:blk%e,2,3) = (dxdxi(1:blk%e,2,1) * dxdxi(1:blk%e,1,3)
     &                - dxdxi(1:blk%e,1,1) * dxdxi(1:blk%e,2,3)) * tmp
        dxidx(1:blk%e,3,1) = (dxdxi(1:blk%e,2,1) * dxdxi(1:blk%e,3,2)
     &                - dxdxi(1:blk%e,2,2) * dxdxi(1:blk%e,3,1)) * tmp
        dxidx(1:blk%e,3,2) = (dxdxi(1:blk%e,3,1) * dxdxi(1:blk%e,1,2)
     &                - dxdxi(1:blk%e,1,1) * dxdxi(1:blk%e,3,2)) * tmp
        dxidx(1:blk%e,3,3) = (dxdxi(1:blk%e,1,1) * dxdxi(1:blk%e,2,2)
     &                - dxdxi(1:blk%e,1,2) * dxdxi(1:blk%e,2,1)) * tmp
c
c        wght=Qwt(lcsyst,intp)  ! may be different now
        wght=Qwtf(intp)
        WdetJ = wght / tmp
c
      fresli=zero
c
      if(matflg(1,1).eq.0) then ! compressible
         do i=1,blk%s
            fresli(1:blk%e,22) = fresli(1:blk%e,22)+shp(i,intp)*ycl(1:blk%e,i,1) !density at qpt
         enddo
      else   ! incompressible, set density
         fresli(1:blk%e,22)= one ! reduce comp2incompr regardless of rho  datmat(1,1,1)
      endif
c
      do n = 1,blk%s
        shg(1:blk%e,n,1) = (shgl(1,n,intp) * dxidx(1:blk%e,1,1)
     &              + shgl(2,n,intp) * dxidx(1:blk%e,2,1)
     &              + shgl(3,n,intp) * dxidx(1:blk%e,3,1))
        shg(1:blk%e,n,2) = (shgl(1,n,intp) * dxidx(1:blk%e,1,2)
     &              + shgl(2,n,intp) * dxidx(1:blk%e,2,2)
     &              + shgl(3,n,intp) * dxidx(1:blk%e,3,2))
        shg(1:blk%e,n,3) = (shgl(1,n,intp) * dxidx(1:blk%e,1,3)
     &              + shgl(2,n,intp) * dxidx(1:blk%e,2,3)
     &              + shgl(3,n,intp) * dxidx(1:blk%e,3,3))
      enddo

      do j=10,12  ! normal strainrate u_{i,i} no sum on i
       ig=j-9
       iv=j-8
       do i=1,blk%s
        fresli(1:blk%e,j) = fresli(1:blk%e,j)+shg(1:blk%e,i,ig)*ycl(1:blk%e,i,iv)
       enddo
      enddo

c shear stresses  NOTE  there may be faster ways to do this
c                  check agains CM5 code for speed WTP
       
       do i=1,blk%s
        fresli(1:blk%e,13) = fresli(1:blk%e,13)+shg(1:blk%e,i,2)*ycl(1:blk%e,i,2)
     &                             +shg(1:blk%e,i,1)*ycl(1:blk%e,i,3)
        fresli(1:blk%e,14) = fresli(1:blk%e,14)+shg(1:blk%e,i,3)*ycl(1:blk%e,i,2)
     &                             +shg(1:blk%e,i,1)*ycl(1:blk%e,i,4)
        fresli(1:blk%e,15) = fresli(1:blk%e,15)+shg(1:blk%e,i,3)*ycl(1:blk%e,i,3)
     &                             +shg(1:blk%e,i,2)*ycl(1:blk%e,i,4)
       enddo

      fresli(1:blk%e,13) = pt5 * fresli(1:blk%e,13)
      fresli(1:blk%e,14) = pt5 * fresli(1:blk%e,14)
      fresli(1:blk%e,15) = pt5 * fresli(1:blk%e,15)

      strnrm(1:blk%e,intp) = fresli(1:blk%e,22) * sqrt(
     &   two * (fresli(1:blk%e,10)**2 + fresli(1:blk%e,11)**2 + fresli(1:blk%e,12)**2)
     &  + four * ( fresli(1:blk%e,13)**2 + fresli(1:blk%e,14)**2 + 
     &    fresli(1:blk%e,15)**2 ) )

c
c S_ij
c

      fresli(1:blk%e,10) = fresli(1:blk%e,10) * WdetJ ! u_{1,1}*WdetJ
      fresli(1:blk%e,11) = fresli(1:blk%e,11) * WdetJ ! u_{2,2}*WdetJ
      fresli(1:blk%e,12) = fresli(1:blk%e,12) * WdetJ ! u_{3,3}*WdetJ
      fresli(1:blk%e,13) = fresli(1:blk%e,13) * WdetJ ! (1/2)*(u_{1,2}+u_{2,1})*WdetJ
      fresli(1:blk%e,14) = fresli(1:blk%e,14) * WdetJ ! (1/2)*(u_{1,3}+u_{3,1})*WdetJ
      fresli(1:blk%e,15) = fresli(1:blk%e,15) * WdetJ ! (1/2)*(u_{2,3}+u_{3,2})*WdetJ

      fresli(1:blk%e,22) = fresli(1:blk%e,22) * WdetJ   !rho * WdetJ
c     fresli(:,24) = fresli(:,24) * WdetJ
     
      u1=zero
      u2=zero
      u3=zero
      do i=1,blk%s
       u1 = u1 + shp(i,intp)*ycl(1:blk%e,i,2)
       u2 = u2 + shp(i,intp)*ycl(1:blk%e,i,3)
       u3 = u3 + shp(i,intp)*ycl(1:blk%e,i,4)
      enddo

      fresli(1:blk%e,1) = fresli(1:blk%e,22) * u1   !rho u1 * WdetJ
      fresli(1:blk%e,2) = fresli(1:blk%e,22) * u2   !rho u2 * WdetJ
      fresli(1:blk%e,3) = fresli(1:blk%e,22) * u3   !rho u3 * WdetJ

      fresli(1:blk%e,4) = fresli(1:blk%e,1) * u1    !rho u1 u1 *WdetJ
      fresli(1:blk%e,5) = fresli(1:blk%e,2) * u2    !rho u2 u2 *WdetJ
      fresli(1:blk%e,6) = fresli(1:blk%e,3) * u3    !rho u3 u3 *WdetJ
      fresli(1:blk%e,7) = fresli(1:blk%e,1) * u2    !rho u1 u2 *WdetJ
      fresli(1:blk%e,8) = fresli(1:blk%e,1) * u3    !rho u1 u3 *WdetJ
      fresli(1:blk%e,9) = fresli(1:blk%e,2) * u3    !rho u2 u3 *WdetJ

      fresli(1:blk%e,16) = strnrm(1:blk%e,intp) * fresli(1:blk%e,10) ! rho *|Eps| *Eps11 *WdetJ
      fresli(1:blk%e,17) = strnrm(1:blk%e,intp) * fresli(1:blk%e,11) ! rho *|Eps| *Eps22 *WdetJ
      fresli(1:blk%e,18) = strnrm(1:blk%e,intp) * fresli(1:blk%e,12) ! rho *|Eps| *Eps33 *WdetJ
      fresli(1:blk%e,19) = strnrm(1:blk%e,intp) * fresli(1:blk%e,13) ! rho *|Eps| *Eps12 *WdetJ
      fresli(1:blk%e,20) = strnrm(1:blk%e,intp) * fresli(1:blk%e,14) ! rho *|Eps| *Eps13 *WdetJ
      fresli(1:blk%e,21) = strnrm(1:blk%e,intp) * fresli(1:blk%e,15) ! rho *|Eps| *Eps23 *WdetJ

      fresli(1:blk%e,23) = WdetJ   !    Integral of 1 over the element
c
      do i = 1, 23
         fresl(1:blk%e,i) = fresl(1:blk%e,i) + fresli(1:blk%e,i)
      enddo
   
      enddo !end of loop over integration points
c
      do j = 1,blk%s
      do nel = 1,blk%e
        fres(ien(nel,j),:) = fres(ien(nel,j),:) + fresl(nel,:) 
      enddo
      enddo

      return
      end









