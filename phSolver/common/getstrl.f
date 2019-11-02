      subroutine getstrl(blk, y, xl, ien, strnrm, shgl, shp )

      use eblock

      include "common.h"
      type (LocalBlkData) blk

      dimension y(nshg,ndofl)
      dimension x(numnp,nsd),            xl(blk%e,blk%n,nsd)
      dimension ien(blk%e,blk%s),        yl(blk%e,blk%s,ndofl),
     &          u1(blk%e),              u2(blk%e),
     &          u3(blk%e),              dxdxi(blk%e,nsd,nsd),
     &          strnrm(blk%e,blk%g),    dxidx(blk%e,nsd,nsd),
     &          shgl(nsd,blk%s,blk%g),       shg(blk%e,blk%s,nsd),
     &          shp(blk%s,blk%g)         
      dimension tmp(blk%e),             fresli(blk%e,24)

      call localy (blk, y,      yl,     ien,    ndofl,  'gather  ')
!      call localx (blk, x,      xl,     ien,    nsd,  'gather  ')
c

      if(matflg(1,1).eq.0) then ! compressible
         yl (:,:,1) = yl(:,:,1) / (Rgas * yl(:,:,5)) 
      else
         yl (:,:,1) = one
      endif

      do intp = 1, ngauss

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

      fresli=zero
      do i=1,blk%s
        fresli(1:blk%e,22) = fresli(1:blk%e,22)+shp(i,intp)*yl(1:blk%e,i,1)  ! density at qpt
c       fresli(:,24) = fresli(:,24)+shp(i,intp)*yl(:,i,5)  !temperature at qpt
      enddo
c
c
c     fresli(:,22)=fresli(:,22)*wght
c     fresli(:,24)=fresli(:,24)*wght


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
        fresli(1:blk%e,j) = fresli(1:blk%e,j)+shg(1:blk%e,i,ig)*yl(1:blk%e,i,iv)
       enddo
      enddo

c shear stresses  NOTE  there may be faster ways to do this
c                  check agains CM5 code for speed WTP
       
       do i=1,blk%s
        fresli(1:blk%e,13) = fresli(1:blk%e,13)+shg(1:blk%e,i,2)*yl(1:blk%e,i,2)
     &                             +shg(1:blk%e,i,1)*yl(1:blk%e,i,3)
        fresli(1:blk%e,14) = fresli(1:blk%e,14)+shg(1:blk%e,i,3)*yl(1:blk%e,i,2)
     &                             +shg(1:blk%e,i,1)*yl(1:blk%e,i,4)
        fresli(1:blk%e,15) = fresli(1:blk%e,15)+shg(1:blk%e,i,3)*yl(1:blk%e,i,3)
     &                             +shg(1:blk%e,i,2)*yl(1:blk%e,i,4)
       enddo


      fresli(1:blk%e,13) = pt5 * fresli(1:blk%e,13)
      fresli(1:blk%e,14) = pt5 * fresli(1:blk%e,14)
      fresli(1:blk%e,15) = pt5 * fresli(1:blk%e,15)

      strnrm(1:blk%e,intp) = fresli(1:blk%e,22) * sqrt(
     &   two * (fresli(1:blk%e,10)**2 + fresli(1:blk%e,11)**2 + fresli(1:blk%e,12)**2)
     &  + four * ( fresli(1:blk%e,13)**2 + fresli(1:blk%e,14)**2 + 
     &    fresli(1:blk%e,15)**2 ) )

      
      enddo !end of loop over integration points

      
      return
      end
