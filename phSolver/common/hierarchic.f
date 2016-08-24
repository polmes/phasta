c------------------------------------------------------------------------
c
c This file contains functions for dealing with higher order shape 
c functions at the element level.
c
c Christian Whiting, Winter 1999
c------------------------------------------------------------------------

      subroutine getsgn(ien, sgn)
c------------------------------------------------------------------------
c     returns the matrix of mode signs used for negating higher order
c     basis functions. Connectivity array is assumed to have negative
c     signs on all modes to be negated.
c------------------------------------------------------------------------
      include "common.h"

      dimension ien(ibksiz,blk%s),  sgn(ibksiz,blk%s)
      
      do i=nenl+1,blk%s
         where ( ien(1:ibksiz,i) < 0 )
            sgn(1:ibksiz,i) = -one
         elsewhere
            sgn(1:ibksiz,i) = one
         endwhere
      enddo
      
      return 
      end
      
      subroutine getshp(ith,shp, shgl, sgn, shape, shdrv)
c------------------------------------------------------------------------
c     returns the matrix of element shape functions with the higher
c     order modes correctly negated at the current quadrature point.
c------------------------------------------------------------------------
      include "common.h"
      
      dimension shp(blk%s,ngauss),   shgl(nsd,blk%s,ngauss),
     &          sgn(ibksiz,blk%s),     shape(ibksiz,blk%s),
     &          shdrv(ibksiz,nsd,blk%s)
      
      
      do i=1,nenl
         shape(1:ibksiz,i) = shp(i,ith)
         do j=1,3
            shdrv(1:ibksiz,j,i) = shgl(j,i,ith)
         enddo
      enddo
      if ( ipord > 1 ) then
         do i=nenl+1,blk%s
            shape(1:ibksiz,i) = sgn(1:ibksiz,i) * shp(i,ith)
            do j=1,3
               shdrv(1:ibksiz,j,i) = shgl(j,i,ith)*sgn(1:ibksiz,i) 
            enddo
         enddo
      endif
      
      return 
      end
      
      subroutine getshpb(shp, shgl, sgn, shape, shdrv)
c------------------------------------------------------------------------
c     returns the matrix of element shape functions with the higher
c     order modes correctly negated at the current quadrature point.
c------------------------------------------------------------------------
      include "common.h"

      dimension shp(blk%s,ngaussb),  shgl(nsd,blk%s,ngaussb),
     &          sgn(ibksiz,blk%s),     shape(ibksiz,blk%s),
     &          shdrv(ibksiz,nsd,blk%s)
      
      
      do i=1,nenl
         shape(1:ibksiz,i) = shp(i,intp)
         do j=1,3
            shdrv(1:ibksiz,j,i) = shgl(j,i,intp)
         enddo
      enddo
      if ( ipord > 1 ) then
         do i=nenl+1,blk%s
            shape(1:ibksiz,i) = sgn(1:ibksiz,i) * shp(i,intp)
            do j=1,3
               shdrv(1:ibksiz,j,i) = shgl(j,i,intp)*sgn(1:ibksiz,i) 
            enddo
         enddo
      endif
      
      return 
      end
      
      subroutine getbnodes(lnode)
c------------------------------------------------------------------------
c     compute the higher order modes that lie on the boundary of the 
c     element.
c------------------------------------------------------------------------
      include "common.h"
      
      dimension lnode(27)

c
c.... boundary triangle of tet element
c
      if (lcsyst .eq. 1) then
         do n = 1, nenbl
            lnode(n) = n
         enddo
         if( ipord>1 ) then
            nem = ipord-1
            do n=1,3*nem
               lnode(nenbl+n) = nenbl+1+n
            enddo
         endif
         if( ipord>2 ) then
            nfm = (ipord-2)*(ipord-1)/2
            do n=1,nfm
               lnode(3+3*nem+n) = 4+6*nem+n
            enddo
         endif
c
c.....boundary quadrilateral for a hex element
c
      else if(lcsyst .eq. 2) then
         do n = 1, nenbl
            lnode(n) = n
         enddo
         if( ipord > 1) then
            nem = ipord -1
            do n=1,4*nem
               lnode(nenbl+n) = 8+n
            enddo
         endif
         
         if( ipord > 3) then
            nfm = (ipord-2)*(ipord-3)/2
            do n=1,nfm
               lnode(4+4*nem+n) = 8+12*nem+n
            enddo
         endif         
c
c
c.... This renumbers the boundary nodes for a wedge element when the
c     boundary element is a quad face.  From lnode = [1 2 3 4]
c                                        to   lnode = [1 4 5 2]
c
      else if(lcsyst .eq. 3) then
         do n = 1, nenbl
            lnode(n) = n
         enddo
c
c     Need to implement for cubic, this valid only for ipord=2
c 
         if( ipord>1 ) then
            nem = ipord-1
            do n=1,3*nem
               lnode(nenbl+n) = 6+n
            enddo
         endif
c
c     Boundary quad of wedge element
c
      else if(lcsyst .eq. 4) then
         lnode(1) = 1
         lnode(2) = 4
         lnode(3) = 5
         lnode(4) = 2
c$$$c     
c$$$c     Need to implement for cubic, this valid only for ipord=2
c$$$c     
c$$$         if( ipord > 1) then
c$$$            lnode(5) = 9
c$$$            lnode(6) = 15
c$$$            lnode(7) = 12
c$$$            lnode(8) = 13
c$$$               nem = ipord -1
c$$$               do n=1,4*nem
c$$$                  lnode(nenbl+n) = 6+n
c$$$               enddo
c$$$         endif         
c     
c     Boundary quad of pyramid element
c
      else if(lcsyst .eq. 5) then
         lnode(1) = 1
         lnode(2) = 2
         lnode(3) = 3
         lnode(4) = 4
c$$$  c
c$$$  c     Need to implement for cubic, this valid only for ipord=2
c$$$  c          
c$$$            if( ipord > 1) then
c$$$               lnode(5) = 9
c$$$               lnode(6) = 15
c$$$               lnode(7) = 12
c$$$               lnode(8) = 13
c$$$               nem = ipord -1
c$$$               do n=1,4*nem
c$$$                  lnode(nenbl+n) = 6+n
c$$$               enddo
c$$$            endif         
c     
c     Boundary triangle of pyramid element
c
      else if(lcsyst .eq. 6) then
         lnode(1) = 1
         lnode(2) = 5
         lnode(3) = 2
c$$$c
c$$$c     Need to implement for cubic, this valid only for ipord=2
c$$$c          
c$$$            if( ipord > 1) then
c$$$               lnode(5) = 9
c$$$               lnode(6) = 15
c$$$               lnode(7) = 12
c$$$               lnode(8) = 13
c$$$               nem = ipord -1
c$$$               do n=1,4*nem
c$$$                  lnode(nenbl+n) = 6+n
c$$$               enddo
c$$$            endif         
c     
c.... other element types need to be implemented
c
      else
         write (*,*) 'Boundary element not implemented for lcyst='
     &        ,lcsyst
         stop
      endif
      
      return 
      end
      
c-----------------------------------------------------------------------
c
c  Evaluate coefficient vector at its interpolation points
c
c-----------------------------------------------------------------------
      subroutine evalAtInterp( ycoeff,  yvals,  x,   nvars, npts )

      use     pointer_data
      include "common.h"
      
      integer nvars, npts, nHits(nshg)
      
      real*8  ycoeff(nshg,ndof),   yvals(nshg,nvars),
     &        shp(blk%s,npts),      shgl(nsd,blk%s,npts),
     &        intpnt(3,npts),      x(numnp,nsd)
      
      real*8, allocatable :: ycl(1:ibksiz,:,:)
      real*8, allocatable :: xl(1:ibksiz,:,:)
      real*8, allocatable :: yvl(1:ibksiz,:,:)
      real*8, allocatable :: sgn(1:ibksiz,:)

      yvals = zero
c
c.... generate the shape functions at the interpolation points
c
      call getIntPnts(intpnt,npts)
      do i=1,npts
         call shpTet(ipord,intpnt(1:ibksiz,i),shp(1:ibksiz,i),shgl(1:ibksiz,:,i))
      enddo
c
c.... loop over element blocks
c
      nHits = 0
      do iblk = 1, nelblk
         iel    = lcblk(1,iblk)
         lcsyst = lcblk(3,iblk)
         nenl   = lcblk(5,iblk) ! no. of vertices per element
         blk%s   = lcblk(10,iblk)
         ndofl  = lcblk(8,iblk)
         ibksiz   = lcblk(1,iblk+1) - iel 

         allocate ( ycl(ibksiz,blk%s,ndof ) )
         allocate ( yvl(ibksiz,blk%s,nvars) )
         allocate ( xl(ibksiz,nenl,nsd   ) )
         allocate ( sgn(ibksiz,blk%s)       )
         
         call getsgn(mien(iblk)%p,sgn)
         
         call localy( ycoeff, ycl, mien(iblk)%p, ndof,  'gather  ')
         call localx( x,      xl,  mien(iblk)%p, nsd,   'gather  ')

         call eval  ( xl,       ycl,      yvl,      
     &                shp,      shgl,     sgn,      
     &                nvars,    npts    )

c
c.... average coefficients since stresses may be discontinuous
c         
         call localSum( yvals,    yvl,    mien(iblk)%p,  
     &                  nHits,    nVars)  
         
         
         deallocate ( ycl )
         deallocate ( yvl )
         deallocate ( sgn )
         deallocate ( xl  )
c
      enddo

c
c.... average the global values
c
      do i = 1, nshg
         do j = 1, nvars
            yvals(i,j) = yvals(i,j)/nHits(i) !(real(nHits(i),8))
         enddo
      enddo
      
      return
      end

c-----------------------------------------------------------------------
c
c  evaluate in element coordinate system
c
c-----------------------------------------------------------------------
      subroutine eval( xl,      ycl,     yvl,     
     &                 shp,     shgl,    sgn,
     &                 nvars,   npts ) 
      
      include "common.h"
      
      integer nvars
c
      real*8  ycl(ibksiz,blk%s,ndof),   yvl(ibksiz,blk%s,nvars),
     &        sgn(ibksiz,blk%s),        shape(ibksiz,blk%s),
     &        shdrv(ibksiz,nsd,blk%s),  shp(blk%s,npts),
     &        shgl(nsd,blk%s,npts),   xl(ibksiz,nenl,nsd),
     &        shg(ibksiz,blk%s,nsd),    gradV(ibksiz,nsd,nsd),
     &        dxidx(ibksiz,nsd,nsd),   tmp(ibksiz), wtmp
      
      yvl = zero
c
c.... loop over interpolation points
c
      do intp = 1, npts
         call getshp(intp, shp,          shgl,      sgn, 
     &               shape,        shdrv)
      
c
c.... pressure and velocity
c
         do i = 1, blk%s
            do j = 1, 4
               yvl(1:ibksiz,intp,j) = yvl(1:ibksiz,intp,j) + shape(1:ibksiz,i) * ycl(1:ibksiz,i,j)
            enddo
         enddo
c
c.... viscous stress
c
         call e3metric(intp, xl,         shdrv,      dxidx,  
     &                  shg,        tmp)

         gradV = zero
         do n = 1, blk%s
            gradV(1:ibksiz,1,1) = gradV(1:ibksiz,1,1) + shg(1:ibksiz,n,1) * ycl(1:ibksiz,n,2)
            gradV(1:ibksiz,2,1) = gradV(1:ibksiz,2,1) + shg(1:ibksiz,n,1) * ycl(1:ibksiz,n,3)
            gradV(1:ibksiz,3,1) = gradV(1:ibksiz,3,1) + shg(1:ibksiz,n,1) * ycl(1:ibksiz,n,4)
c     
            gradV(1:ibksiz,1,2) = gradV(1:ibksiz,1,2) + shg(1:ibksiz,n,2) * ycl(1:ibksiz,n,2)
            gradV(1:ibksiz,2,2) = gradV(1:ibksiz,2,2) + shg(1:ibksiz,n,2) * ycl(1:ibksiz,n,3)
            gradV(1:ibksiz,3,2) = gradV(1:ibksiz,3,2) + shg(1:ibksiz,n,2) * ycl(1:ibksiz,n,4)
c     
            gradV(1:ibksiz,1,3) = gradV(1:ibksiz,1,3) + shg(1:ibksiz,n,3) * ycl(1:ibksiz,n,2)
            gradV(1:ibksiz,2,3) = gradV(1:ibksiz,2,3) + shg(1:ibksiz,n,3) * ycl(1:ibksiz,n,3)
            gradV(1:ibksiz,3,3) = gradV(1:ibksiz,3,3) + shg(1:ibksiz,n,3) * ycl(1:ibksiz,n,4)
         enddo

         rmu = datmat(1,2,1)
            
         yvl(1:ibksiz,intp,6 ) = two * rmu * gradV(1:ibksiz,1,1)
         yvl(1:ibksiz,intp,7 ) = two * rmu * gradV(1:ibksiz,2,2)
         yvl(1:ibksiz,intp,8 ) = two * rmu * gradV(1:ibksiz,3,3)

         yvl(1:ibksiz,intp,9 ) = rmu * ( gradV(1:ibksiz,1,2) + gradV(1:ibksiz,2,1) )
         yvl(1:ibksiz,intp,10) = rmu * ( gradV(1:ibksiz,1,3) + gradV(1:ibksiz,3,1) )
         yvl(1:ibksiz,intp,11) = rmu * ( gradV(1:ibksiz,2,3) + gradV(1:ibksiz,3,2) )

c
c.... loop over interpolation points
c         
      enddo
      
      return
      end

         

