        subroutine getelsize (x,  shp,  shgl,  elem_size,
     &                           shpb, shglb,  elemb_size)
c
c----------------------------------------------------------------------
c
c 
c----------------------------------------------------------------------
c
      use pvsQbi  ! brings in NABI
      use stats   !  
      use pointer_data  ! brings in the pointers for the blocked arrays
      use local_mass
      use eblock
c
      include "common.h"
      type (LocalBlkData) blk
      include "mpif.h"


c
        dimension x(numnp,nsd)               
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT),
     &            shpb(MAXTOP,maxsh,MAXQPT),  
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT)  
c
	real*8 elem_size(numel), elemb_size(numelb)
        integer numel_array(numpe), numel_tot, index(numpe)
	real*8, allocatable :: elem_size_tot(:), elem_vol_tot(:)
	real*8 elem_vol(numel), elemb_vol(numelb)
        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)


c
c initialize array
c
        elem_size(:) = zero
        elem_vol(:) = zero
        numel_tot = 0
c
c loop over element blocks to compute volume
c
!              allocate (tmpshp(4,125))
!              allocate (tmpshgl(nsd,4,125))
        
           do iblk = 1, nelblk
              blk%i    = lcblk(1,iblk)
              blk%l = lcblk(3,iblk)
              blk%o = lcblk(4,iblk)
              blk%n   = lcblk(5,iblk) ! no. of vertices per element
              blk%s   = lcblk(10,iblk)
              blk%e   = lcblk(1,iblk+1) - blk%i 
              blk%g = nint(blk%l)
              allocate (tmpshp(blk%s,MAXQPT))
              allocate (tmpshgl(nsd,blk%s,MAXQPT))
              tmpshp(1:blk%s,:) = shp(blk%l,1:blk%s,:)
              tmpshgl(:,1:blk%s,:) = shgl(blk%l,:,1:blk%s,:)
c
c compute local element volume
c
               call e3elsize (blk,tmpshp, tmpshgl,
     &         mxl(iblk)%p,mien(iblk)%p,elem_size(blk%i),
     &         elem_vol(blk%i) )

               deallocate(tmpshp)
               deallocate(tmpshgl)
           enddo
c
c... Set Characteristic Size
c    i_spat_var_eps_flag = 0 - no spatially varying eps.  size = 1.0
c                          1 - based on vol**1/3
c                          2 - based on max edge length
c
        if (i_spat_var_eps_flag .eq. 0) then
          elem_size(:) = 1.0
        elseif (i_spat_var_eps_flag .eq. 1) then
          elem_size(:) = elem_vol(:)**(1.0/3.0)
        else
          elem_size(:) = elem_size(:)
        endif
c
c Now compute element volumes for boundary elements
c
c
c initialize array
c
        elemb_vol(:) = zero
        elemb_size(:) = zero	
c
c loop over element blocks to compute volume
c        
        do iblk = 1, nelblb
          blk%i    = lcblkb(1,iblk)
          blk%l = lcblkb(3,iblk)
          blk%o = lcblkb(4,iblk)
          blk%n   = lcblkb(5,iblk)  ! no. of vertices per element
          blk%s   = lcblkb(9,iblk)
          blk%e   = lcblkb(1,iblk+1) - blk%i 
          blk%g = nint(blk%l)
          allocate (tmpshp(blk%s,MAXQPT))
          allocate (tmpshgl(nsd,blk%s,MAXQPT))
          tmpshp(1:blk%s,:) = shp(blk%l,1:blk%s,:)
          tmpshgl(:,1:blk%s,:) = shgl(blk%l,:,1:blk%s,:)
c
c compute local element volume
c
           call e3elsize (blk,tmpshp,tmpshgl,
     &         mxlb(iblk)%p,mienb(iblk)%p,elemb_size(blk%i),
     &         elemb_vol(blk%i) )

          deallocate (tmpshp)
          deallocate (tmpshgl)
        enddo
!               deallocate(tmpshp)
!               deallocate(tmpshgl)
c
c... Set Characteristic Size
c    i_spat_var_eps_flag = 0 - no spatially varying eps.  size = 1.0
c                          1 - based on vol**1/3
c                          2 - based on max edge length
c
        if (i_spat_var_eps_flag .eq. 0) then
          elemb_size(:) = 1.0
        elseif (i_spat_var_eps_flag .eq. 1) then
          elemb_size(:) = elemb_vol(:)**(1.0/3.0)
        else
          elemb_size(:) = elemb_size(:)
        endif
c
      return
      end
      
      
c..***********************************************************
c....this small routine just for calculate volume of the element
c*************************************************************
      subroutine e3elsize (blk,shp,shgl,xl,ien,loc_el_size,loc_el_vol)

c                                                                      
c----------------------------------------------------------------------
c This routine calculates characteristic size of each element 
c
c input: 
c  shp    (nen,blk%g) 		: element shape-functions
c  shgl   (nsd,nen,blk%g)	: element local-grad-shape-functions
c  x     (numnp,nsd)		: nodal coordinates at current step
c  ien (blk%e,blk%s)
c output:
c  loc_el_vol(blk%e)		: volume of the element
c  loc_el_size(blk%e)		: maximum node-to-node distance for element
c  
c
c----------------------------------------------------------------------
c
      use eblock
      include "common.h"
      type (LocalBlkData) blk


c....Passed arrays
        dimension shp(blk%s,blk%g), shgl(nsd,blk%s,blk%g),
     &            ien(blk%e,blk%s)

c
c local arrays
c
        real*8 loc_el_vol(blk%e)
	real*8 loc_el_size(blk%e)

        dimension dxidx(blk%e,nsd,nsd),WdetJ(blk%e)

        real*8 tmp(blk%e), disttmp(blk%e)
        dimension dxdxi(blk%e,nsd,nsd)

        dimension sgn(blk%e,blk%s),  shape(blk%e,blk%s),
     &  shdrv(blk%e,nsd,blk%s),   xl(blk%e,blk%n,nsd)

c
c Initialize array
c
        loc_el_vol(:) = zero
c
c Set sign
c
        do i=1,blk%s
          where ( ien(:,i) < 0 )
            sgn(:,i) = -one
          elsewhere
            sgn(:,i) = one
          endwhere
        enddo
c
c*************Localizing coordinates*************
c
!        call localx (blk,x,xl,ien,nsd,'gather  ')
 
c
c.... loop through the integration points
c
        
        do ith = 1, blk%g   

          call getshp(blk,ith, shp,          shgl,      sgn, 
     &                shape,        shdrv)

c
c*************get WdetJ*********************************
c
c
c.... compute the deformation gradient
c
          dxdxi = zero
c
          do n = 1, blk%n
            dxdxi(:,1,1) = dxdxi(:,1,1) + xl(1:blk%e,n,1) * shdrv(:,1,n)
            dxdxi(:,1,2) = dxdxi(:,1,2) + xl(1:blk%e,n,1) * shdrv(:,2,n)
            dxdxi(:,1,3) = dxdxi(:,1,3) + xl(1:blk%e,n,1) * shdrv(:,3,n)
            dxdxi(:,2,1) = dxdxi(:,2,1) + xl(1:blk%e,n,2) * shdrv(:,1,n)
            dxdxi(:,2,2) = dxdxi(:,2,2) + xl(1:blk%e,n,2) * shdrv(:,2,n)
            dxdxi(:,2,3) = dxdxi(:,2,3) + xl(1:blk%e,n,2) * shdrv(:,3,n)
            dxdxi(:,3,1) = dxdxi(:,3,1) + xl(1:blk%e,n,3) * shdrv(:,1,n)
            dxdxi(:,3,2) = dxdxi(:,3,2) + xl(1:blk%e,n,3) * shdrv(:,2,n)
            dxdxi(:,3,3) = dxdxi(:,3,3) + xl(1:blk%e,n,3) * shdrv(:,3,n)
          enddo
c
c.... compute the inverse of deformation gradient
c
          dxidx(:,1,1) =   dxdxi(:,2,2) * dxdxi(:,3,3) 
     &                   - dxdxi(:,3,2) * dxdxi(:,2,3)
          dxidx(:,1,2) =   dxdxi(:,3,2) * dxdxi(:,1,3) 
     &                   - dxdxi(:,1,2) * dxdxi(:,3,3)
          dxidx(:,1,3) =   dxdxi(:,1,2) * dxdxi(:,2,3) 
     &                   - dxdxi(:,1,3) * dxdxi(:,2,2)
          tmp          = one / ( dxidx(:,1,1) * dxdxi(:,1,1) 
     &                         + dxidx(:,1,2) * dxdxi(:,2,1)  
     &                         + dxidx(:,1,3) * dxdxi(:,3,1) )
          dxidx(:,1,1) = dxidx(:,1,1) * tmp
          dxidx(:,1,2) = dxidx(:,1,2) * tmp
          dxidx(:,1,3) = dxidx(:,1,3) * tmp
          dxidx(:,2,1) = (dxdxi(:,2,3) * dxdxi(:,3,1) 
     &                  - dxdxi(:,2,1) * dxdxi(:,3,3)) * tmp
          dxidx(:,2,2) = (dxdxi(:,1,1) * dxdxi(:,3,3) 
     &                  - dxdxi(:,3,1) * dxdxi(:,1,3)) * tmp
          dxidx(:,2,3) = (dxdxi(:,2,1) * dxdxi(:,1,3) 
     &                  - dxdxi(:,1,1) * dxdxi(:,2,3)) * tmp
          dxidx(:,3,1) = (dxdxi(:,2,1) * dxdxi(:,3,2) 
     &                  - dxdxi(:,2,2) * dxdxi(:,3,1)) * tmp
          dxidx(:,3,2) = (dxdxi(:,3,1) * dxdxi(:,1,2) 
     &                  - dxdxi(:,1,1) * dxdxi(:,3,2)) * tmp
          dxidx(:,3,3) = (dxdxi(:,1,1) * dxdxi(:,2,2) 
     &                  - dxdxi(:,1,2) * dxdxi(:,2,1)) * tmp
c
          WdetJ = Qwt(blk%l,ith)/ tmp
c
          do i=1,blk%s
            loc_el_vol(:) = loc_el_vol(:) + abs(shape(:,i)*WdetJ)
          enddo
c
c.... end of the loop over integration points
c
        enddo
c
c... Find maximum distance across nodes for element
c
        loc_el_size(:) = zero
        do n1 = 1, blk%n-1
          do n2 = 2, blk%n
            disttmp(:) = ( (xl(1:blk%e,n1,1) - xl(1:blk%e,n2,1))**2 +
     &                     (xl(1:blk%e,n1,2) - xl(1:blk%e,n2,2))**2 +
     &                     (xl(1:blk%e,n1,3) - xl(1:blk%e,n2,3))**2 )**0.5
            where (disttmp(:) .gt. loc_el_size(:))
              loc_el_size(:) = disttmp(:)
            endwhere
          enddo
        enddo
c
c      write(*,*) "e3elsize:elem vol value=",loc_el_vol(:)  

       return
       end

