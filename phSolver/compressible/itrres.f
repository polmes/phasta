        subroutine ItrRes (yp,        yc,        x,
     &                     shp,       shgl,      
     &                     iBC,       BC,        shpb,
     &                     shglb,     rmes,      iper,
     &                     ilwork,    ac)
c
c----------------------------------------------------------------------
c
c This routine calculates the modified residual vector.
c
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        use pointer_data
        use eblock
c
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

      type (LocalBlkData) blk
c
        dimension yp(nshg,nflow),             yc(nshg,ndof),
     &            x(numnp,nsd),             ac(nshg,ndof), 
     &            iBC(nshg),                BC(nshg,ndofBC), 
     &            rmes(nshg,nflow),          ilwork(nlwork),
     &            iper(nshg)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT)
    
        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)

        ttim(81) = ttim(81) - secs(0.0)

c
c.... -------------------->   interior elements   <--------------------
c
        jump  = 0
        ires  = 2
        iprec = 0
c
c.... loop over the element-blocks
c
        do iblk = 1, nelblk
c
c.... set up the parameters
c$$$c
c$$$          iel    = lcblk(1,iblk)
c$$$          nenl   = lcblk(5,iblk)
c$$$          mattyp = lcblk(7,iblk)
c$$$          ndofl  = lcblk(8,iblk)
c$$$          npro   = lcblk(1,iblk+1) - iel 
c
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel
          ngauss = nint(lcsyst)
          blk%n   = lcblk(5,iblk) ! no. of vertices per element
          blk%s   = lcblk(10,iblk)
          blk%e   = lcblk(1,iblk+1) - iel  
          blk%l = lcblk(3,iblk)
          blk%g = nint(blk%l)
          blk%o = lcblk(4,iblk)
          blk%i = lcblk(1,iblk)

c
c
c.... compute and assemble the residuals and the preconditioner
c
          allocate (tmpshp(nshl, MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))

          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)

          call AsIRes (blk, yp,                      yc,
     &                 mxl(iblk)%p, mdwl(iblk)%p,    mxmudmi(iblk)%p,
     &                 tmpshp,                  tmpshgl,
     &                 mien(iblk)%p,            mmat(iblk)%p,
     &                 rmes,                    ac)

          deallocate (tmpshp)
          deallocate (tmpshgl)

c
c.... end of interior element loop
c
        enddo
c
c.... -------------------->   boundary elements   <--------------------
c
        if (Navier .eq. 1 .and. Jactyp.ne.0) then
c
c.... loop over the elements
c
        do iblk = 1, nelblb
c
c.... set up the parameters
c$$$c
c$$$          iel    = lcblkb(1,iblk)
c$$$          nenl   = lcblkb(5,iblk)
c$$$          nenbl  = lcblkb(6,iblk)
c$$$          mattyp = lcblkb(7,iblk)
c$$$          ndofl  = lcblkb(8,iblk)
c$$$          npro   = lcblkb(1,iblk+1) - iel 
c$$$c
c
          iel    = lcblkb(1,iblk)
          lelCat = lcblkb(2,iblk)
          lcsyst = lcblkb(3,iblk)
          iorder = lcblkb(4,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          mattyp = lcblkb(7,iblk)
          ndofl  = lcblkb(8,iblk)
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          npro   = lcblkb(1,iblk+1) - iel 
          if(lcsyst.eq.3) lcsyst=nenbl
          ngaussb = nintb(lcsyst)
          blk%n   = lcblkb(5,iblk) ! no. of vertices per element
          blk%s   = lcblkb(10,iblk)
          blk%e   = lcblkb(1,iblk+1) - iel  
          blk%l = lcblkb(3,iblk)
          blk%g = nintb(blk%l)
          blk%o = lcblkb(4,iblk)
          blk%i = lcblkb(1,iblk)

c
          allocate (tmpshpb(nshl,MAXQPT))
          allocate (tmpshglb(nsd,nshl,MAXQPT))
          
          tmpshpb(1:nshl,:) = shpb(lcsyst,1:nshl,:)
          tmpshglb(:,1:nshl,:) = shglb(lcsyst,:,1:nshl,:)
          
c
c.... compute and assemble the residuals
c

          call AsBRes (blk, yp,   yc,     mxlb(iblk)%p, mdwl(iblk)%p,
     &                 tmpshpb,           tmpshglb,
     &                 mienb(iblk)%p,     mmatb(iblk)%p,
     &                 miBCB(iblk)%p,     mBCB(iblk)%p,
     &                 rmes)
c

          deallocate (tmpshpb)
          deallocate (tmpshglb)
c.... end of boundary element loop
c
        enddo

        endif

        ttim(81) = ttim(81) + secs(0.0)

c
c.... ---------------------->   communications  <-----------------------
c
        if((iabc==1)) !are there any axisym bc's
     &       call rotabc(rmes(1,2), iBC, 'in ')
c
        if (numpe > 1) then
           call commu (rmes, ilwork, nflow, 'in ')
        endif
        
c
c.... ---------------------->   post processing  <----------------------
c
c.... satisfy the BCs on the modified residual
c
        call bc3Res (yc,  iBC,  BC,  rmes, iper, ilwork)
c
c.... return
c
        return
        end
