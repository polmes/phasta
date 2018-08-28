      subroutine elmStatsLhs( x,  iBC,   iper,  ilwork )
c-----------------------------------------------------------------------
c     compute the necessary terms for the statistics projection
c     matrices.
c-----------------------------------------------------------------------
      use     stats
      use     pointer_data
      use eblock
      
      include "common.h"
      type (LocalBlkData) blk

      real*8  x(numnp,3)
      integer iBC(nshg), iper(nshg), ilwork(nlwork)
      
      real*8, allocatable :: xl(:,:,:)
      real*8, allocatable :: lStsVec(:,:,:)

c
c.... loop over element blocks
c
      stsVec = zero
      
      do iblk = 1, nelblk
         iel    = lcblk(1,iblk)
         lcsyst = lcblk(3,iblk)
         nenl   = lcblk(5,iblk) ! no. of vertices per element
         nshl   = lcblk(10,iblk)
         ndofl  = lcblk(8,iblk)
         npro   = lcblk(1,iblk+1) - iel 
          blk%n   = lcblk(5,iblk) ! no. of vertices per element
          blk%s   = lcblk(10,iblk)
          blk%e   = lcblk(1,iblk+1) - iel
          blk%g = nint(lcsyst)
          blk%l = lcblk(3,iblk)
          blk%o = lcblk(4,iblk)

         allocate ( xl(bsz,nenl,3)             )
         allocate ( lStsVec(bsz,nshl,nResDims) )
c
c.... localize needed data
c
         call localx (blk, x,    xl,  mien(iblk)%p, nsd,   'gather  ' )
c
c.... form the Lhs
c
         call e3StsLhs(blk, xl, lStsVec )
c
c.... assemble
c
         call local (blk,stsVec, lStsVec, mien(iblk)%p,
     &               nResDims, 'scatter ' ) 

         deallocate ( xl       )
         deallocate ( lStsVec  )
c
c.... end loop over element blocks
c
      enddo

      if (numpe > 1) then
        call commu (stsVec, ilwork, nResDims  , 'in ')
      endif
c
c.... local periodic boundary conditions (no communications)
c
      do j = 1,nshg
         if (btest(iBC(j),10)) then
            i = iper(j)
            stsVec(i,:) = stsVec(i,:) + stsVec(j,:)
         endif
      enddo
c
      do i = 1,nshg
         stsVec(i,:) = stsVec(iper(i),:)
      enddo
      if (numpe > 1) then
        call commu (stsVec, ilwork, nResDims  , 'out')
      endif

      return
      end
      
c-----------------------------------------------------------------------
c  Assemble the residual for the statistics
c-----------------------------------------------------------------------
      subroutine elmStatsRes( y,        ac,    u, x,      shp,     shgl, 
     &                        shpb,     shglb,       iBC,     BC, 
     &                        iper,     ilwork,      rowp,    colm)
      
      use     stats
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      
      
      real*8  y(nshg,ndof),             ac(nshg,ndof), x(numnp,nsd),
     &        u(nshg,nsd),
     &        shp(MAXTOP,maxsh,MAXQPT),  shgl(MAXTOP,nsd,maxsh,MAXQPT),
     &        shpb(MAXTOP,maxsh,MAXQPT),
     &        shglb(MAXTOP,nsd,maxsh,MAXQPT),
     &        BC(nshg,ndofBC), res(nshg,ndof) 

      integer iBC(nshg),                iper(nshg),
     &        ilwork(nlwork),           rowp(nshg*nnz),
     &        colm(nshg+1)
      dimension GradV(nshg,nsdsq) 
      

      lhs    = 0
      stsVec = zero
      GradV=zero     
 
      stsResFlg = 1
      ierrcalctmp=ierrcalc ! save the current value of ierrcalc
      ierrcalc=0           ! set it to zero so that we don't calculate error
                           ! here (which would overflow memory around rjunk)
      !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myrank.eq.master) write(*,*) 'calling ElmGMR'
      call ElmGMR (u,         y,     ac,    x,
     &             shp,       shgl,       iBC,       
     &             BC,        shpb,       shglb,
     &             res,       iper,       ilwork,   
     &             rowp,      colm,      
#ifdef HAVE_PETSC
     &             lhsPETSc,
#endif
     &             rerr,       GradV   )      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if(myrank.eq.master) write(*,*) 'after ElmGMR in elmStatsRes'
      stsResFlg = 0
      ierrcalc=ierrcalctmp  ! reset it back to original value
      return 
      end

