      subroutine ElmGMR (u,         y,         ac,        x,     
     &                     shp,       shgl,      iBC,
     &                     BC,        shpb,      shglb,
     &                     res,       iper,      ilwork,
     &                     rowp,      colm,     
#ifdef HAVE_PETSC
     &                     lhsPETSc, ! cannot call it lhsP here because leslib uses that
#endif
     &                     rerr,     GradV)
c
c----------------------------------------------------------------------
c
c This routine computes the LHS mass matrix, the RHS residual 
c vector, and the preconditioning matrix, for use with the GMRES
c solver.
c
c Zdenek Johan, Winter 1991.      (Fortran 90)
c Chris Whiting, Winter 1998.     (Matrix EBE-GMRES)
c Alberto Figueroa, Winter 2004.  CMM-FSI
c Irene Vignon, Spring 2004.
c----------------------------------------------------------------------
c
      use solvedata  ! brings in lhs16
      use pvsQbi  ! brings in NABI
      use stats   !  
      use pointer_data  ! brings in the pointers for the blocked arrays
      use local_mass
      use timedata
      use eblock
c
      include "common.h"

c
        dimension y(nshg,ndof),         ac(nshg,ndof),
     &            u(nshg,nsd),
     &            x(numnp,nsd),               
     &            iBC(nshg),           
     &            BC(nshg,ndofBC),  
     &            res(nshg,nflow),
     &            iper(nshg)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
        dimension qres(nshg,idflx),     rmass(nshg)
        dimension GradV(nshg,nsdsq)
c
        dimension ilwork(nlwork)

        integer rowp(nshg*nnz),         colm(nshg+1)


        real*8, allocatable, dimension(:,:,:,:,:) :: xlhs
        real*8, allocatable, dimension(:,:,:,:) :: rl, rerrl,StsVecl

        real*8  rerr(nshg,10)

        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)

        real*8 spmasstot(20),  ebres(nshg)
        real*8 cfl(nshg), CFLfl_maxtmp
        integer icflhits(nshg)
c
c.... set up the timer
c

c
c.... -------------------->   diffusive flux   <--------------------
c
c.... set up parameters
c
        ires   = 1

        if (idiff==1 .or. idiff==3 .or. isurf==1) then ! global reconstruction
                                                       ! of qdiff
c
c loop over element blocks for the global reconstruction
c of the diffusive flux vector, q, and lumped mass matrix, rmass
c
        qres = zero
        if(iter == nitr .and. icomputevort == 1 ) then
          GradV = zero
        endif
        rmass = zero

      nshlc=lcblk(10,1) ! set to first block and maybe all blocks if monotop.
      allocate (tmpshp(nshlc,MAXQPT))
      allocate (tmpshgl(nsd,nshlc,MAXQPT))
      lcsyst= lcblk(3,1)
      tmpshp(1:nshlc,:) = shp(lcsyst,1:nshlc,:)
      tmpshgl(:,1:nshlc,:) = shgl(lcsyst,:,1:nshlc,:)
      
        do iblk = 1, nelblk
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk) ! no. of vertices per element
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
          if(blk%s.ne.nshlc) then  ! never true in monotopology but makes code 
            nshlc=blk%s
            deallocate (tmpshp)
            deallocate (tmpshgl)
            allocate (tmpshp(blk%s,MAXQPT))
            allocate (tmpshgl(nsd,blk%s,MAXQPT))
            tmpshp(1:blk%s,:) = shp(blk%l,1:blk%s,:)
            tmpshgl(:,1:blk%s,:) = shgl(blk%l,:,1:blk%s,:)
          endif
c     
c.... compute and assemble diffusive flux vector residual, qres,
c     and lumped mass matrix, rmass


          if(iter == nitr .and. icomputevort == 1 ) then
            call AsIqGradV (blk,y,                x,
     &                   tmpshp,
     &                   tmpshgl,
     &                   mien(iblk)%p,
     &                   GradV)
          endif
          call AsIq (blk,y,                x,                       
     &                   tmpshp,
     &                   tmpshgl,
     &                   mien(iblk)%p,     mxmudmi(iblk)%p,  
     &                   qres,             rmass )
        enddo
       
c
c.... form the diffusive flux approximation
c
        call qpbc( rmass, qres, iBC, iper, ilwork )       
        if(iter == nitr .and. icomputevort == 1 ) then
          call solveGradV( rmass, GradV, iBC, iper, ilwork )
        endif
        deallocate (tmpshp)
        deallocate (tmpshgl)
c
      endif 
c
c.... -------------------->   interior elements   <--------------------
c
      res    = zero
      if (stsResFlg .ne. 1) then
        flxID = zero
      endif

      if ((usingpetsc.eq.0).and.(lhs .eq. 1)) then
        lhs16   = zero
      endif
c
c.... loop over the element-blocks
c
c
c.... allocate the element matrices
c
      
#ifdef HAVE_OMP
      BlockPool=8
#else
      BlockPool=1
#endif
      nshlc=lcblk(10,1) ! set to first block and maybe all blocks if monotop.
      allocate ( rl (bsz,nshlc,ndof,BlockPool) )
      allocate (tmpshp(nshlc,MAXQPT))
      allocate (tmpshgl(nsd,nshlc,MAXQPT))
      lcsyst= lcblk(3,1)
      tmpshp(1:nshlc,:) = shp(lcsyst,1:nshlc,:)
      tmpshgl(:,1:nshlc,:) = shgl(lcsyst,:,1:nshlc,:)
      if (lhs .eq. 1) then
        allocate ( xlhs(bsz,16,nshlc,nshlc,BlockPool) )
      endif
      if ( ierrcalc .eq. 1 ) allocate ( rerrl (bsz,nshlc,6,BlockPool) )
      if ( stsResFlg .eq. 1 ) allocate ( StsVecl (bsz,nshlc,nResDims,BlockPool) )
#ifdef HAVE_OMP
      do iblko = 1, nelblk, BlockPool
        rdelta=TMRC() 
!$OMP parallel do
!$OMP& private (ith,iblk,blk,nshc)
        do iblk = iblko,iblko+BlockPool-1
         if(iblk.le.nelblk) then
          ith=1+iblk-iblko
# else
      do iblk = 1, nelblk
          rdelta=TMRC() 
          ith=1
#endif
          iblock = iblk         ! used in local mass inverse (p>2)
          iblkts = iblk         ! used in timeseries
!          iel    = lcblk(1,iblk)
!          lelCat = lcblk(2,iblk)
!          lcsyst = lcblk(3,iblk)
!          iorder = lcblk(4,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
!          inum   = iel + npro - 1
          blk%b   = iblk
          blk%t   = ith
          blk%n   = lcblk(5,iblk) ! no. of vertices per element
          blk%s   = lcblk(10,iblk)
          blk%e   = lcblk(1,iblk+1) - lcblk(1,iblk) 
          blk%l = lcblk(3,iblk)
          blk%g = nint(blk%l)
          blk%o = lcblk(4,iblk)
          blk%i = lcblk(1,iblk)
          if(blk%s.ne.nshlc) then  ! never true in monotopology but makes code 
            nshlc=blk%s
            deallocate (rl)
            allocate ( rl (bsz,blk%s,ndof,BlockPool) )
            deallocate (tmpshp)
            deallocate (tmpshgl)
            allocate (tmpshp(blk%s,MAXQPT))
            allocate (tmpshgl(nsd,blk%s,MAXQPT))
            tmpshp(1:blk%s,:) = shp(blk%l,1:blk%s,:)
            tmpshgl(:,1:blk%s,:) = shgl(blk%l,:,1:blk%s,:)
            if (lhs .eq. 1) then
              deallocate (xlhs)   ! below (local) easier if blk%s is correct size
              allocate ( xlhs(bsz,16,blk%s,blk%s,BlockPool) )
            endif
            if ( ierrcalc .eq. 1 ) then
              deallocate (rerrl)
              allocate ( rerrl (bsz,blk%s,6,BlockPool) )
            endif
            if ( stsResFlg .eq. 1 ) then
              deallocate(StsVecl)
              allocate ( StsVecl (bsz,blk%s,nResDims,BlockPool) )
            endif
          endif   ! different topology endif
c
c.... compute and assemble the residual and tangent matrix
c

          call AsIGMR (blk, y,                   ac,
     &                 x,                   mxmudmi(iblk)%p,      
     &                 tmpshp,
     &                 tmpshgl,
     &                 mien(iblk)%p,
     &                 rl(:,:,:,ith),
     &                 qres,              
     &                 xlhs(:,:,:,:,ith),   rerrl(:,:,:,ith), 
     &                 StsVecl(:,:,:,ith),
     &                 cfl,                 icflhits )
#ifdef HAVE_OMP
         endif ! this is the skip if threads available but blocks finished
        enddo !threaded loop closes here
      rdelta = TMRC() - rdelta
      rthreads = rthreads + rdelta
      rdelta = TMRC() 

        iblkStop=min(nelblk, iblko+BlockPool-1)
        do iblk = iblko,iblkStop
          ith=1+iblk-iblko
          blk%n   = lcblk(5,iblk) ! no. of vertices per element
          blk%s   = lcblk(10,iblk)
          blk%e   = lcblk(1,iblk+1) - lcblk(1,iblk) 
          blk%g = nint(lcsyst)
          blk%l = lcblk(3,iblk)
          blk%o = lcblk(4,iblk)
#else
      rdelta = TMRC() - rdelta
      rthreads = rthreads + rdelta
      rdelta = TMRC() 
#endif
c
c.... assemble the residual
c
          call local (blk,res,    rl(:,:,:,ith),     mien(iblk)%p,    nflow,  'scatter ')
c
c.... assemble the statistics residual
c
          if ( stsResFlg .eq. 1 ) then 
            call local( blk, stsVec, StsVecl(:,:,:,ith), mien(iblk)%p, nResDims, 'scatter ')
          endif
          if ( ierrcalc .eq. 1 ) then
            call local (blk, rerr, rerrl(:,:,:,ith),  
     &                  mien(iblk)%p, 6, 'scatter ')
          endif
          nshl=blk%s  !nshl still used in these routines...temp solution
          npro=blk%e  !npro still used in these routines...temp solution
          if (impl(1) .ne. 9 .and. lhs .eq. 1) then
            if(ipord.eq.1) 
     &        call bc3lhs (iBC, BC,mien(iblk)%p, xlhs(:,:,:,:,ith))
            if(usingpetsc.eq.1) then
#ifdef HAVE_PETSC
              call fillsparsecpetsci (mieng(iblk)%p, 
     &                                xlhs(:,:,:,:,ith),lhsPETSc) 
#else
              write(*,*) 'requested unavailable PETSc'
              call error('elmgmr', 'no PETSc', usingpetc)
#endif
            else
              call fillsparseI16 (mien(iblk)%p, 
     &                 xlhs(:,:,:,:,ith) ,            lhs16,
     &                 rowp,                      colm)
            endif
          endif
c
c.... end of interior element loop assembly
c
#ifdef HAVE_OMP
        enddo ! inner thread assembly loop
#endif
        rdelta = TMRC() - rdelta
        rassembly = rassembly + rdelta

      enddo ! outer loop (only if flat MPI
      deallocate ( rl  )
      deallocate (tmpshp)
      deallocate (tmpshgl)
      if(lhs.eq.1) then
        deallocate ( xlhs )
      endif
      if ( ierrcalc .eq. 1 )   deallocate ( rerrl  )
      if ( stsResFlg .eq. 1 )          deallocate ( StsVecl  )
c
c.... add in lumped mass contributions if needed
c
      if((flmpr.ne.0).or.(flmpl.ne.0)) then
        write(*,*) 'not checked for blk'
        call lmassadd(ac,res,rowp,colm,lhs16,gmass)
      endif

       have_local_mass = 1

c  Divide CFL number by the number of contributors to get
c  the average CFL number.
       cfl = cfl / icflhits

c
c Find element with worst (largest) CFL number
c
        CFLfl_max = cfl(1)
        iCFLfl_maxelem = 1
        do i = 1, nshg
          if (cfl(i) .gt. CFLfl_max) then
            CFLfl_max = cfl(i)
            iCFLfl_maxelem = i
          endif
        enddo
c
c Communicate with other processors to get maximum CFL number over
c all processors
c
        if(numpe. gt. 1) then
           call MPI_ALLREDUCE (CFLfl_max, CFLfl_maxtmp, 1,
     &       MPI_DOUBLE_PRECISION,MPI_MAX, MPI_COMM_WORLD,ierr)
        else
           CFLfl_maxtmp = CFLfl_max
        endif
        CFLfl_max = CFLfl_maxtmp 
c
c.... time average statistics
c       
      if ( stsResFlg .eq. 1 ) then

          if (numpe > 1) then
             call commu (stsVec, ilwork, nResDims  , 'in ')
          endif
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
          
       endif
c
c.... -------------------->   boundary elements   <--------------------
c
c.... loop over the boundary elements
c
        do iblk = 1, nelblb
c
c.... set up the parameters
c
          iel    = lcblkb(1,iblk)
          lelCat = lcblkb(2,iblk)
          lcsyst = lcblkb(3,iblk)
          iorder = lcblkb(4,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          mattyp = lcblkb(7,iblk)
          ndofl  = lcblkb(8,iblk)
          npro   = lcblkb(1,iblk+1) - iel 


          if(lcsyst.eq.3) lcsyst=nenbl
c
          if(lcsyst.eq.3 .or. lcsyst.eq.4) then
             ngaussb = nintb(lcsyst)
          else
             ngaussb = nintb(lcsyst)
          endif
! note the following is the volume element characteristics as needed for 
! routines like getdiff, getdiffsclr, getshpb which were already converted for
! threading and thus need to dimension based on this data.  
! debug to confirm

          blk%n   = lcblkb(5,iblk) ! no. of vertices per element
          blk%s   = lcblkb(9,iblk)
          blk%e   = lcblkb(1,iblk+1) - iel 
          blk%g = nintb(lcsyst)
          blk%l = lcblkb(3,iblk)
          blk%o = lcblkb(4,iblk)
          blk%i = lcblkb(1,iblk)
c
c.... allocate the element matrices
c
!disable          allocate ( xKebe(npro,9,nshl,nshl) )
!disable          allocate ( xGoC (npro,4,nshl,nshl) )
          
c
c.... compute and assemble the residuals corresponding to the 
c     boundary integral
c
          allocate (tmpshpb(nshl,MAXQPT))
          allocate (tmpshglb(nsd,nshl,MAXQPT))
          
          tmpshpb(1:nshl,:) = shpb(lcsyst,1:nshl,:)
          tmpshglb(:,1:nshl,:) = shglb(lcsyst,:,1:nshl,:)

          call AsBMFG (blk,u,                       y,
     &                 ac,                      x,
     &                 tmpshpb,
     &                 tmpshglb,
     &                 mienb(iblk)%p,           mmatb(iblk)%p,
     &                 miBCB(iblk)%p,           mBCB(iblk)%p,
     &                 res,                     xlhsdisabled)

c
c.... satisfy (again, for the vessel wall contributions) the BC's on the implicit LHS
c
c.... first, we need to make xGoC zero, since it doesn't have contributions from the 
c.... vessel wall elements

!disable          xGoC = zero
!disable
!disable          if (impl(1) .ne. 9 .and. lhs .eq. 1) then
!disable             if(ipord.eq.1)
!disable     &         call bc3lhs (iBC, BC,mienb(iblk)%p, xlhsdisabled)
!disable             call fillsparseI (mienb(iblk)%p,
!disable     &                 xlhsdisabled,           lhsK,
!disable     &                 xGoC,             lhsP,
!disable     &                 rowp,                      colm)
!disable          endif

!disable          deallocate ( xKebe )
!disable          deallocate ( xGoC )
          deallocate (tmpshpb)
          deallocate (tmpshglb)
c
c.... end of boundary element loop
c
       enddo

       if(ipvsq.ge.1) then
c
c....  pressure vs. resistance boundary condition sets pressure at
c      outflow to linearly increase as flow through that face increases
c      (routine is at bottom of this file)
c
          call ElmpvsQ (res,y,-1.0d0)     
       endif
           
c
c before the commu we need to rotate the residual vector for axisymmetric
c boundary conditions (so that off processor periodicity is a dof add instead
c of a dof combination).  Take care of all nodes now so periodicity, like
c commu is a simple dof add.
c
       if(iabc==1)              !are there any axisym bc's
     &       call rotabc(res, iBC,  'in ')
c
c
c.... -------------------->   communications <-------------------------
c

       if (numpe > 1) then
          call commu (res  , ilwork, nflow  , 'in ')
       endif

c
c.... ---------------------->   post processing  <----------------------
c
c.... satisfy the BCs on the residual
c
      call bc3Res (iBC,  BC,  res,  iper, ilwork)
c
c.... return
c
c      call timer ('Back    ')
      return
      end


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!********************************************************************
!--------------------------------------------------------------------

      subroutine ElmGMRSclr (y,         ac,        x,     
     &                       shp,       shgl,      iBC,
     &                       BC,        shpb,      shglb,
     &                       res,       iper,      ilwork,
     &                       rowp,      colm,
     &                       cfl      
#ifdef HAVE_PETSC
     &                       ,lhsPs)
#else
     &                       )
#endif
c
c----------------------------------------------------------------------
c
c This routine computes the LHS mass matrix, the RHS residual 
c vector, and the preconditioning matrix, for use with the GMRES
c solver.
c
c----------------------------------------------------------------------
c
      use solvedata
      use pointer_data
      use local_mass
      use eblock
c
      include "common.h"
        include "mpif.h"
c
        dimension y(nshg,ndof),         ac(nshg,ndof),
     &            x(numnp,nsd),         iBC(nshg),           
     &            BC(nshg,ndofBC),      res(nshg),
     &            iper(nshg)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
        dimension qres(nshg,nsd),     rmass(nshg)
c
        integer ilwork(nlwork), rowp(nshg*nnz),   colm(nshg+1)


        real*8, allocatable, dimension(:,:,:) :: xSebe
        real*8  cfl(nshg), CFLls_maxtmp, cflold(nshg)
        integer icflhits(nshg)
        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)
c
c.... set up the timer
c

c
c.... -------------------->   diffusive flux   <--------------------
c
      ires   = 1

      if (idiff==1 .or. idiff==3) then ! global reconstruction of qdiff
c
c loop over element blocks for the global reconstruction
c of the diffusive flux vector, q, and lumped mass matrix, rmass
c
        qres = zero
        rmass = zero
        icflhits = 0
        cfl = zero

      nshlc=lcblk(10,1) ! set to first block and maybe all blocks if monotop.
      allocate (tmpshp(nshlc,MAXQPT))
      allocate (tmpshgl(nsd,nshlc,MAXQPT))
      lcsyst= lcblk(3,1)
      tmpshp(1:nshlc,:) = shp(lcsyst,1:nshlc,:)
      tmpshgl(:,1:nshlc,:) = shgl(lcsyst,:,1:nshlc,:)
        
        do iblk = 1, nelblk
          iel    = lcblk(1,iblk)
          lcsyst = lcblk(3,iblk)
          nenl   = lcblk(5,iblk) ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          npro   = lcblk(1,iblk+1) - iel 
          ngauss = nint(lcsyst)
              
          blk%n   = lcblk(5,iblk) ! no. of vertices per element
          blk%s   = lcblk(10,iblk)
          blk%e   = lcblk(1,iblk+1) - iel 
          blk%g = nint(lcsyst)
          blk%l = lcblk(3,iblk)
          blk%o = lcblk(4,iblk)
          if(blk%s.ne.nshlc) then  ! never true in monotopology but makes code 
            nshlc=blk%s
            deallocate (tmpshp)
            deallocate (tmpshgl)
            allocate (tmpshp(blk%s,MAXQPT))
            allocate (tmpshgl(nsd,blk%s,MAXQPT))
            tmpshp(1:blk%s,:) = shp(blk%l,1:blk%s,:)
            tmpshgl(:,1:blk%s,:) = shgl(blk%l,:,1:blk%s,:)
          endif
c     
c.... compute and assemble diffusive flux vector residual, qres,
c     and lumped mass matrix, rmass

              call AsIqSclr (blk,y,                   x,                       
     &                       tmpshp, 
     &                       tmpshgl, 
     &                       mien(iblk)%p,     qres,                   
     &                       rmass, cfl, icflhits )
       
           enddo

c  Divide CFL number by the number of contributors to get
c  the average CFL number.
           cfl = cfl / icflhits

c
c.... form the diffusive flux approximation
c
           call qpbcSclr ( rmass, qres, iBC, iper, ilwork )       

            deallocate (tmpshp)
            deallocate (tmpshgl)
c
        endif 
c
c.... -------------------->   interior elements   <--------------------
c
        res    = zero
        spmass = zero

        if (lhs .eq. 1) then
           lhsS   = zero
        endif

        if ((impl(1)/10) .eq. 0) then   ! no flow solve so flxID was not zeroed
           flxID = zero
        endif

        icflhits = 0
        cflold = cfl
        cfl = zero
      nshlc=lcblk(10,1) ! set to first block and maybe all blocks if monotop.
      allocate (tmpshp(nshlc,MAXQPT))
      allocate (tmpshgl(nsd,nshlc,MAXQPT))
      lcsyst= lcblk(3,1)
      tmpshp(1:nshlc,:) = shp(lcsyst,1:nshlc,:)
      tmpshgl(:,1:nshlc,:) = shgl(lcsyst,:,1:nshlc,:)
c
c.... loop over the element-blocks
c
        do iblk = 1, nelblk
          iblock = iblk         ! used in local mass inverse (p>2)
          iel    = lcblk(1,iblk)
          lcsyst = lcblk(3,iblk)
          nenl   = lcblk(5,iblk) ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          ndofl  = lcblk(8,iblk)
          npro   = lcblk(1,iblk+1) - iel

          ngauss = nint(lcsyst)
          blk%n   = lcblk(5,iblk) ! no. of vertices per element
          blk%s   = lcblk(10,iblk)
          blk%e   = lcblk(1,iblk+1) - iel 
          blk%g = nint(lcsyst)
          blk%l = lcblk(3,iblk)
          blk%o = lcblk(4,iblk)
          if(blk%s.ne.nshlc) then  ! never true in monotopology but makes code 
            nshlc=blk%s
            deallocate (tmpshp)
            deallocate (tmpshgl)
            allocate (tmpshp(blk%s,MAXQPT))
            allocate (tmpshgl(nsd,blk%s,MAXQPT))
            tmpshp(1:blk%s,:) = shp(blk%l,1:blk%s,:)
            tmpshgl(:,1:blk%s,:) = shgl(blk%l,:,1:blk%s,:)
          endif
c
c.... allocate the element matrices
c
          allocate ( xSebe(bsz,nshl,nshl) )
c
c.... compute and assemble the residual and tangent matrix
c
          call AsIGMRSclr(blk,y,                   ac,
     &                 x,
     &                 tmpshp,
     &                 tmpshgl,
     &                 mien(iblk)%p,        res,
     &                 qres,                xSebe, mxmudmi(iblk)%p,
     &                 cfl,  icflhits, cflold )
c
c.... satisfy the BC's on the implicit LHS
c     
          if (impl(1) .ne. 9 .and. lhs .eq. 1) then
            if(usingpetsc.eq.1) then
#ifdef HAVE_PETSC
              if(ipord.eq.1) 
     &          call bc3LHSSclr (iBC, BC,mien(iblk)%p, xSebe)
              call fillsparsecpetscs( mieng(iblk)%p, xSebe, lhsPs)
#else
               write(*,*) 'requested unavailable PETSc'
               call error('elmgmrsclr', 'no PETSc', usingpetc)
#endif 
            else
              call fillsparseSclr (mien(iblk)%p, 
     &                 xSebe,             lhsS,
     &                 rowp,              colm)
            endif
          endif

          deallocate ( xSebe )
c
c.... end of interior element loop
c
       enddo
       deallocate (tmpshp)
       deallocate (tmpshgl)

c
c.... add in lumped mass contributions if needed
c
       if(usingpetsc.eq.0) then
       if((flmpr.ne.0).or.(flmpl.ne.0)) then
          call lmassaddSclr(ac(:,isclr), res,rowp,colm,lhsS,gmass)
       endif
       endif

       have_local_mass = 1

c  Divide CFL number by the number of contributors to get
c  the average CFL number.
       cfl = cfl / icflhits

c
c Find element with worst (largest) CFL number
c
        CFLls_max = cfl(1)
        iCFLls_maxelem = 1
        do i = 1, nshg
          if (cfl(i) .gt. CFLls_max) then
            CFLls_max = cfl(i)
            iCFLls_maxelem = i
          endif
        enddo
c
c Communicate with other processors to get maximum CFL number over
c all processors
c
        if(numpe. gt. 1) then
           call MPI_ALLREDUCE (CFLls_max, CFLls_maxtmp, 1,
     &       MPI_DOUBLE_PRECISION,MPI_MAX, MPI_COMM_WORLD,ierr)
        else
           CFLls_maxtmp = CFLls_max
        endif
        CFLls_max = CFLls_maxtmp

c
c
c  call DtN routine which updates the flux to be consistent with the
c  current solution values.  We will put the result in the last slot of
c  BC (we added a space in input.f).  That way we can localize this
c  value to the boundary elements.  This is important to keep from calling
c  the DtN evaluator more than once per node (it can be very expensive).
c
         if(idtn.eq.1)  call DtN(iBC,BC,y)
c
c.... -------------------->   boundary elements   <--------------------
c
c
c.... loop over the boundary elements
c
        do iblk = 1, nelblb
c
c.... set up the parameters
c
          iel    = lcblkb(1,iblk)
          lcsyst = lcblkb(3,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          ndofl  = lcblkb(8,iblk)
          npro   = lcblkb(1,iblk+1) - iel

          if(lcsyst.eq.3) lcsyst=nenbl
          if(lcsyst.eq.3 .or. lcsyst.eq.4) then
             ngaussb = nintb(lcsyst)
          else
             ngaussb = nintb(lcsyst)
          endif
! note the following is the volume element characteristics as needed for 
! routines like getdiff, getdiffsclr, getshpb which were already converted for
! threading and thus need to dimension based on this data.  
! debug to confirm

          blk%n   = lcblkb(5,iblk) ! no. of vertices per element
          blk%s   = lcblkb(9,iblk)
          blk%e   = lcblkb(1,iblk+1) - iel 
          blk%g = nintb(lcsyst)
          blk%l = lcblkb(3,iblk)
          blk%o = lcblkb(4,iblk)
c
c localize the dtn boundary condition
c

          if(idtn.eq.1)   call dtnl(   iBC, BC, mienb(iblk)%p,
     &              miBCB(iblk)%p,  mBCB(iblk)%p)

c
c.... compute and assemble the residuals corresponding to the 
c     boundary integral
c
          allocate (tmpshpb(nshl,MAXQPT))
          allocate (tmpshglb(nsd,nshl,MAXQPT))
          
          tmpshpb(1:nshl,:) = shpb(lcsyst,1:nshl,:)
          tmpshglb(:,1:nshl,:) = shglb(lcsyst,:,1:nshl,:)

          call AsBSclr (blk,y,                       x,
     &                  tmpshpb, 
     &                  tmpshglb, 
     &                  mienb(iblk)%p,           mmatb(iblk)%p,
     &                  miBCB(iblk)%p,           mBCB(iblk)%p,
     &                  res)

          deallocate(tmpshpb)
          deallocate(tmpshglb)
c
c.... end of boundary element loop
c
        enddo
c
c
c.... -------------------->   communications <-------------------------
c

       if (numpe > 1) then
        call commu (res  , ilwork, 1  , 'in ')
       endif

c
c.... ---------------------->   post processing  <----------------------
c
c.... satisfy the BCs on the residual
c
      call bc3ResSclr (iBC,  res,  iper, ilwork)
c
c.... return
c
      return
      end

        
c
c....routine to compute and return the flow rates for coupled surfaces of a given type
c        
      subroutine GetFlowQ (qsurf,y,srfIdList,numSrfs)
        
      use pvsQbi  ! brings in NABI
c
      include "common.h"
      include "mpif.h"

      real*8  y(nshg,3)
      real*8  qsurf(0:MAXSURF),qsurfProc(0:MAXSURF)
      integer numSrfs, irankCoupled, srfIdList(0:MAXSURF)

c note we only need the first three entries (u) from y

      qsurfProc=zero
      
      do i = 1,nshg
      
        if(numSrfs.gt.zero) then
          do k = 1,numSrfs
            irankCoupled = 0
            if (srfIdList(k).eq.ndsurf(i)) then
              irankCoupled=k
              do j = 1,3              
                 qsurfProc(irankCoupled) = qsurfProc(irankCoupled)
     &                            + NABI(i,j)*y(i,j)
              enddo
              exit
            endif      
          enddo       
        endif
      
      enddo
c      
c     at this point, each qsurf has its "nodes" contributions to Q
c     accumulated into qsurf. Note, because NABI is on processor this
c     will NOT be Q for the surface yet
c
c.... reduce integrated Q for each surface, push on qsurf
c
       npars=MAXSURF+1
       if(impistat.eq.1) then
         iAllR = iAllR+1
       elseif(impistat.eq.2) then
          iAllR = iAllR+1
       endif
       if(impistat2.eq.1) call MPI_BARRIER (MPI_COMM_WORLD, ierr)
       if(impistat.gt.0) rmpitmr = TMRC()
       call MPI_ALLREDUCE (qsurfProc, qsurf(:), npars,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
       if(impistat.eq.1) then 
         rAllR = rAllR+TMRC()-rmpitmr
       elseif(impistat.eq.2) then
         rAllRScal = rAllRScal+TMRC()-rmpitmr
       endif
  
c
c.... return
c
      return
      end


        
c
c... routine to couple pressure with flow rate for each coupled surface
c
      subroutine ElmpvsQ (res,y,sign)     

      use pvsQbi  ! brings in NABI
      use convolImpFlow !brings in the current part of convol coef for imp BC
      use convolRCRFlow !brings in the current part of convol coef for RCR BC

      include "common.h"
      include "mpif.h"

      real*8 res(nshg,ndof),y(nshg,3)
      real*8 p(0:MAXSURF)
      integer irankCoupled

c
c... get p for the resistance BC
c           
      if(numResistSrfs.gt.zero) then
        call GetFlowQ(p,y,nsrflistResist,numResistSrfs)  !Q pushed into p but at this point 
                          ! p is just the full Q for each surface
        p(:)=sign*p(:)*ValueListResist(:) ! p=QR  now we have the true pressure on each
                                        ! outflow surface.  Note sign is -1
                                        ! for RHS, +1 for LHS
c
c....  multiply it by integral NA n_i
c     
       do i = 1,nshg
          do k = 1,numResistSrfs
              irankCoupled = 0
              if (nsrflistResist(k).eq.ndsurf(i)) then 
                  irankCoupled=k
                  res(i,1:3)=res(i,1:3) + p(irankCoupled)*NABI(i,1:3)     
                  exit 
              endif
          enddo   
       enddo
       
      endif !end of coupling for Resistance BC

      
c
c... get p for the impedance BC
c     
      if(numImpSrfs.gt.zero) then
        call GetFlowQ(p,y,nsrflistImp,numImpSrfs)  !Q pushed into p but at this point 
                          ! p is just the full Q for each surface
        do j = 1,numImpSrfs
            if(sign.lt.zero) then ! RHS so -1
               p(j)= sign*(poldImp(j) + p(j)*ImpConvCoef(ntimeptpT+2,j))  !pressure p=pold+ Qbeta
            elseif(sign.gt.zero) then ! LHS so sign is positive
                p(j)= sign*p(j)*ImpConvCoef(ntimeptpT+2,j)
            endif
        enddo
             
c
c....  multiply it by integral NA n_i
c     
       do i = 1,nshg
          do k = 1,numImpSrfs
              irankCoupled = 0
              if (nsrflistImp(k).eq.ndsurf(i)) then 
                  irankCoupled=k
                  res(i,1:3)=res(i,1:3) + p(irankCoupled)*NABI(i,1:3)      
                  exit
              endif
          enddo   
       enddo
       
      endif !end of coupling for Impedance BC
c
c... get p for the RCR BC
c     
      if(numRCRSrfs.gt.zero) then
        call GetFlowQ(p,y,nsrflistRCR,numRCRSrfs)  !Q pushed into p but at this point 
                          ! p is just the full Q for each surface
        do j = 1,numRCRSrfs
            if(sign.lt.zero) then ! RHS so -1
                p(j)= sign*(poldRCR(j) + p(j)*RCRConvCoef(lstep+2,j)) !pressure p=pold+ Qbet
                p(j)= p(j) - HopRCR(j) ! H operator contribution 
            elseif(sign.gt.zero) then ! LHS so sign is positive
                p(j)= sign*p(j)*RCRConvCoef(lstep+2,j)
            endif
        enddo
             
c
c....  multiply it by integral NA n_i
c     
       do i = 1,nshg
          do k = 1,numRCRSrfs
              irankCoupled = 0
              if (nsrflistRCR(k).eq.ndsurf(i)) then 
                  irankCoupled=k
                  res(i,1:3)=res(i,1:3) + p(irankCoupled)*NABI(i,1:3)
                  exit
              endif
          enddo   
       enddo
       
      endif !end of coupling for RCR BC

      return
      end







