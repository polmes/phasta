      subroutine hello
      use omp_lib
      use eblock
      include "common.h"
      type (LocalBlkData) blk

      integer idshared(4)
#ifdef HAVE_OMP
!$OMP  parallel do 
!$OMP& private (ith,id)
!$OMP& shared (idshared)
      do ith = 1, 4
          id = omp_get_thread_num ( )
!$OMP critical
          write (*,*) 'HELLO from rank, process, ith ', myrank,id,ith
!$OMP end critical
          idshared(ith)=id
      enddo
!     do ith = 1, 4
!         write (*,*) 'HELLO from rank, process, ith ', myrank,idshared(ith),ith
!     enddo
#endif
      write (*,*) 'Back to MPI ', myrank


!   stuff that is local to a block is moving to a derived data type
!   nenl,nshl,npro,ngauss -> blkDi

!
! next three lines are in eblock.h
!       type localBdata
!           integer :: nenl(MAXBLK),nshl(MAXBLK),npro(MAXBLK),ngauss(MAXBLK)
!       end type LocalBdata
! next line is in the routine that  does any operation on the data structure
!      type (localBdata) blkDi      ! then, for a given private value of iblk in all routines below elmgmr
!  next few lines describe access to the shared data structure for a block
!      nenl -> blkDi%nenl
!      nshl -> blkDi%nshl     !including
!
!      real*8 u1(blkDi%nenl)  
! assuming that iblk and blkDi where passed to the routine doing this allocation...! I also assume the compiler is not going to choke on
!      nproL=blkDi(%npro)
!      u1(1:nproL)= tmp(1:nproL)*0.5 !assuming tmp is dimensioned like u1

! TESING  First fill the new data structure
       do iblk = 1, nelblk
          iblock = iblk         ! used in local mass inverse (p>2)
          iblkts = iblk         ! used in timeseries
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
          inum   = iel + npro - 1
          ngauss = nint(lcsyst)
          blk%n   = lcblk(5,iblk) ! no. of vertices per element
          blk%s   = lcblk(10,iblk)
          blk%g = nint(lcsyst)
          blk%e   = lcblk(1,iblk+1) - iel
       enddo

! now see if it is accessible when passed
#ifdef HAVE_OMP
!$OMP  parallel do 
!$OMP& private (iblk,blk)
#endif
      do iblk = 1, nelblk
          blk%n   = lcblk(5,iblk) ! no. of vertices per element
          blk%s   = lcblk(10,iblk)
          blk%g = nint(lcsyst)
          blk%e   = lcblk(1,iblk+1) - lcblk(1,iblk)
       call testPass2(iblk, blk)
      enddo


      return
      end


!DEC$ ATTRIBUTES NOINLINE :: testpass2
      subroutine testPass2(iblk,blk)
      use eblock
      include "common.h"
      type (LocalBlkData) blk
      real*8 ui(blk%e,blk%g)
      real*8 yl(blk%e,blk%s,nflow)  
      nelbk=blk%e
      ngqpt=blk%g
      do i=1,ngqpt
          ui(1:nelbk,ngqpt)=0.25*(yl(1:nelbk,1,1)+yl(1:nelbk,2,1)
     &                           +yl(1:nelbk,3,1)+yl(1:nelbk,4,1))
      enddo
      return
      end

      subroutine testPass(iblk,blk)
      use eblock
      include "common.h"
      type (LocalBlkData) blk
      do iblk = 1, nelblk
          iblock = iblk         ! used in local mass inverse (p>2)
          iblkts = iblk         ! used in timeseries
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
          inum   = iel + npro - 1
          ngauss = nint(lcsyst)
          write(*,*) 'nenl ', nenl, blk%n  
          write(*,*) 'nshl ', nshl, blk%s 
          write(*,*) 'ngauss ', ngauss, blk%g
          write(*,*) 'npro ', npro, blk%e 
       enddo

       return
       end
      
