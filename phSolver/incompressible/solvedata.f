      module solvedata

      integer nsolflow,npermDims, nTmpDims, nPermDimsS, nTmpDimsS
#ifdef SP_LHS
      real*4,  allocatable :: lhs16(:,:), lhsP(:,:), lhsS(:,:)
#else
      real*8,  allocatable :: lhs16(:,:), lhsP(:,:), lhsS(:,:)
#endif
#ifdef SP_Solve
      real*4,  allocatable :: aperm(:,:), atemp(:,:)
      real*4,  allocatable :: apermS(:,:,:), atempS(:,:)
#else
      real*8,  allocatable :: aperm(:,:), atemp(:,:)
      real*8,  allocatable :: apermS(:,:,:), atempS(:,:)
#endif

      end module

c-----------------------------------------------------------------------
c allocate the arrays
c-----------------------------------------------------------------------
      subroutine aSDf 
      use solvedata
      include "common.h"
      allocate (lhs16(16,nnz_tot))
      allocate (lhsP(4,nnz_tot))
      if(leslib.eq.1) then
        allocate (aperm(nshg,nPermDims))
        allocate (atemp(nshg,nTmpDims))
        aperm=zero
        atemp=zero
      endif
      return
      end
   
      subroutine aSDs(nsclrsol) 
      use solvedata
      include "common.h"
      allocate (lhsS(nnz_tot,nsclrsol))
      if(leslib.eq.1) then
        allocate (apermS(nshg,nPermDimsS,nsclrsol))
        allocate (atempS(nshg,nTmpDimsS))
        apermS=zero
        atempS=zero
      endif
      return
      end
c-----------------------------------------------------------------------
c delete the arrays
c-----------------------------------------------------------------------
      subroutine dSDf 
      use solvedata
      include "common.h"
      if(usingpetsc.ne.1) then
         deallocate (lhs16)
         deallocate (lhsP)
      endif
      if(leslib.eq.1) then
        deallocate (aperm)
        deallocate (atemp)
        nPermDims=0
      endif
      return
      end

      subroutine dSDs 
      use solvedata
      include "common.h"
      deallocate (lhsS)
      if(leslib.eq.1) then
        deallocate (apermS)
        deallocate (atempS)
      endif
      return
      end
