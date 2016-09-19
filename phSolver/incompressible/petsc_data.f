      module petsc_data
#ifdef USE_PETSC

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscpc.h"
#include "finclude/petscksp.h"
!      PetscOffset poff
!      PetscOffset poffs
      Mat,save :: lhsP,lhsPs
!      PC,save :: pc,pcs
      !PC pc
!      KSP,save :: ksp,ksps
!      Vec,save :: DyP, resP, DyPLocal, resPGlobal
!      Vec,save :: DyPs, resPs, DyPLocals, resPGlobals
!      PCType ptype
!      PetscErrorCode ierr
!      PetscInt, save:: PetscOne, PetscRow, PetscCol, LocalRow, LocalCol
!      IS,save :: LocalIndexSet, GlobalIndexSet
!      IS,save :: LocalIndexSets, GlobalIndexSets
!      ISLocalToGlobalMapping,save :: VectorMapping
!      ISLocalToGlobalMapping,save :: GblVectorMapping
!      ISLocalToGlobalMapping,save :: VectorMappings
!      ISLocalToGlobalMapping,save :: GblVectorMappings
!      VecScatter,save :: scatter
!      VecScatter,save :: scatters
      integer, save :: firstpetsccall = 1
!      integer, save :: firstpetsccalls = 1
#endif
      end module

