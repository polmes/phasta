#define PETSC_USE_FORTRAN_MODULES 1
      module petsc_data
#ifdef HAVE_PETSC

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscviewer.h>
#include <petsc/finclude/petscviewer.h90>
!      PetscOffset poff
!      PetscOffset poffs
       PetscInt, save:: petsc_bs,petsc_o, petsc_t, petsc_PA 

       PetscViewer :: viewer

      Mat,save :: lhsPETSc,lhsPs
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

