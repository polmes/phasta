c-----------------------------------------------------------------------
c
c     module for time averaged statistics (conservative projection).
c
c-----------------------------------------------------------------------
      module stats
      
      integer nResDims, nSolDims, nLhsDims, nTimeStep, stsResFlg
      integer stsCompFreq, stsWriteFreq, stsResetFreq, step1,
     &        stsType
      
      real*8, allocatable :: stsVec(:,:)
      
      real*8, allocatable :: stsReg(:)
      real*8, allocatable :: stsMInv(:,:)
      real*8, allocatable :: stsB(:,:)
      real*8, allocatable :: stsDInv(:,:)
      real*8, allocatable :: stsCInv(:,:)
      
      real*8, allocatable :: stsPres(:), stsPresSqr(:), stsVel(:,:),
     &                       stsVelSqr(:,:), stsVelReg(:,:),
     &                       stsStress(:,:), stsVelSqInst(:,:)

      end module stats

c-----------------------------------------------------------------------
c
c     module for spanwise averaged statistics (conservative projection).
c
c-----------------------------------------------------------------------
      module spanStats
      
      real*8, allocatable, dimension (:,:) :: stsBar
      real*8, allocatable, dimension (:,:) :: stsBarKeq


      end module spanStats


