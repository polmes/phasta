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
     &                       stsStress(:,:), stsVelSqInst(:,:),
     &                       stsVelInst(:,:)

      end module stats

c-----------------------------------------------------------------------
c
c     module for spanwise averaged statistics (conservative projection).
c
c-----------------------------------------------------------------------
      module spanStats
      
      real*8, allocatable, dimension (:,:) :: velbar
      real*8, allocatable, dimension (:,:) :: stsBar
      real*8, allocatable, dimension (:,:) :: stsBarKeq
      real*8, allocatable, dimension (:,:) :: velf
      real*8, allocatable, dimension (:,:) :: velftG
      real*8, allocatable, dimension (:,:) :: tmpStatsf
      real*8, allocatable, dimension (:,:) :: tmpStatsftG
      real*8, allocatable, dimension (:,:) :: tmpKeqf
      real*8, allocatable, dimension (:,:) :: tmpKeqftG
      integer, allocatable, dimension (:) :: locifath
      integer, allocatable, dimension (:) :: ifathG
      integer, allocatable, dimension (:) :: rcounts
      integer, allocatable, dimension (:) :: displs

      integer locnfath, stacksz

      end module spanStats


