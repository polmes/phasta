      module spanStats
      real*8, allocatable, dimension (:,:) :: stsBar

      end module spanStats


      subroutine getStsBar(ilwork, nsons, ifath, iBC)

      use stats
      use spanStats
 
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      dimension nsons(nfath), rinvsons(nfath), tmpStats(nshg,10),
     &          ifath(numnp), ilwork(nlwork), iBC(numnp),
     &          tmpStatsf(nfath,10), tmpStatsft(nfath,10)
      integer periodicity
      
c     Assign the conservative statistics to a temporary array
      tmpStats(:,1) = stsPres(:)
      tmpStats(:,2:4) = stsVel(:,:)
      tmpStats(:,5:10) = stsVelSqr(:,:) ! 11,22,33,12,23,31

c     Zero on processor periodic nodes so will not be added twice
      if(myrank.eq.master) then
         do i=1, numnp
            if(btest(iBC(i),10)) then 
              periodicity = 1
              exit
            endif
         enddo
      endif
      call MPI_BCAST(periodicity,1,MPI_INTEGER,master,
     &               MPI_COMM_WORLD,ierr)

      where(btest(iBC,10).or.btest(iBC,12))
          tmpStats(:,1)=zero
          tmpStats(:,2)=zero
          tmpStats(:,3)=zero
          tmpStats(:,4)=zero
          tmpStats(:,5)=zero
          tmpStats(:,6)=zero
          tmpStats(:,7)=zero
          tmpStats(:,8)=zero
          tmpStats(:,9)=zero
          tmpStats(:,10)=zero
      endwhere

      if (numpe.gt.1) then
            
            numtask = ilwork(1)
            itkbeg = 1
            
c     zero the nodes that are "solved" on the other processors  
            do itask = 1, numtask
               
               iacc   = ilwork (itkbeg + 2)
               numseg = ilwork (itkbeg + 4)
               
               if (iacc .eq. 0) then
                  do is = 1,numseg
                     isgbeg = ilwork (itkbeg + 3 + 2*is)
                     lenseg = ilwork (itkbeg + 4 + 2*is)
                     isgend = isgbeg + lenseg - 1
                     tmpStats(isgbeg:isgend,:) = zero
                  enddo
               endif
               
               itkbeg = itkbeg + 4 + 2*numseg
               
            enddo !itask
            
      endif  ! numpe.gt.1

      tmpStatsf = zero

c     accumulate sum of sons to the fathers
      do i = 1,numnp
         ifathi=ifath(i)
         tmpStatsf(ifathi,1:10) = tmpStatsf(ifathi,1:10) 
     &                             + tmpStats(i,1:10)            
      enddo

c     Now  the true fathers and serrogates combine results and update
c     each other.
c     
      if(numpe .gt. 1) then
         call drvAllreduce(tmpStatsf, tmpStatsft,nfath*10)
      else
         tmpStatsft=tmpStatsf
      endif
c     divide by # of sons to get average father for this step
c
      if (periodicity.eq.1.and.nohomog.eq.1) then
         rinvsons = one/(nsons-one)   ! division is expensive
      else
         rinvsons = one/nsons   ! division is expensive
      endif
      tmpStatsft(:,1) = tmpStatsft(:,1) * rinvsons(:) 
      tmpStatsft(:,2) = tmpStatsft(:,2) * rinvsons(:)
      tmpStatsft(:,3) = tmpStatsft(:,3) * rinvsons(:)
      tmpStatsft(:,4) = tmpStatsft(:,4) * rinvsons(:)
      tmpStatsft(:,5) = tmpStatsft(:,5) * rinvsons(:)
      tmpStatsft(:,6) = tmpStatsft(:,5) * rinvsons(:)
      tmpStatsft(:,7) = tmpStatsft(:,5) * rinvsons(:)
      tmpStatsft(:,8) = tmpStatsft(:,5) * rinvsons(:)
      tmpStatsft(:,9) = tmpStatsft(:,5) * rinvsons(:)
      tmpStatsft(:,10) = tmpStatsft(:,5) * rinvsons(:)


c     Add to running time average
      den = max(1,lstep-istartSpanAvg)
      tfact = one/den
      stsBar(:,:)=tfact*tmpStatsft(:,:)+(one-tfact)*stsBar(:,:) 


      end subroutine getStsBar


      subroutine wstsBar(ifail)

      use spanStats
      include "common.h"
  
      character*20 fname1,fmt1
      character*5 cname
      integer itmp, irstou

      itmp = 1
      if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
      write (fmt1,"('(''stsbar.'',i',i1,',1x)')") itmp
      write (fname1,fmt1) lstep
      fname1 = trim(fname1) // cname(myrank+1)
c     
      open (unit=irstou, file=fname1, status='unknown',
     &     err=996)
              
      write (irstou,*) nfath, lstep
      if((itwmod.gt.0) .or. (irscale.ge.0)
     &   .or. (ispanAvg.eq.1)) then
           do i=1,nfath            
               write (irstou,*) stsBar(i,1),stsBar(i,2),stsBar(i,3),
     &                          stsBar(i,4),stsBar(i,5),stsBar(i,6),
     &                          stsBar(i,7),stsBar(i,8),stsBar(i,9),
     &                          stsBar(i,10)
           enddo
      endif
      close (irstou)

      return

996   call error ('stsBar  ','opening ', irstou)

      end subroutine wstsBar
