      subroutine rerun_check()

      include "common.h"
      include "mpif.h"

      logical rr_exts


      if(myrank.eq.master) then
        inquire(file='rerun-check',exist=rr_exts)
        if(rr_exts) then
          open(unit=772,file='rerun-check',status='old')
          read(772,*,iostat=ierr2)stopjob
          if(ierr2.gt.0) then
            write(*,*) "read from rerun-check failed in unknown way."
            write(*,*) " Try again to set as you have not yet set a stop step number"
            stopjob=-1
            close(772)
!            close(772,status='delete')  
          else if(ierr2.lt.0) then
            write(*,*) "read from rerun-check found empty file."
            write(*,*) " Try again to set as you have not yet set a stop step number"
            stopjob=-1
            close(772)
!            close(772,status='delete') 
          else
            write(*,*) "read from rerun-check found stop step number=",stopjob
            close(772,status='delete')  !cool,  google taught me a new fortran trick
            if(stopjob.lt.0) then  !Another job wants me to stop 
              istepsPerPhase=-1*stopjob
              iphaseDone=lstep/istepsPerPhase
              stopjob=(iphaseDone+1)*istepsPerPhase
            endif
!            close(772,status='delete')  !cool,  google taught me a new fortran trick
          endif ! end of ierr2 on read
        endif ! end of exist
      endif ! master

      call MPI_BCAST(stopjob,1,MPI_INTEGER,master,
     .               MPI_COMM_WORLD,ierr)

      return
      end


      subroutine doubleruncheck()

      include "common.h"
      include "mpif.h"
  
      logical dr_exts

      stopjob = 0
      

      if(myrank.eq.master) then
        inquire(file='doubleRun-check',exist=dr_exts)
        if(dr_exts) then
          open(unit=772,file='doubleRun-check',status='old')
          read(772,*,iostat=ierr2)seconds_remaining_old
          if(ierr2.gt.0) then
            write(*,*) "read from doubleRun-check failed in unknown way."
            stopjob=0
            close(772)
          else if(ierr2.lt.0) then
            write(*,*) "read from doubleRun-check found empty file."
            stopjob=0
            close(772)
          else ! good inputs were read
            tRemainMe=allocated_seconds 
            if(tRemainMe/seconds_remaining_old .lt. 1.1) then  ! kill me
              stopjob=1
              write(*,*) "old job can do more than new ",stopjob
              close(772)
            else ! get old job to checkpoint and then I will run from there
              stopjob=0
              open(unit=773,file='rerun-check',status='unknown')
              modStop=1
              if(nphasesincycle.ne.0) modStop=nstepsincycle/nphasesincycle
              write(773,*) -1*modStop
              close(773)
!this completes the signal to the old run to run until the step number
!mod with modStop = 0 (e.g., stop on a phase as phase averaging requires
!currently.  Now new program has to put into a wait loop to let the old
!run run those steps and write a checkpoint.  We will know we can proceed
!when doublerun-check no longer exists because the file is deleted at 
!the end of PHASTAS stepping loop (in the original job).
              close(772)
              write(*,*) 'entering do loop'
              do i = 1,1000
#ifdef __bgq__
! FAILS but is supposed to be what new XL compilers want #ifdef __ibmxl__
!WORKS but requires us to maintain flags in CMakesLists.txt #ifdef HAVE_XLCOMPILER
                call sleep_(1)  ! sleep for second between checks for change
#else
                call sleep(1)  ! sleep for second between checks for change
#endif
                inquire(file='doubleRun-check',exist=dr_exts)
                write(*,*) i, dr_exts, TMRC() 
                if(.not.dr_exts) exit
              enddo
              if(dr_exts) then
                 write(*,*) 'first run never cleared dr-check'
              else
                 write(*,*) 'first cleared dr-check'
              endif  ! dr_exts
            endif  ! who gets to continue
          endif  ! good inputs
        endif ! the file existed ?
      endif   ! myrank=0

      call MPI_BCAST(stopjob,1,MPI_INTEGER,master,
     .               MPI_COMM_WORLD,ierr)
      if(myrank.eq.master) write(*,*) ' cleared MPI_BCAST in doublerunCheck'

      return
      end



      subroutine TimeRemainingIntoDoubleRunCheck()

      include "common.h"

      if(myrank.eq.master) then
        open(unit=772,file='doubleRun-check',status='unknown')
        tRunMe=TMRC() -ttim(100)
        timeLeft=allocated_seconds - tRunMe
        write(772,*) timeLeft
!this completes the signal to any new run how much time I can still run
! so that they can decide to to quit or tell me to quit through
! rerun-check
        close(772)
      endif

      return
      end

      subroutine ClearDoubleRunCheck()

      include "common.h"

      if(myrank.eq.master) then
        open(unit=772,file='doubleRun-check',status='unknown')
        close(772,status='delete')  
        open(unit=772,file='rerun-check',status='unknown')
        close(772,status='delete')  
      endif

      return
      end
