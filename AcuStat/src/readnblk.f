c  readnblk.f (pronounce "Reed and Block Dot Eff") contains:
c
c    module readarrays ("Red Arrays") -- contains the arrays that
c     are read in from binary files but not immediately blocked 
c     through pointers.
c
c    subroutine readnblk ("Reed and Block") -- allocates space for
c     and reads data to be contained in module readarrays.  Reads
c     all remaining data and blocks them with pointers.
c


      module readarrays
      
      real*8, allocatable :: point2x(:,:)
      real*8, allocatable :: qold(:,:)
      real*8, allocatable :: dwal(:)
      real*8, allocatable :: errors(:,:)
      real*8, allocatable :: ybar(:,:)
      real*8, allocatable :: uold(:,:)
      real*8, allocatable :: acold(:,:)
      integer, allocatable :: iBCtmp(:)
      real*8, allocatable :: BCinp(:,:)

      integer, allocatable :: point2ilwork(:)
      integer, allocatable :: nBC(:)
      integer, allocatable :: point2iper(:)
      integer, allocatable :: point2ifath(:)
      integer, allocatable :: point2nsons(:)
      
      end module



      subroutine readnblk
c
      use readarrays
      include "commonAcuStat.h"
      include "mpif.h"
c
      real*8, allocatable :: xread(:,:),qread(:,:)
      real*8, allocatable :: qread1(:) 
      real*8, allocatable :: uread(:,:), acread(:,:)
      real*8, allocatable :: BCinpread(:,:)
      real*8 globmax,globmin
      integer, allocatable :: iperread(:), iBCtmpread(:)
      integer, allocatable :: ilworkread(:), nBCread(:)
      character*10 cname2
      character*30 fmt1
      character*255 fname1,fnamer,fnamelr
      character*255 warning
      
      integer :: descriptor, color, nfiles, nfields
      integer ::  numparts, nppf
      character*255 fname2, temp2
      character*64 temp1
      integer :: igeom, ibndc, irstin, ierr
      integer :: ndof, ndoferrors, ndofybar, ndofyphbar
      integer :: itmp, itmp2

      integer intfromfile(50) ! integers read from headers
      logical exinput

      integer :: numts, numphavg, its, isize, sumts, totsumts, debug
      integer :: sumtsQ, totsumtsQ
      integer :: ierror, idwal, iphavg, ivort, ndofvort
      integer :: descriptorG,GPID
      integer :: timestep, ithree
      integer :: nshg2,ndof2
      integer, allocatable :: tablets(:,:)
      real*8, allocatable :: qybarcumul(:,:),qerrorcumul(:,:)
      real*8, allocatable :: qyphbarcumul(:,:,:), qvorticity(:,:)

      rQInit=0
      rInit=0
      rOpen=0
      rRead=0
      rFnlz=0
      rWrite=0

      debug = 0

      numts=-1

      if(myrank == 0) then 
        fnamer='AcuStat_input.dat'
        fnamer = trim(fnamer) // char(0)
        inquire(file=fnamer,exist=exinput)
        if(exinput) then 
          open(unit=24,file=fnamer,status='old')
          read(24,*) numts
          read(24,*) ierror
          read(24,*) numphavg
          read(24,*) ivort
          read(24,*) idwal
        else
           write(*,*) 'ERROR: Input file ',
     &                 trim(fnamer),' does not exist!'
        endif
      endif

      call mpi_barrier(mpi_comm_world, ierr)
      call mpi_bcast(exinput,1,MPI_LOGICAL,0,mpi_comm_world,ierr)
      if(.not. exinput) then ! AcuStat_input.dat does not exist. Quit
        return
      endif

      call MPI_Bcast(numts,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
      allocate(tablets(numts,3)) ; tablets=-1
      
      if(myrank == 0 ) then
        ! Read the series of time steps
        do its=1,numts
          read(24,*) tablets(its,1), tablets(its,2),tablets(its,3)
        enddo
        close(24)
 
      ! Print some info
        write(*,*) 'There are ',numts, 'restart files to read for ybar'
        write(*,*) 'The solution field from the last time step will be',
     &             ' added'

        if(ierror .gt. 0) then
          write(*,*) 'The error field will also be accumulated'
          ierror = 1 ! security
        else
          write(*,*) 'The error field will NOT be accumulated'
        endif

        if(numphavg .gt. 0) then
          write(*,*) 'The phase average fields (',numphavg,
     &                ') will also be accumulated'
        else
          write(*,*) 'The phase average fields will NOT be accumulated'
        endif

        if(ivort .gt. 0) then
          write(*,*) 'The vorticity field will also be added'
          ivort = 1 ! security
        else
          write(*,*) 'The vorticity field will NOT be added'
        endif

        if(idwal .gt. 0) then
          write(*,*) 'The dwal field will also be added'
          idwal = 1 ! security
        else
          write(*,*) 'The dwal field will NOT be added'
        endif

      endif

      isize=3*numts
      call MPI_Bcast(tablets,isize,MPI_INTEGER,0,MPI_COMM_WORLD, ierr) 
      call MPI_Bcast(ierror,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
      call MPI_Bcast(numphavg,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
      call MPI_Bcast(idwal,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
      call MPI_Bcast(ivort,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)

!      write(*,*) 'my rank is: ',myrank,numts

!      do irank=0,7
!        if (irank == myrank) then
!          do its=1,numts
!            write(*,*) myrank,tablets(its,1:2)
!          enddo
!        endif
!        call mpi_barrier(MPI_COMM_WORLD,ierr)
!      enddo

      nfiles = nsynciofiles
      numparts = numpe
      color = int(myrank/(numparts/nfiles))
      sumts = 0
      sumtsQ = 0
      totsumts = 0
      totsumtsQ = 0

      ithree=3

!     Loop over the restart files to read and accumulate the ybar field
      do its=1,numts

!
!       Read the table 
!
        timestep = tablets(its,1)
        sumts = tablets(its,2)
        sumtsQ = tablets(its,3)
        totsumts = totsumts + sumts
        totsumtsQ=totsumtsQ + sumtsQ

!
!       Get the right restart files
!
        itmp=1
        if (timestep .gt. 0) itmp = int(log10(float(timestep+1)))+1
        write (fmt1,"('(''restart-dat.'',i',i1,',1x)')") itmp
        write (fnamer,fmt1) timestep
        fnamer = trim(fnamer) // cname2(color+1)

        if(debug==1) then
         write(*,*)myrank,timestep,nfiles,numparts,trim(adjustl(fnamer))
        endif

!
!       query, init and open
!
        
        rdelta= TMRC()
        call queryphmpiio(fnamer//char(0), nfields, nppf);
        rQInit=rQInit+ TMRC()-rdelta
        if (myrank == 0) then
          write(*,*) ''
          write(*,*) 'Considering time step',timestep,sumts,sumtsQ
          write(*,*) 'Number of fields in restart-dat: ',nfields
          write(*,*) 'Number of parts per file restart-dat: ',nppf
        endif

        rdelta= TMRC()
        call initphmpiio(nfields,nppf,nfiles,descriptor,'read'//char(0))
        rInit=rInit+ TMRC()-rdelta
        rdelta= TMRC()
        call openfile( fnamer//char(0), 'read'//char(0), descriptor )
        rOpen=rOpen+ TMRC()-rdelta

        rdelta= TMRC()
!
!       Read the ybar field
!
        itmp = int(log10(float(myrank+1)))+1
        write (temp1,"('(''ybar@'',i',i1,',A1)')") itmp
        write (fname1,temp1) (myrank+1),'?'
        fname1 = trim(fname1)

        intfromfile=0
        call readheader(descriptor,fname1//char(0),intfromfile,ithree,
     &                  'integer'//char(0),iotype)
        if(intfromfile(1).ne.0) then   
          nshg2=intfromfile(1) !Do not use nshg and ndof from common.h here!
          ndofybar=intfromfile(2)
          !lstep2=intfromfile(3)
        else
          write(*,*)'ERROR: Did not find the field requested 
     &               in AcuStat_input',myrank
        endif

        ! Allocate some memory for the first ts only
        if(its==1) then
          allocate( qybarcumul(nshg2,ndofybar) ) ; qybarcumul = 0.d0
          ndofybar1=ndofybar
        endif
            

        allocate( qread(nshg2,ndofybar)  ) ; qread = 0.d0
        iqsiz = nshg2*ndofybar
        call readdatablock(descriptor,fname1//char(0),qread,iqsiz,
     &                     'double'//char(0),iotype) 

        ! Weight correctly the field and accumulate for average 
        rRead=rRead+ TMRC()-rdelta
        rdelta= TMRC()
        if(ndofybar.eq.ndofybar1) then ! we are as big as the biggest so post Q
          qybarcumul(:,1:ndofybar-1) = qybarcumul(:,1:ndofybar-1) + qread (:,1:ndofybar-1)*sumts
          qybarcumul(:,ndofybar) = qybarcumul(:,ndofybar) + qread (:,ndofybar)*sumtsQ
        else ! these are pre-Q files 
          qybarcumul(:,1:ndofybar) = qybarcumul(:,1:ndofybar) + qread (:,1:ndofybar)*sumts
        endif
        rStat=rStat+ TMRC()-rdelta
        deallocate(qread)

!
!       Read the errors
!
        if(ierror == 1) then
          rdelta= TMRC()
          write (temp1,"('(''errors@'',i',i1,',a1)')")
     &             itmp
          write (fname1,temp1) (myrank+1),'?'
          fname1 = trim(fname1)
          intfromfile=0
          call readheader(descriptor,fname1//char(0),intfromfile,ithree,
     &                    'integer'//char(0),iotype)
          if(intfromfile(1).ne.0) then   
            nshg2=intfromfile(1)
            ndoferrors=intfromfile(2)
            !lstep=intfromfile(3)
          else
            write(*,*)'ERROR: Did not find the field requested 
     &                 in AcuStat_input',myrank
          endif

          ! Allocate some memory for the first ts only
          if(its==1) then
            allocate( qerrorcumul(nshg2,ndoferrors) ) ; qerrorcumul = 0.d0
          endif

          allocate( qread(nshg2,ndoferrors) ) ; qread = 0.d0
          iqsiz=nshg2*ndoferrors
          call readdatablock(descriptor,fname1//char(0),qread,iqsiz,
     &                         'double'//char(0),iotype)
          rRead=rRead+ TMRC()-rdelta
          rdelta= TMRC()

          ! Weight correctly the field and accumulate for average 
          qerrorcumul(:,:) = qerrorcumul(:,:) + qread(:,:)*sumts
          rStat=rStat+ TMRC()-rdelta
          deallocate(qread)
        endif ! if ierror == 1

!
!       Read the phase_average fields
!
        if(numphavg .gt. 0) then
          do iphavg = 1,numphavg
            rdelta= TMRC()
            itmp = int(log10(float(myrank+1)))+1
            itmp2 = int(log10(float(iphavg)))+1
!            write (temp1,"('(''phase_average@'',i',i1,',A1)')") itmp
!            write(*,*) temp1  !   ('phase_average@',i1,A1)
            write (temp1,
     &              "('(''phase_average'',i',i1,',''@'',i',i1,',A1)')")
     &               itmp2, itmp
!            write(*,*) temp1  !   ('phase_average',i1,'@',i1,A1)
            write (fname1,temp1) iphavg,(myrank+1),'?'
            fname1 = trim(fname1)
            if(debug == 1) then
              if(myrank == 0 .or. myrank  == numpe-1) then
                write(*,*)'Rank: ',myrank,'- phase: ',
     &                     iphavg,' - field: ',trim(fname1)
              endif
            endif
            intfromfile=0
            call readheader(descriptor,fname1//char(0),intfromfile,
     &                  ithree,'integer'//char(0),iotype)
            if(intfromfile(1).ne.0) then   
              nshg2=intfromfile(1) !Do not use nshg and ndof from common.h here!
              ndofyphbar=intfromfile(2)
              !lstep2=intfromfile(3)
            else
              write(*,*)'ERROR: Did not find the field requested 
     &                 in AcuStat_input',myrank
            endif

            ! Allocate some memory for the first ts only
            if(its==1 .and. iphavg==1) then
              allocate( qyphbarcumul(nshg2,ndofyphbar,numphavg) )
              ndofyphbar1=ndofyphbar
              qyphbarcumul = 0.d0
            endif

            allocate( qread(nshg2,ndofyphbar)  ) ; qread = 0.d0
            iqsiz = nshg2*ndofyphbar
            call readdatablock(descriptor,fname1//char(0),qread,iqsiz,
     &                     'double'//char(0),iotype) 
            rRead=rRead+ TMRC()-rdelta
            rdelta= TMRC()
  
            if((myrank.eq.0).and.(iphavg.eq.1)) then
          write(*,*)'ndofybar,ndofyphbar,totsumts,totsumtsQ,
     &sumts,sumtsq'
          write(*,*)ndofybar,ndofyphbar,totsumts,totsumtsQ,sumts,sumtsq
            endif
            ! Weight correctly the field and accumulate for average 
        if(ndofybar.eq.ndofybar1) then ! we are as big as the biggest but last may have different weight
          qyphbarcumul(:,1:ndofyphbar-1,iphavg) = qyphbarcumul(:,1:ndofyphbar-1,iphavg) 
     &                               + qread (:,1:ndofyphbar-1)*sumts
          qyphbarcumul(:,ndofyphbar,iphavg) = qyphbarcumul(:,ndofyphbar,iphavg) 
     &                               + qread (:,ndofyphbar)*sumtsQ
        else ! these are pre-Q files in mixed case..this one is small so sumts covers all 
          qyphbarcumul(:,1:ndofyphbar,iphavg) = qyphbarcumul(:,1:ndofyphbar,iphavg) 
     &                               + qread (:,1:ndofyphbar)*sumts
        endif

            rStat=rStat+ TMRC()-rdelta
            deallocate(qread)
          enddo
        endif
!
!       If last ts, read also the solution, error and dwal field 
!       for reduction with M2N or for convenience
!

        if(its == numts) then
          rdelta= TMRC()
!
!         Read the solution of the last time step 
!
          itmp = int(log10(float(myrank+1)))+1
          write (temp1,"('(''solution@'',i',i1,',A1)')") itmp
          write (fname1,temp1) (myrank+1),'?'
          fname1 = trim(fname1)

          intfromfile=0
          call readheader(descriptor,fname1//char(0) ,intfromfile,
     &                    ithree,'integer'//char(0), iotype)

          if(intfromfile(1).ne.0) then 
            nshg2=intfromfile(1)
            ndof2=intfromfile(2)
            allocate( qold(nshg2,ndof2) )
            qold(:,:) = -9.87654321e32
            lstep=intfromfile(3)
            allocate( qread(nshg2,ndof2) ) ; qread = 0.d0
            iqsiz=nshg2*ndof2
            call readdatablock(descriptor,fname1//char(0),qread,iqsiz,
     &                         'double'//char(0),iotype)
            qold(1:nshg2,1:ndof2) = qread(1:nshg2,1:ndof2)
            deallocate(qread)
          else
            write(*,*)'ERROR: Did not find the field requested 
     &               in AcuStat_input',myrank
          endif

!
!         Read the vorticity field from the last time step
!
          if(ivort == 1) then 
            write (temp1,"('(''vorticity@'',i',i1,',a1)')")
     &             itmp
            write (fname1,temp1) (myrank+1),'?'
            fname1 = trim(fname1)
            intfromfile=0
            call readheader(descriptor,fname1//char(0),intfromfile,
     &                    ithree,'integer'//char(0),iotype)
            if(intfromfile(1).ne.0) then 
              nshg2=intfromfile(1)
              ndofvort=intfromfile(2)
              allocate( qvorticity(nshg2,ndofvort) )
              qvorticity(:,:) = -9.87654321e32
              allocate( qread(nshg2,ndofvort) ) ; qread = 0.d0
              iqsiz=nshg2*ndofvort
              call readdatablock(descriptor,fname1//char(0),qread,
     &                      iqsiz,'double'//char(0),iotype)
              qvorticity(1:nshg2,1:ndofvort) = qread(1:nshg2,1:ndofvort)
              deallocate(qread)
            else
              write(*,*)'ERROR: Did not find the field requested 
     &                 in AcuStat_input',myrank
            endif
          endif

!
!         Read the dwal
!
          if(idwal == 1) then
            write (temp1,"('(''dwal@'',i',i1,',a1)')")
     &             itmp
            write (fname1,temp1) (myrank+1),'?'
            fname1 = trim(fname1)
            intfromfile=0
            call readheader(descriptor,fname1//char(0),intfromfile,
     &                    ithree,'integer'//char(0),iotype)
            if(intfromfile(1).ne.0) then 
              nshg2=intfromfile(1)
              allocate( dwal(nshg2) )
              dwal(:) = -9.87654321e32
              allocate( qread1(nshg2) ) ; qread1 = 0.d0
              iqsiz=nshg2*1
              call readdatablock(descriptor,fname1//char(0),qread1,
     &                      iqsiz,'double'//char(0),iotype)
              dwal(1:nshg2) = qread1(1:nshg2)
              deallocate(qread1)
            else
              write(*,*)'ERROR: Did not find the field requested 
     &                 in AcuStat_input',myrank
            endif
          endif !if idwal == 1

            rRead=rRead+ TMRC()-rdelta

        endif ! First time step
!
!        Close and finalize
!
        rdelta= TMRC()
        call closefile( descriptor, 'read'//char(0) )
        rClose=rClose+ TMRC()-rdelta
  
        rdelta= TMRC()
        call finalizephmpiio( descriptor ) 
        rFnlz=rFnlz+ TMRC()-rdelta

      enddo !numts = number of restart files to read
      if(myrank.eq.0) then
        write(*,*) 'totsumts,totsumtsQ'
        write(*,*) totsumts,totsumtsQ
      endif

!
!     Weight correctly the average field
!
      rdelta= TMRC()
      if(totsumts.eq.totsumsQ) then 
        if(myrank.eq.0) write(*,*) 'no mixed weights'
        qybarcumul(:,1:ndofybar) = qybarcumul(:,1:ndofybar)/totsumts
      else ! more fields in the some files but that means different samples
        if(myrank.eq.0) write(*,*) ' mixed weights',totsumts, totsumtsQ
        qybarcumul(:,1:ndofybar1-1) = qybarcumul(:,1:ndofybar1-1)/totsumts
        qybarcumul(:,ndofybar1) = qybarcumul(:,ndofybar1)/totsumtsQ
      endif
      if (ierror == 1) then
        qerrorcumul(:,:) = qerrorcumul(:,:)/totsumts
      endif
      if(numphavg .gt. 0) then
       do iphavg = 1,numphavg
         if(totsumts.eq.totsumsQ)  then
           qyphbarcumul(:,1:ndofyphbar,iphavg) = qyphbarcumul(:,1:ndofyphbar,iphavg)/totsumts
         else ! more fields in the some files but that means different samples
           qyphbarcumul(:,1:ndofyphbar1-1,iphavg) = qyphbarcumul(:,1:ndofyphbar1-1,iphavg)/totsumts
           qyphbarcumul(:,ndofyphbar1,iphavg) = qyphbarcumul(:,ndofyphbar1,iphavg)/totsumtsQ
         endif
       enddo
      endif
 
!      Write qold, dwal, errors and qybarcumul (= cumulated ybar)

      nsynciofieldswriterestart = 2 + ierror + numphavg + ivort + idwal ! 2 for solution and accumulated ybar
!      write(*,*) 'nsynciofieldswriterestart:',nsynciofieldswriterestart 
      
!      call Write_AcuStat(myrank, tablets(1,1), tablets(numts,1),
!     &          totsumts, nshg2, ndof2, ndoferrors, ndofybar,
!     &          qold, dwal, qerrorcumul, qybarcumul)

      call mpi_barrier(mpi_comm_world, ierr)
      rStat=rStat+ TMRC()-rdelta
      rdelta= TMRC()
      call Write_AcuStat(myrank, tablets(1,1), tablets(numts,1),
     &          totsumts, nshg2, ndof2, qold)

      if(ivort == 1) then 
        call write_field(myrank,'a','vorticity',9,qvorticity,
     &                       'd', nshg2, ndofvort, tablets(numts,1))
      endif

      if (ierror == 1) then
        call write_error(myrank, tablets(numts,1), nshg2, ndoferrors,
     &                   qerrorcumul)
      endif

      call write_field(myrank,'a', 'ybar', 4, 
     &          qybarcumul, 'd', nshg2, ndofybar1, tablets(numts,1))

      if(numphavg .gt. 0) then
        do iphavg = 1,numphavg
          call write_phavg2(myrank,'a','phase_average',13,iphavg,
     &         numphavg,qyphbarcumul(:,:,iphavg),'d',nshg2,ndofyphbar1,
     &         tablets(numts,1))
        enddo
      endif

      if (idwal == 1) then
        call write_field(myrank, 'a', 'dwal', 4, dwal,
     &          'd', nshg2, 1, tablets(numts,1))
      endif
      rWrite=rWrite+ TMRC()-rdelta
      if(myrank.eq.0) then
         write(*,*) 'Time in QueryInit',rQInit
         write(*,*) 'Time in Init',rInit
         write(*,*) 'Time in Open',rOpen
         write(*,*) 'Time in Read',rRead
         write(*,*) 'Time in Stat',rStat
         write(*,*) 'Time in Fnlz',rFnlz
         write(*,*) 'Time in Write',rWrite
      endif
      
!
!     Free memory
!
      deallocate(tablets)
      deallocate(qybarcumul, qold)
      if(ierror == 1) then
        deallocate(qerrorcumul)
      endif
      if(numphavg .gt. 0) then  
        deallocate(qyphbarcumul)
      endif
      if(ivort == 1) then
        deallocate(qvorticity)
      endif
      if(idwal == 1) then
        deallocate(dwal)
      endif
      
      return
 
      end
