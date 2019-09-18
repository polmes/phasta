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
      real*8, allocatable :: uold(:,:)
      real*8, allocatable :: acold(:,:)
      integer, allocatable :: iBCtmp(:)
      real*8, allocatable :: BCinp(:,:)
      integer, allocatable :: point2ilwork(:)
!      integer, allocatable :: fncorp(:)
      integer, allocatable :: twodncorp(:,:)
      integer, allocatable :: nBC(:)
      integer, allocatable :: point2iper(:)
      integer, target, allocatable :: point2ifath(:)
      integer, target, allocatable :: point2nsons(:)
      
      end module

      subroutine readnblk
      use iso_c_binding 
      use readarrays
      use fncorpmod
      use phio
      use phiotimer
      use phstr
      use syncio
      use posixio
      use streamio
      use STG_BC
      use spanStats
      use lesArrs
      include "common.h"
      include "mpif.h"

      real*8, target, allocatable :: xread(:,:), qread(:,:), acread(:,:)
      real*8, target, allocatable :: STGread(:,:), BCrestartRead(:,:)
      real*8, target, allocatable :: uread(:,:)
      real*8, target, allocatable :: BCinpread(:,:)
      real*8, target, allocatable :: cdelsqread(:,:)
      real*8, target, allocatable :: lesnutread(:,:)
      real*8 :: iotime
      integer, target, allocatable :: iperread(:), iBCtmpread(:)
      integer, target, allocatable :: ilworkread(:), nBCread(:)
      integer, target, allocatable :: fncorpread(:)
      integer, target, allocatable :: nsonsread(:)
      integer fncorpsize
      integer irandVars, iSTGsiz
      integer iVelbsiz, nfath2, nflow2, istssiz
      character*5 cname
      character*10 cname2, cname2nd
      character*8 mach2
      character*30 fmt1
      character*25 fmtv,fmts,fmtk,fmtI
      character*255 fname1,fnamer,fnamelr
      character*25 fnamev,fnames,fnamek,fnameI
      character*255 warning
      integer igeomBAK, ibndc, irstin, ierr
      integer, target :: intfromfile(50) ! integers read from headers
      integer :: descriptor, descriptorG, GPID, color, nfields
      integer ::  numparts, nppf
      integer :: ierr_io, numprocs, itmp, itmp2
      integer :: ignored
      integer :: fileFmt
      integer :: readstep
      character*255 fname2, temp2
      character*64 temp1
      type(c_ptr) :: handle
      character(len=1024) :: dataInt, dataDbl

      integer ny, nx
      real*8, allocatable :: ypoints(:)
      logical exlog

      real*8, allocatable :: tpW(:,:),tpW2(:,:),tpC(:,:),tpC2(:,:)

      dataInt = c_char_'integer'//c_null_char
      dataDbl = c_char_'double'//c_null_char
c
c.... determine the step number to start with
c
      irstart=0
      if (myrank.eq.master) then
        open(unit=72,file='numstart.dat',status='old')
        read(72,*) irstart
        close(72)
      endif
      call drvAllreduceMaxInt(irstart,readstep)
      irstart=readstep
      lstep=irstart ! in case restart files have no fields

      fname1='geombc.dat'
      fname1= trim(fname1)  // cname2(myrank+1)

      itmp=1
      if (irstart .gt. 0) itmp = int(log10(float(irstart)))+1
      write (fmt1,"('(''restart.'',i',i1,',1x)')") itmp
      write (fnamer,fmt1) irstart
      fnamer = trim(fnamer) // cname2(myrank+1)
c
c.... open input files
c.... input the geometry parameters
c
      numparts = numpe !This is the common settings. Beware if you try to compute several parts per process

      itwo=2
      ione=1
      ieleven=11
      itmp = int(log10(float(myrank+1)))+1

      iotime = TMRC()
      if( input_mode .eq. -1 ) then
        call streamio_setup_read(fhandle, geomRestartStream)
      else if( input_mode .eq. 0 ) then
        call posixio_setup(fhandle, c_char_'r')
      else if( input_mode .ge. 1 ) then
        call syncio_setup_read(nsynciofiles, fhandle)
      end if
      call phio_constructName(fhandle, 
     &        c_char_'geombc' // char(0), fname1)
      call phastaio_setfile(GEOMBC_READ)
      call phio_openfile(fname1, fhandle);

      call phio_readheader(fhandle,c_char_'number of nodes' // char(0),
     & c_loc(numnp),ione, dataInt, iotype)

      call phio_readheader(fhandle,c_char_'number of modes' // char(0),
     & c_loc(nshg),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     &  c_char_'number of interior elements' // char(0),
     &  c_loc(numel),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     &  c_char_'number of boundary elements' // char(0),
     &  c_loc(numelb),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     &  c_char_'maximum number of element nodes' // char(0),
     &  c_loc(nen),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     &  c_char_'number of interior tpblocks' // char(0),
     &  c_loc(nelblk),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     & c_char_'number of boundary tpblocks' // char(0),
     & c_loc(nelblb),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     & c_char_'number of nodes with Dirichlet BCs' // char(0),
     & c_loc(numpbc),ione, dataInt, iotype)

      call phio_readheader(fhandle,
     & c_char_'number of shape functions' // char(0),
     & c_loc(ntopsh),ione, dataInt, iotype)
c
c.... calculate the maximum number of boundary element nodes
c     
      nenb = 0
      do i = 1, melCat
         if (nen .eq. nenCat(i,nsd)) nenb = max(nenCat(i,nsd-1), nenb)
      enddo
      if (myrank == master) then
         if (nenb .eq. 0) call error ('input   ','nen     ',nen)
      endif
c
c.... setup some useful constants
c
      I3nsd  = nsd / 3          ! nsd=3 integer flag
      E3nsd  = float(I3nsd)     ! nsd=3 real    flag
      if(matflg(1,1).lt.0) then
         nflow = nsd + 1
      else
         nflow = nsd + 2
      endif 
      ndof   = nsd + 2
      nsclr=impl(1)/100
      ndof=ndof+nsclr           ! number of sclr transport equations to solve
      
      ndofBC = ndof + I3nsd     ! dimension of BC array
      ndiBCB = 2                ! dimension of iBCB array
      ndBCB  = ndof + 1         ! dimension of BCB array
      nsymdf = (ndof*(ndof + 1)) / 2 ! symm. d.o.f.'s
c
c.... ----------------------> Communication tasks <--------------------
c
      if(numpe > 1) then
         call phio_readheader(fhandle,
     &    c_char_'size of ilwork array' // char(0),
     &    c_loc(nlwork),ione, dataInt, iotype)

         call phio_readheader(fhandle,
     &    c_char_'ilwork' //char(0),
     &    c_loc(nlwork),ione, dataInt, iotype)

         allocate( point2ilwork(nlwork) )
         allocate( ilworkread(nlwork) )
         call phio_readdatablock(fhandle, c_char_'ilwork' // char(0),
     &      c_loc(ilworkread), nlwork, dataInt , iotype)

         point2ilwork = ilworkread
         call ctypes (point2ilwork)

       if((usingPETSc.eq.1).or.(svLSFlag.eq.1)) then
         fncorpsize = nshg
         allocate(fncorp(fncorpsize))
         call gen_ncorp(fncorp, ilworkread, nlwork, fncorpsize)
!
! the  following code finds the global range of the owned nodes
!
         maxowned=0
         minowned=maxval(fncorp)
         do i = 1,nshg      
            if(fncorp(i).gt.0) then  ! don't consider remote copies
               maxowned=max(maxowned,fncorp(i))
               minowned=min(minowned,fncorp(i))
            endif
         enddo
!
!  end of global range code
!
         call commuInt(fncorp, point2ilwork, 1, 'out')
         ncorpsize = fncorpsize 
       endif
       if(svLSFlag.eq.1) then
         allocate(ltg(ncorpsize))
         ltg(:)=fncorp(:)
       endif
      else  !serial
       if((usingPETSc.eq.1)) then
         fncorpsize = nshg
         allocate(fncorp(fncorpsize))
         maxowned=nshg
         minowned=1
         iownnodes=nshg
           do i =1,nshg
             fncorp(i)=i
           enddo
       endif
       if((svLSFlag.eq.1)) then
           allocate(ltg(nshg))
           do i =1,nshg
             ltg(i)=i
           enddo
       endif
       nlwork=1
       allocate( point2ilwork(1))
       nshg0 = nshg
      endif


      itwo=2

      call phio_readheader(fhandle,
     & c_char_'co-ordinates' // char(0),
     & c_loc(intfromfile),itwo, dataDbl, iotype)
      numnp=intfromfile(1)
      allocate( point2x(numnp,nsd) )
      allocate( xread(numnp,nsd) )
      ixsiz=numnp*nsd
      call phio_readdatablock(fhandle,
     & c_char_'co-ordinates' // char(0),
     & c_loc(xread),ixsiz, dataDbl, iotype)
      point2x = xread

      if (iMoveTopWall.eq.1) then
        open (unit=123,file="TopWallWrong.dat",status="old")
        read(123,*) nx
        allocate(tpW(nx,2))
        do i=1,nx
          read(123,*) (tpW(i,j),j=1,2)
        enddo
        close(123)
        open (unit=124,file="TopWallWrong2.dat",status="old")
        read(124,*) nx
        allocate(tpW2(nx,2))
        do i=1,nx
           read(124,*) (tpW2(i,j),j=1,2)
        enddo
        close(124)
        open (unit=125,file="TopWallCorrect.dat",status="old")
        read(125,*) nx
        allocate(tpC(nx,2))
        do i=1,nx
           read(125,*) (tpC(i,j),j=1,2)
        enddo
        close(125)
        open (unit=126,file="TopWallCorrect2.dat",status="old")
        read(126,*) nx
        allocate(tpC2(nx,2))
        do i=1,nx
           read(126,*) (tpC2(i,j),j=1,2)
        enddo
        close(126)
        ic = 0
        do n=1,nshg
           do j=1,nx 
             if (abs(point2x(n,1)-tpW(j,1)).lt.1.0d-5.and.
     &         abs(point2x(n,2)-tpW(j,2)).lt.1.0d-5) then
               point2x(n,2) = tpC(j,2)
               ic = ic+1
               exit
             endif
             if (abs(point2x(n,1)-tpW2(j,1)).lt.1.0d-5.and.
     &         abs(point2x(n,2)-tpW2(j,2)).lt.1.0d-5) then
               point2x(n,2) = tpC2(j,2)
               ic = ic+1
               exit
             endif
           enddo
        enddo
        deallocate(tpW,tpW2,tpC,tpC2)
        !write(*,*) 'Moved ',ic,' nodes'
      endif



c..............................for Duct
      if(istretchOutlet.eq.1)then
         
c...geometry6
        if(iDuctgeometryType .eq. 6) then
          xmaxn = 1.276
          xmaxo = 0.848
          xmin  = 0.42
c...geometry8
        elseif(iDuctgeometryType .eq. 8)then
          xmaxn=1.6*4.5*0.0254+0.85*1.5
          xmaxo=1.6*4.5*0.0254+0.85*1.0
          xmin =1.6*4.5*0.0254+0.85*0.5
        endif
c...
        alpha=(xmaxn-xmaxo)/(xmaxo-xmin)**2
        where (point2x(:,1) .ge. xmin)
c..... N=# of current elements from .42 to exit(~40)
c..... (x_mx-x_mn)/N=.025
c..... alpha=3    3*.025=.075
           point2x(:,1)=point2x(:,1)+
     &     alpha*(point2x(:,1)-xmin)**2
c..... ftn to stretch x at exit
        endwhere
      endif

c
c.... read in and block out the connectivity
c
      if( input_mode .eq. -1 ) then
        call genblk 
      else if( input_mode .eq. 0 ) then
        call genblkPosix 
      else if( input_mode .ge. 1 ) then
        call genblkSyncIO 
      end if
c
c.... read the boundary condition mapping array
c
      ione=1
      call phio_readheader(fhandle,
     & c_char_'bc mapping array' // char(0),
     & c_loc(nshg),ione, dataInt, iotype)

      allocate( nBC(nshg) )

      allocate( nBCread(nshg) )

      call phio_readdatablock(fhandle,
     & c_char_'bc mapping array' // char(0),
     & c_loc(nBCread), nshg, dataInt, iotype)

      nBC=nBCread
c
c.... read the temporary iBC array
c
      ione=1
      call phio_readheader(fhandle,
     & c_char_'bc codes array' // char(0),
     & c_loc(numpbc),ione, dataInt, iotype)

      if ( numpbc > 0 ) then
        allocate( iBCtmp(numpbc) )
        allocate( iBCtmpread(numpbc) )
      else
        allocate( iBCtmp(1) )
        allocate( iBCtmpread(1) )
      endif
      call phio_readdatablock(fhandle,
     & c_char_'bc codes array' // char(0),
     & c_loc(iBCtmpread), numpbc, dataInt, iotype)

      if ( numpbc > 0 ) then
         iBCtmp=iBCtmpread
      else  ! sometimes a partition has no BC's
         deallocate( iBCtmpread)
         iBCtmp=0
      endif
c
c.... read boundary condition data
c
      ione=1
      call phio_readheader(fhandle,
     & c_char_'boundary condition array' // char(0),
     & c_loc(intfromfile),ione, dataDbl, iotype)

      if ( numpbc > 0 ) then
         allocate( BCinp(numpbc,ndof+7) )
         nsecondrank=intfromfile(1)/numpbc
         allocate( BCinpread(numpbc,nsecondrank) )
         iBCinpsiz=intfromfile(1)
      else
         allocate( BCinp(1,ndof+7) )
         allocate( BCinpread(0,0) ) !dummy
         iBCinpsiz=intfromfile(1)
      endif

      call phio_readdatablock(fhandle,
     & c_char_'boundary condition array' // char(0),
     & c_loc(BCinpread), iBCinpsiz, dataDbl, iotype)

      if ( numpbc > 0 ) then
         BCinp(:,1:(ndof+7))=BCinpread(:,1:(ndof+7))
      else  ! sometimes a partition has no BC's
         deallocate(BCinpread)
         BCinp=0
      endif
c
c.... read periodic boundary conditions
c
      ione=1
      call phio_readheader(fhandle,
     & c_char_'periodic masters array' // char(0),
     & c_loc(nshg), ione, dataInt, iotype)
      allocate( point2iper(nshg) )
      allocate( iperread(nshg) )
      call phio_readdatablock(fhandle,
     & c_char_'periodic masters array' // char(0),
     & c_loc(iperread), nshg, dataInt, iotype)
      point2iper=iperread
c
c.... generate the boundary element blocks
c
      if( input_mode .eq. -1 ) then
        if(numelb.gt.0) then
           call genbkb 
        else
           nelblb=0
        endif
      else if( input_mode .eq. 0 ) then
        call genbkbPosix 
      else if( input_mode .ge. 1 ) then
        call genbkbSyncIO 
      end if
c
c  Read in the nsons and ifath arrays if needed
c
c  There is a fundamental shift in the meaning of ifath based on whether
c  there exist homogenous directions in the flow.  
c
c  HOMOGENOUS DIRECTIONS EXIST:  Here nfath is the number of inhomogenous
c  points in the TOTAL mesh.  That is to say that each partition keeps a 
c  link to  ALL inhomogenous points.  This link is furthermore not to the
c  sms numbering but to the original structured grid numbering.  These 
c  inhomogenous points are thought of as fathers, with their sons being all
c  the points in the homogenous directions that have this father's 
c  inhomogeneity.  The array ifath takes as an arguement the sms numbering
c  and returns as a result the father.
c
c  In this case nsons is the number of sons that each father has and ifath
c  is an array which tells the 
c
c  NO HOMOGENOUS DIRECTIONS.  In this case the mesh would grow to rapidly
c  if we followed the above strategy since every partition would index its
c  points to the ENTIRE mesh.  Furthermore, there would never be a need
c  to average to a node off processor since there is no spatial averaging.
c  Therefore, to properly account for this case we must recognize it and
c  inerrupt certain actions (i.e. assembly of the average across partitions).
c  This case is easily identified by noting that maxval(nsons) =1 (i.e. no
c  father has any sons).  Reiterating to be clear, in this case ifath does
c  not point to a global numbering but instead just points to itself.
c
      nfath=1  ! some architectures choke on a zero or undeclared
                 ! dimension variable.  This sets it to a safe, small value.
      if(((iLES .lt. 20) .and. (iLES.gt.0))
     &                   .or. (itwmod.gt.0)
     &                   .or. (ispanAvg.eq.1) ) then ! don't forget same
                                                    ! conditional in proces.f
                                                    ! needed for alloc
         ione=1
         if(nohomog.gt.0) then
            call phio_readheader(fhandle,
     &       c_char_'number of father-nodes' // char(0),
     &       c_loc(nfath), ione, dataInt, iotype)

            call phio_readheader(fhandle,
     &       c_char_'number of son-nodes for each father' // char(0),
     &       c_loc(nfath), ione, dataInt, iotype)

            allocate(nsonsread(nfath))

            call phio_readdatablock(fhandle,
     &       c_char_'number of son-nodes for each father' // char(0),
     &       c_loc(nsonsread),nfath, dataInt, iotype)
      
            nsonmax=maxval(nsonsread)
            nsonmin=minval(nsonsread)
            if (nsonmax.eq.nsonmin) then
                ndistsons = 1
                allocate(point2nsons(ndistsons))
                point2nsons = nsonmax
            else
                ndistsons = nfath
                allocate(point2nsons(ndistsons))
                point2nsons = nsonsread
            endif
            deallocate(nsonsread)

            call phio_readheader(fhandle,
     &       c_char_'keyword ifath' // char(0),
     &       c_loc(nshg), ione, dataInt, iotype);

            allocate (point2ifath(nshg))

            call phio_readdatablock(fhandle,
     &       c_char_'keyword ifath' // char(0),
     &       c_loc(point2ifath), nshg, dataInt, iotype)
     
            nsonmax=maxval(point2nsons)

            if(myrank.eq.master) write(*,*) 'Number of fathers is: ',nfath
            if(myrank.eq.master) write(*,*) 'Number of sons is: ',nsonmax
            inquire(file="dynSmagY.dat",exist=exlog)
            if(exlog) then
              if(myrank.eq.master) write(*,*) 
     &                        "Setting ifath from dynSmagY.dat"
              open (unit=123,file="dynSmagY.dat",status="old")
              read(123,*) ny
              nfath = ny
              ndistsons = 1
              point2nsons = 4800
              nsonmax=maxval(point2nsons)
              allocate(ypoints(ny))
              do i=1,ny
                read(123,*) ypoints(i)
              enddo
              do i=1,nshg
                do j=1,ny
                  if (abs(point2x(i,2)-ypoints(j)).lt.1.0e-4) then
                    point2ifath(i) = j            
                  endif
                enddo
              enddo
              deallocate(ypoints)
              close(123)
              if(myrank.eq.master) write(*,*) 'Number of fathers is: ',nfath
              if(myrank.eq.master) write(*,*) 'Number of sons is: ',nsonmax
            else
               if(myrank.eq.master) 
     &               write(*,*) 'Did not read file dynSmagY.dat'
            endif
c            do i=1,nshg
c               write(*,*) point2x(i,1),point2x(i,2),point2x(i,3),point2ifath(i)
c            enddo
         else  ! this is the case where there is no homogeneity
               ! therefore ever node is a father (too itself).  sonfath
               ! (a routine in NSpre) will set this up but this gives
               ! you an option to avoid that.
            nfath=nshg
            ndistsons = nshg
            allocate (point2nsons(ndistsons))
            point2nsons=1
            allocate (point2ifath(nshg))
            do i=1,nshg
               point2ifath(i)=i
            enddo
            nsonmax=1
         endif
      else
         ndistsons = 1
         allocate (point2nsons(ndistsons))
         allocate (point2ifath(1))
      endif

      call phio_closefile(fhandle);
      iotime = TMRC() - iotime
      if (myrank.eq.master) then
        write(*,*) 'time to read geombc (seconds)', iotime
      endif

c.... Read restart files
      iotime = TMRC()
      if( input_mode .eq. -1 ) then
        call streamio_setup_read(fhandle, geomRestartStream)
      else if( input_mode .eq. 0 ) then
        call posixio_setup(fhandle, c_char_'r')
      else if( input_mode .ge. 1 ) then
        call syncio_setup_read(nsynciofiles, fhandle)
      end if
      call phio_constructName(fhandle,
     &        c_char_'restart' // char(0), fnamer)
      call phstr_appendInt(fnamer, irstart)
      call phstr_appendStr(fnamer, c_char_'.'//c_null_char)
      call phastaio_setfile(RESTART_READ)
      call phio_openfile(fnamer, fhandle);

      ithree=3

      itmp = int(log10(float(myrank+1)))+1

      intfromfile=0
      call phio_readheader(fhandle,
     & c_char_'solution' // char(0), 
     & c_loc(intfromfile), ithree, dataInt, iotype)
c
c.... read the values of primitive variables into q
c
      allocate( qold(nshg,ndof) )
      if(intfromfile(1).ne.0) then
         nshg2=intfromfile(1)
         ndof2=intfromfile(2)
         lstep=intfromfile(3)
         if(ndof2.ne.ndof) then

         endif
        if (nshg2 .ne. nshg)
     &        call error ('restar  ', 'nshg   ', nshg)
         allocate( qread(nshg,ndof2) )
         iqsiz=nshg*ndof2
         call phio_readdatablock(fhandle,
     &    c_char_'solution' // char(0),
     &    c_loc(qread),iqsiz, dataDbl,iotype)
         qold(:,1:ndof)=qread(:,1:ndof)
         deallocate(qread)
      else
         if (myrank.eq.master) then
            if (matflg(1,1).eq.0) then ! compressible
               warning='Solution is set to zero (with p and T to one)'
            else
               warning='Solution is set to zero'
            endif
            write(*,*) warning
         endif
         qold=zero
         if (matflg(1,1).eq.0) then ! compressible
            qold(:,1)=one ! avoid zero pressure
            qold(:,nflow)=one ! avoid zero temperature
         endif
      endif

      intfromfile=0
      call phio_readheader(fhandle,
     & c_char_'time derivative of solution' // char(0),
     & c_loc(intfromfile), ithree, dataInt, iotype)
      allocate( acold(nshg,ndof) )
      if(intfromfile(1).ne.0) then
         nshg2=intfromfile(1)
         ndof2=intfromfile(2)
         lstep=intfromfile(3)

         if (nshg2 .ne. nshg)
     &        call error ('restar  ', 'nshg   ', nshg)
         allocate( acread(nshg,ndof2) )
         acread=zero
         iacsiz=nshg*ndof2
         call phio_readdatablock(fhandle,
     &    c_char_'time derivative of solution' // char(0),
     &    c_loc(acread), iacsiz, dataDbl,iotype)
         acold(:,1:ndof)=acread(:,1:ndof)
         deallocate(acread)
      else
         if (myrank.eq.master) then
            warning='Time derivative of solution is set to zero (SAFE)'
            write(*,*) warning
         endif
         acold=zero
      endif


cc
cc.... Read the STG  fields if lstep not equal to iSTGStart
cc
      if (lstep.ne.iSTGStart.and.iSTG.eq.1) then
cc....   First read the STG random variables
         intfromfile=0
         call phio_readheader(fhandle,
     &   c_char_'STG_rnd' // char(0), 
     &   c_loc(intfromfile), ithree, dataInt, iotype)

         if(intfromfile(1).ne.0) then
           nKWave=intfromfile(1)
           irandVars=intfromfile(2)
           if (myrank.eq.master) then
             write(*,*) 'nKWave for STG inflow is set to : ',nKWave
             write(*,*) 'Number of random variables is : ', irandVars
           endif
           allocate(STGread(nKWave,irandVars))
           allocate(STGrnd(nKWave,irandVars))
           iSTGsiz=nKWave*irandVars
           call phio_readdatablock(fhandle,
     &       c_char_'STG_rnd' // char(0),
     &       c_loc(STGread),iSTGsiz, dataDbl,iotype)
      
           STGrnd=STGread
           deallocate(STGread)

         else
           if (myrank.eq.master) then
             warning='The random variables for the STG 
     &                   inflow were not read and set to zero'
             write(*,*) warning
           endif
           STGrnd=zero
         endif

cc....   Now read the BC array. This was a quick fix to problem with inflow BC value set
cc....   in SimModeler and present in the geombc files
         intfromfile=0
         call phio_readheader(fhandle,
     &   c_char_'BCs' // char(0), 
     &   c_loc(intfromfile), ithree, dataInt, iotype)

         if(intfromfile(1).ne.0) then
           allocate(BCrestartRead(nshg,ndofBC))
           allocate(BCrestart(nshg,ndofBC))
           iSTGsiz=nshg*ndofBC
           call phio_readdatablock(fhandle,
     &       c_char_'BCs' // char(0),
     &       c_loc(BCrestartRead),iSTGsiz, dataDbl,iotype)
      
           BCrestart=BCrestartRead
           deallocate(BCrestartRead)

         else
           if (myrank.eq.master) then
             warning='BCs from restart were not read'
             write(*,*) warning
           endif
           !BCrestart=zero
         endif

      endif  ! end of the STG fields

cc
cc.... Read the velbar array if want spanwise average
cc
      if (ispanAvg.eq.1) then
         if (myrank.eq.master) then 
           itmp = 1
           if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
           write (fmtv,"('(''velbar.'',i',i1,',1x)')") itmp
           write (fnamev,fmtv) lstep
           fnamev = trim(fnamev) // cname(myrank+1)
           inquire(file=fnamev,exist=exlog)
           if (exlog) then
             allocate(velbar(nfath,ndof))
             open(unit=123,file=fnamev,status="old")
             do i=1,nfath
                read(123,*) (velbar(i,j),j=1,ndof)
             enddo
             close(123)
             write(*,*) 'Read velbar from file ',fnamev
           else ! did not find velbar for current time step
             write(*,*) 'Velbar not read from file, setting to zero'
             allocate(velbar(nfath,ndof),STAT=IERR1)
             if(IERR1.gt.0) write(*,*) 
     &                 'Not enough space to allocate velbar'
             velbar = zero
           endif
         endif
         if (numpe .gt. 1) call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      endif ! end of ispanAvg for velbar

cc
cc.... Read the stsbar array if want spanwise average and stats
cc
      if (ispanAvg.eq.1) then
         if (myrank.eq.master) then 
           itmp = 1
           if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
           write (fmts,"('(''stsbar.'',i',i1,',1x)')") itmp
           write (fnames,fmts) lstep
           fnames = trim(fnames) // cname(myrank+1)
           inquire(file=fnames,exist=exlog)
           if (exlog) then
             allocate(stsBar(nfath,iConsStressSz))
             open(unit=123,file=fnames,status="old")
             do i=1,nfath
                read(123,*) (stsBar(i,j),j=1,iConsStressSz)
             enddo
             close(123)
             write(*,*) 'Read stsbar from file ',fnames
           else ! did not find velbar for current time step
             write(*,*) 'stsbar not read from file, setting to zero'
             allocate(stsBar(nfath,iConsStressSz),STAT=IERR1)
             if(IERR1.gt.0) write(*,*)
     &                    'Not enough space to allocate stsBar'
             stsBar = zero
           endif
         endif
         if (numpe .gt. 1) call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      endif ! end of ispanAvg for stsbar
cc
cc.... Read the stsbarKeq array if want spanwise average and K eq stats
cc
      if (ispanAvg.eq.1.and.iKeq.eq.1) then
         if (myrank.eq.master) then 
           itmp = 1
           if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
           write (fmtk,"('(''stsbarKeq.'',i',i1,',1x)')") itmp
           write (fnamek,fmtk) lstep
           fnamek = trim(fnamek) // cname(myrank+1)
           inquire(file=fnamek,exist=exlog)
           if (exlog) then
             allocate(stsBarKeq(nfath,10))
             open(unit=123,file=fnamek,status="old")
             do i=1,nfath
                read(123,*) (stsBarKeq(i,j),j=1,10)
             enddo
             close(123)
             write(*,*) 'Read stsbarKeq from file ',fnamek
           else ! did not find stsbarKeq for current time step
             write(*,*) 'stsBarKeq not read from file, setting to zero'
             allocate(stsBarKeq(nfath,10),STAT=IERR1)
             if(IERR1.gt.0) write(*,*) 
     &                 'Not enough space to allocate stsBarKeq'
             stsBarKeq = zero
           endif
         endif
         if (numpe .gt. 1) call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      endif ! end of if for ispanAvg and iKeq for K eq terms
cc
cc.... Read the IDDESbar array if want spanwise average of IDDES functions
cc
      if (ispanAvg.eq.1.and.ispanIDDES.eq.1) then
         if (myrank.eq.master) then 
           itmp = 1
           if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
           write (fmtI,"('(''IDDESbar.'',i',i1,',1x)')") itmp
           write (fnameI,fmtI) lstep
           fnameI = trim(fnameI) // cname(myrank+1)
           inquire(file=fnameI,exist=exlog)
           if (exlog) then
             allocate(IDDESbar(nfath,nfun))
             open(unit=123,file=fnameI,status="old")
             do i=1,nfath
                read(123,*) (IDDESbar(i,j),j=1,nfun)
             enddo
             close(123)
             write(*,*) 'Read IDDESbar from file ',fnameI
           else ! did not find IDDESbar for current time step
             write(*,*) 'IDDESbar not read from file, setting to zero'
             allocate(IDDESbar(nfath,nfun),STAT=IERR1)
             if(IERR1.gt.0) write(*,*) 
     &                 'Not enough space to allocate IDDESbar'
             IDDESbar = zero
           endif
         endif
         if (numpe .gt. 1) call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      endif ! end of if for ispanAvg and ispanIDDES for IDDES functions
cc
cc.... Read the cdelsq array if want running spanwise average and doing LES
cc
      if (iLES.gt.0.and.irunTave.eq.1) then
         intfromfile=0
         call phio_readheader(fhandle,
     &   c_char_'cdelsq' // char(0), 
     &   c_loc(intfromfile), ithree, dataInt, iotype)

         if(intfromfile(1).ne.0) then
           nshg2=intfromfile(1)
           if (nshg.ne.nshg2)
     &          call error ('restar  ', 'cdelsq   ', cdelsq)
           allocate(cdelsqread(nshg2,3))
           allocate(cdelsq(nshg2,3))
           isiz=nshg2*3
           call phio_readdatablock(fhandle,
     &       c_char_'cdelsq' // char(0),
     &       c_loc(cdelsqread),isiz, dataDbl,iotype)
           cdelsq=cdelsqread
           deallocate(cdelsqread)
           if (myrank.eq.master) 
     &            write(*,*) 'cdelsq field found in restart file'
         else
           if (myrank.eq.master) then
             warning='cdelsq not read from restart and set to zero'
             write(*,*) warning
           endif
           allocate(cdelsq(nshg,3),STAT=IERR2)
           if(IERR2.gt.0)write(*,*)'Not enough space to allocate cdelsq'
           cdelsq=zero
         endif
      else if (iLES.gt.0) then
         allocate(cdelsq(nshg,3))
         cdelsq=zero
      else if (iLES.eq.0) then
         allocate(cdelsq(1,1))
      endif ! end of iLES and irunTave for cdelsq
cc
cc.... Read the lesnut array if want running spanwise average and doing LES
cc
      if (iLES.gt.0.and.irunTave.eq.1) then
         intfromfile=0
         call phio_readheader(fhandle,
     &   c_char_'lesnut' // char(0), 
     &   c_loc(intfromfile), ithree, dataInt, iotype)

         if(intfromfile(1).ne.0) then
           numel2=intfromfile(1)
           if (numel.ne.numel2)
     &          call error ('restar  ', 'lesnut   ', lesnut)
           allocate(lesnutread(numel2,2))
           allocate(lesnut(numel2,2))
           isiz=numel2*2
           call phio_readdatablock(fhandle,
     &       c_char_'lesnut' // char(0),
     &       c_loc(lesnutread),isiz, dataDbl,iotype)
           lesnut=lesnutread
           deallocate(lesnutread)
           if (myrank.eq.master) 
     &            write(*,*) 'lesnut field found in restart file'
         else
           if (myrank.eq.master) then
             warning='lesnut not read from restart and set to zero'
             write(*,*) warning
           endif
           allocate(lesnut(numel,2),STAT=IERR2)
           if(IERR2.gt.0)write(*,*)'Not enough space to allocate lesnut'
           lesnut=zero
         endif
      else if (iLES.gt.0) then
         allocate(lesnut(numel,2))
         lesnut=zero
      else if (iLES.eq.0) then
         allocate(lesnut(1,1))
      endif ! end of iLES and irunTave for lesnut

cc
cc.... read the header and check it against the run data
cc
      if (ideformwall.eq.1) then

          intfromfile=0
          call phio_readheader(fhandle,
     &     c_char_'displacement' // char(0),
     &     c_loc(intfromfile), ithree, dataInt, iotype)

         nshg2=intfromfile(1)
         ndisp=intfromfile(2)
         lstep=intfromfile(3)
         if(ndisp.ne.nsd) then
            warning='WARNING ndisp not equal nsd'
            write(*,*) warning , ndisp
         endif
         if (nshg2 .ne. nshg) 
     &        call error ('restar  ', 'nshg   ', nshg)
c
c.... read the values of primitive variables into uold
c
         allocate( uold(nshg,nsd) )
         allocate( uread(nshg,nsd) )
         
         iusiz=nshg*nsd

         call phio_readdatablock(fhandle,
     &    c_char_'displacement' // char(0),
     &    c_loc(uread), iusiz, dataDbl, iotype)

         uold(:,1:nsd)=uread(:,1:nsd)
       else
         allocate( uold(nshg,nsd) )
         uold(:,1:nsd) = zero
       endif
c
c.... close c-binary files
c
      call phio_closefile(fhandle)
      iotime = TMRC() - iotime
      if (myrank.eq.master) then
        write(*,*) 'time to read restart (seconds)', iotime
      endif

      deallocate(xread)
      if ( numpbc > 0 )  then
         deallocate(bcinpread)
         deallocate(ibctmpread)
      endif
      deallocate(iperread)
      if(numpe.gt.1)
     &     deallocate(ilworkread)
      deallocate(nbcread)

      return
 994  call error ('input   ','opening ', igeomBAK)
 995  call error ('input   ','opening ', igeomBAK)
 997  call error ('input   ','end file', igeomBAK)
 998  call error ('input   ','end file', igeomBAK)
      end
c
c No longer called but kept around in case....
c
      subroutine genpzero(iBC)

      use pointer_data
      include "common.h"
      integer iBC(nshg)
c
c....  check to see if any of the nodes have a dirichlet pressure
c
      pzero=1
      if (any(btest(iBC,2))) pzero=0  
      do iblk = 1, nelblb
         npro = lcblkb(1,iblk+1)-lcblkb(1,iblk)
         do i=1, npro
            iBCB1=miBCB(iblk)%p(i,1)
c     
c.... check to see if any of the nodes have a Neumann pressure 
c     but not periodic (note that 
c     
            if(btest(iBCB1,1)) pzero=0
         enddo
c     
c.... share results with other processors
c     
         pzl=pzero
         if (numpe .gt. 1)
     &        call MPI_ALLREDUCE (pzl, pzero, 1,
     &        MPI_DOUBLE_PRECISION,MPI_MIN, MPI_COMM_WORLD,ierr)
      enddo
      return
      end
