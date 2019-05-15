      subroutine getvel (y,ilwork, iBC,
     &                   nsons, ifath)

      use spanStats

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

c
      dimension nsons(ndistsons),        rinvsons(ndistsons),
     &          velo(nshg,nflow),    !velf(nfath,nflow),
     &          velft(nfath,nflow), 
     &          y(nshg,ndof), 
     &          ifath(nshg),
     &          ilwork(nlwork),        iBC(nshg)
      real*8 tmp, den, tfact 
      integer periodicity, minifath, maxifath, ifathrng,
     &        ifathi, ind, c
      logical, allocatable, dimension (:) :: isonpart 
      peroidicity = 0

c
c  for now keep the compressible numbering in velbar
c
      velo(:,1)=y(:,4)
      velo(:,2)=y(:,1)      
      velo(:,3)=y(:,2)
      velo(:,4)=y(:,3)
      if(nflow.eq.5) velo(:,5)=y(:,5)
      nsonmax=maxval(nsons)
      if ((nsonmax.eq.1)) then  ! we are doing local clipping -no homog dir
         velft=velo
      else
c     
c     zero on processor periodic nodes so that they will not be added twice
c    
         if(myrank.eq.master) then
            do i=1, nshg
               if(btest(iBC(i),10)) then 
                 periodicity = 1
                 exit
               endif
            enddo
         endif
         call MPI_BCAST(periodicity,1,MPI_INTEGER,master,
     &               MPI_COMM_WORLD,ierr)

         where(btest(iBC,10).or.btest(iBC,12))
            velo(:,1)=zero
            velo(:,2)=zero
            velo(:,3)=zero
            velo(:,4)=zero
         endwhere      
         if(nflow.eq.5) then
            where(btest(iBC,10).or.btest(iBC,12))
                velo(:,5)=zero
            endwhere
         endif         
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
                     velo(isgbeg:isgend,:) = zero
                  enddo
               endif
               
               itkbeg = itkbeg + 4 + 2*numseg
               
            enddo !itask
            
         endif  ! numpe.gt.1

         if (.not.allocated(velf)) then
            minifath = minval(ifath)
            maxifath = maxval(ifath)
            ifathrng = maxifath-minifath+1
            allocate(isonpart(ifathrng))
            isonpart = .false.
            locnfath = 0
            do ii=1,nshg
               ifathi=ifath(ii)
               do jj=minifath,maxifath
                  if (ifathi.eq.jj) then
                     ind = jj-minifath+1
                     if (.not.isonpart(ind)) then
                        isonpart(ind) = .true.
                        locnfath = locnfath+1
                     endif
                  endif
               enddo
            enddo
            allocate(velf(locnfath,nflow))
            allocate(locifath(locnfath))
            if (numpe.gt.1) then
                call drvAllreduceSumInt(locnfath,stacksz)
                if (.not.allocated(rcounts)) allocate(rcounts(numpe))
                if (.not.allocated(displs)) allocate(displs(numpe))
                call MPI_GATHER(locnfath,1,MPI_INT,rcounts,1,
     &                      MPI_INT,master,
     &                      MPI_COMM_WORLD,ierr)
                displs(1) = 0
                do n=2,numpe
                   displs(n) = displs(n-1)+rcounts(n-1)
                enddo
            endif
            c = 1
            isonpart = .false.
            do ii=1,nshg
               ifathi=ifath(ii)
               do jj=minifath,maxifath
                  if (ifathi.eq.jj) then
                     ind = jj-minifath+1
                     if (.not.isonpart(ind)) then
                        isonpart(ind) = .true.
                        locifath(c) = ifathi
                        c = c +1
                     endif
                  endif
               enddo
            enddo
            deallocate(isonpart)
         endif
         velf = zero
         velft = zero
c     
c     accumulate sum of sons to the fathers
c     
         do i = 1,nshg
            ifathi=ifath(i)
            do jj=1,locnfath
               if (ifathi.eq.locifath(jj)) then
                  velf(jj,1:nflow) = velf(jj,1:nflow)
     &                               + velo(i,1:nflow)
                  exit
               endif
            enddo
!            velf(ifathi,1:nflow) = velf(ifathi,1:nflow) 
!     &                             + velo(i,1:nflow)            
         enddo
         
c     
c     Now  the true fathers and serrogates combine results and update
c     each other.
c     
         if(numpe .gt. 1) then
            if (.not.allocated(ifathG)) then
               if (myrank.eq.master) allocate(ifathG(stacksz))
            endif
            if (.not.allocated(velftG)) then
               if (myrank.eq.master) allocate(velftG(stacksz,nflow))
            endif
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call MPI_GATHERV(locifath,locnfath,MPI_INT,ifathG,rcounts,
     &                       displs,MPI_INT,master,
     &                       MPI_COMM_WORLD,ierr)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            do i=1,nflow
               call MPI_BARRIER(MPI_COMM_WORLD,ierr)
               call MPI_GATHERV(velf(:,i),locnfath,MPI_DOUBLE_PRECISION,
     &                          velftG(:,i),rcounts,displs,
     &                          MPI_DOUBLE_PRECISION,master,
     &                          MPI_COMM_WORLD,ierr)
               call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            enddo
            if (myrank.eq.master) then
               do i=1,stacksz
                  ifathi = ifathG(i)
                  velft(ifathi,1:nflow) = velft(ifathi,1:nflow) 
     &                                    + velftG(i,1:nflow)
               enddo
            endif

            !call drvAllreduce(velf, velft,nfath*nflow)
         else
            velft=velf
         endif
c     
c     xvelf is the sum of the sons for each father on this processor
c
c     xvelft is the sum of the sons for each father on all processor combined
c     (the same as if we had not partitioned the mesh for each processor)
c     
c     divide by # of sons to get average father for this step
c
         if (periodicity.eq.1.and.nohomog.eq.1) then
           rinvsons = one/(nsons-one)   ! division is expensive
         else
           rinvsons = one/nsons   ! division is expensive
         endif
         if (myrank.eq.master) then
           if (ndistsons.eq.nfath) then 
             velft(:,1) = velft(:,1) * rinvsons(:) !  / nsons(:)
             velft(:,2) = velft(:,2) * rinvsons(:) !  / nsons(:)
             velft(:,3) = velft(:,3) * rinvsons(:) !  / nsons(:)
             velft(:,4) = velft(:,4) * rinvsons(:) !  / nsons(:)
             if(nflow.eq.5) velft(:,5) = velft(:,5) * rinvsons
           else if (ndistsons.eq.1) then
             tmp =  rinvsons(1)
             velft = velft * tmp
           endif
         endif
      endif  ! end of homog direction averaging

      if (myrank.eq.master) then 
         den = max(1,lstep-istartSpanAvg)
         tfact = one/den
         velbar(:,:)=tfact*velft(:,:)+(one-tfact)*velbar(:,:)
      endif

      
      
      return
      end


