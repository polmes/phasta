      subroutine getvel (y,ilwork, iBC,
     &                   nsons, ifath)

      use spanStats

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

c
      dimension nsons(ndistsons),        rinvsons(ndistsons),
     &          velo(nshg,ndof),    velf(nfath),
     &          velft(nfath), 
     &          y(nshg,ndof), 
     &          ifath(nshg),
     &          ilwork(nlwork),        iBC(nshg)
      real*8 tmp, den, tfact 
      integer periodicity, sumper
      
      peroidicity = 0
      den = max(1,lstep-istartSpanAvg)
      tfact = one/den

c
c  for now keep the compressible numbering in velbar
c
      velo = zero
      velo(:,1)=y(:,4)
      velo(:,2)=y(:,1)      
      velo(:,3)=y(:,2)
      velo(:,4)=y(:,3)
      if(nflow.eq.5) velo(:,5)=y(:,5)
      if (nsclr.gt.0) then
         if (iRANS.eq.-1) velo(:,6) = y(:,6)
         if (iRANS.eq.-5) velo(:,6:7) = y(:,6:7)
      endif
      nsonmax=maxval(nsons)
      if ((nsonmax.eq.1)) then  ! we are doing local clipping -no homog dir
         !velft=velo
      else
c     
c     zero on processor periodic nodes so that they will not be added twice
c    
         do i=1, nshg
            if(btest(iBC(i),10)) then 
              periodicity = 1
              exit
            endif
         enddo
         call drvAllreduceSumInt(periodicity,sumper)
!         call MPI_BCAST(periodicity,1,MPI_INTEGER,master,
!     &               MPI_COMM_WORLD,ierr)
         if (sumper.gt.0) periodicity = 1

         if (periodicity.eq.1.and.nohomog.eq.1) then
             rinvsons = one/(nsons-one)   ! division is expensive
         else
             rinvsons = one/nsons   ! division is expensive
         endif

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
         if (nsclr.gt.0) then  
            if (iRANS.eq.-1) then
               where(btest(iBC,10).or.btest(iBC,12))
                  velo(:,6)=zero
               endwhere
            endif
            if (iRANS.eq.-5) then
               where(btest(iBC,10).or.btest(iBC,12))
                  velo(:,6)=zero
                  velo(:,7)=zero
               endwhere
            endif
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

c     
c     accumulate sum of sons to the fathers
c    
         do n=1,ndof 
            velf = zero
            do i = 1,nshg
               ifathi=ifath(i)
               velf(ifathi) = velf(ifathi) 
     &                             + velo(i,n)            
            enddo
         
c     
c     Now  the true fathers and serrogates combine results and update
c     each other.
c     
            if(numpe .gt. 1) then
               call drvAllreduceDP(velf, velft,nfath)
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
            if (ndistsons.eq.nfath) then 
              velft(:) = velft(:) * rinvsons(:) !  / nsons(:)
            else if (ndistsons.eq.1) then
              tmp =  rinvsons(1)
              velft = velft * tmp
            endif
          
            if (myrank.eq.master) then 
               velbar(:,n)=tfact*velft+(one-tfact)*velbar(:,n)
            endif

         enddo
            

      endif  ! end of homog direction averaging

      
      
      return
      end subroutine getvel

c ===============================================================================
c    getvel2
c ===============================================================================
      subroutine getvel2 (y, ilwork, iBC, ifath)

      use spanStats

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      dimension velo(nshg,ndof),  velf(locnfath),
     &          y(nshg,ndof), 
     &          ifath(nshg),
     &          ilwork(nlwork), iBC(nshg)
      real*8 tmp, den, tfact 
      real*4 tmpcount(nshg), locnsons(locnfath), c(locnfath)

      den = max(1,lstep-istartSpanAvg)
      tfact = one/den
      tmpcount = one

c
c  for now keep the compressible numbering in velbar
c
      velo = zero
      velo(:,1)=y(:,4)
      velo(:,2)=y(:,1)      
      velo(:,3)=y(:,2)
      velo(:,4)=y(:,3)
      if(nflow.eq.5) velo(:,5)=y(:,5)
      if (nsclr.gt.0) then
         if (iRANS.eq.-1) velo(:,6) = y(:,6)
         if (iRANS.eq.-5) velo(:,6:7) = y(:,6:7)
      endif

c     
c     zero on processor periodic nodes so that they will not be added twice
c    
      where(btest(iBC,10).or.btest(iBC,12))
         velo(:,1)=zero
         velo(:,2)=zero
         velo(:,3)=zero
         velo(:,4)=zero
         tmpcount(:) = zero
      endwhere      
      if(nflow.eq.5) then
         where(btest(iBC,10).or.btest(iBC,12))
             velo(:,5)=zero
         endwhere
      endif
      if (nsclr.gt.0) then  
         if (iRANS.eq.-1) then
            where(btest(iBC,10).or.btest(iBC,12))
               velo(:,6)=zero
            endwhere
         endif
         if (iRANS.eq.-5) then
            where(btest(iBC,10).or.btest(iBC,12))
               velo(:,6)=zero
               velo(:,7)=zero
            endwhere
         endif
      endif     

c     zero the nodes that are "solved" on the other processors  
      if (numpe.gt.1) then
         numtask = ilwork(1)
         itkbeg = 1
         do itask = 1, numtask
            iacc   = ilwork (itkbeg + 2)
            numseg = ilwork (itkbeg + 4)
            if (iacc .eq. 0) then
               do is = 1,numseg
                  isgbeg = ilwork (itkbeg + 3 + 2*is)
                  lenseg = ilwork (itkbeg + 4 + 2*is)
                  isgend = isgbeg + lenseg - 1
                  velo(isgbeg:isgend,:) = zero 
                  tmpcount(isgbeg:isgend) = zero
               enddo
            endif
            itkbeg = itkbeg + 4 + 2*numseg
         enddo !itask
      endif  ! numpe.gt.1

c     
c     accumulate sum of sons to the fathers and accumulate onto velbar
c
      do n=1,ndof 
         velf = zero
         locnsons = zero
         c = zero
         do i = 1,nshg
            ifathi = ifath(i)
            ind = locifath(i)
            velf(ind) = velf(ind) + velo(i,n)
            if (velo(i,1).ne.zero.and.ifathi.ne.1.and.ifathi.ne.70) then
               c(ind) = c(ind)+1
            endif
            locnsons(ind) = locnsons(ind) + tmpcount(i)
            if (istep.eq.1) then
              velbar(ind,1) = ifathi
            endif
         enddo
!         velf(:) = velf(:) / locnsons(:)
!         where (locnsons.eq.zero)
!            velf(:) = zero
!         endwhere
         velbar(:,1+n)=tfact*velf+(one-tfact)*velbar(:,1+n)
      enddo


      end subroutine getvel2


