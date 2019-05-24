      subroutine getStsBar(y, GradV, ilwork, nsons, ifath, iBC)

      use stats
      use spanStats
 
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      dimension nsons(ndistsons), rinvsons(ndistsons), tmpStats(nshg,iConsStressSz),
     &          ifath(nshg), ilwork(nlwork), iBC(nshg), tmpStatsf(nfath),
     &          tmpStatsft(nfath),
     &          y(nshg,ndof), tmpKeq(nshg,10), tmpKeqf(nfath),
     &          tmpKeqft(nfath), GradV(nshg,nsdsq)
      real*8 tmp(nshg),S11(nshg),S22(nshg),S33(nshg),
     &       S12(nshg),S13(nshg),S23(nshg), invtmp
      integer periodicity


      den = max(1,lstep-istartSpanAvg)
      tfact = one/den
      
c     Assign the conservative statistics to a temporary array
      if (iConsStress.eq.1) then 
         tmpStats(:,1:3) = stsVelInst(:,:)
         tmpStats(:,4:9) = stsVelSqInst(:,:) ! 11,22,33,12,23,31
      else
         tmpStats(:,1) = y(:,1)*y(:,1)
         tmpStats(:,2) = y(:,2)*y(:,2)
         tmpStats(:,3) = y(:,3)*y(:,3)
         tmpStats(:,4) = y(:,1)*y(:,2)
         tmpStats(:,5) = y(:,1)*y(:,3)
         tmpStats(:,6) = y(:,2)*y(:,3)
      endif
      if (iKeq.eq.1) then
         tmpKeq(:,1) = y(:,4)*y(:,1) ! p*u1
         tmpKeq(:,2) = y(:,4)*y(:,2) ! p*u2
         tmpKeq(:,3) = y(:,4)*y(:,3) ! p*u3
         tmp = y(:,1)**2 + y(:,2)**2 + y(:,3)**2
         tmpKeq(:,4) = y(:,1)*tmp ! u_i*u_j*u_j, i=1
         tmpKeq(:,5) = y(:,2)*tmp ! u_i*u_j*u_j, i=2
         tmpKeq(:,6) = y(:,3)*tmp ! u_i*u_j*u_j, i=3
         if (ipord.eq.1) then
           S11 = GradV(:,1)
           S22 = GradV(:,5)
           S33 = GradV(:,9)
           S12 = 0.50d0*(GradV(:,2)+GradV(:,4))
           S13 = 0.50d0*(GradV(:,3)+GradV(:,7))
           S23 = 0.50d0*(GradV(:,6)+GradV(:,8))
         endif
         tmpKeq(:,7) = y(:,1)*S11+y(:,2)*S12+y(:,3)*S13 ! u_j*S_ij, i=1
         tmpKeq(:,8) = y(:,1)*S12+y(:,2)*S22+y(:,3)*S23 ! u_j*S_ij, i=2
         tmpKeq(:,9) = y(:,1)*S13+y(:,2)*S23+y(:,3)*S33 ! u_j*S_ij, i=3
         tmpKeq(:,10) = S11**2 + S22**2 + S33**2 + 
     &                  two*(S12**2 + S13**2 + S23**2) ! S_ij*S_ij
      endif

c     Zero on processor periodic nodes so will not be added twice
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

      if (periodicity.eq.1.and.nohomog.eq.1) then
         rinvsons = one/(nsons-one)   ! division is expensive
      else
         rinvsons = one/nsons   ! division is expensive
      endif
 
      if (iConsStress.eq.1) then
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
        endwhere
      else
        where(btest(iBC,10).or.btest(iBC,12))
          tmpStats(:,1)=zero
          tmpStats(:,2)=zero
          tmpStats(:,3)=zero
          tmpStats(:,4)=zero
          tmpStats(:,5)=zero
          tmpStats(:,6)=zero
        endwhere
      endif

      if (iKeq.eq.1) then
        where(btest(iBC,10).or.btest(iBC,12))
          tmpKeq(:,1)=zero
          tmpKeq(:,2)=zero
          tmpKeq(:,3)=zero
          tmpKeq(:,4)=zero
          tmpKeq(:,5)=zero
          tmpKeq(:,6)=zero
          tmpKeq(:,7)=zero
          tmpKeq(:,8)=zero
          tmpKeq(:,9)=zero
          tmpKeq(:,10)=zero
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
                     tmpStats(isgbeg:isgend,:) = zero
                     if (iKeq.eq.1) then
                       tmpKeq(isgbeg:isgend,:) = zero
                     endif
                  enddo
               endif
               
               itkbeg = itkbeg + 4 + 2*numseg
               
            enddo !itask
            
      endif  ! numpe.gt.1


c     accumulate sum of sons to the fathers
      maxsz = iConsStressSz
      if (iKeq.eq.1) maxsz = max(10,iConsStressSz)
      do n=1,maxsz
         tmpStatsf = zero
         tmpKeqf = zero
         do i = 1,nshg
            ifathi=ifath(i)
            if (n.le.iConsStressSz) then
               tmpStatsf(ifathi) = tmpStatsf(ifathi) 
     &                                       + tmpStats(i,n)
            endif
            if (iKeq.eq.1) then
               tmpKeqf(ifathi) = tmpKeqf(ifathi) 
     &                             + tmpKeq(i,n)
            endif          
         enddo

c     Now  the true fathers and serrogates combine results and update
c     each other.
c     
         if(numpe .gt. 1) then
            call drvAllreduce(tmpStatsf, tmpStatsft,nfath)
            if (iKeq.eq.1) then
               call drvAllreduce(tmpKeqf, tmpKeqft,nfath)
            endif
         else
            tmpStatsft=tmpStatsf
            if (iKeq.eq.1) then
               tmpKeqft=tmpKeqf  
            endif
         endif
c     divide by # of sons to get average father for this step
c
         if (ndistsons.eq.nfath) then 
            tmpStatsft = tmpStatsft * rinvsons
         else if (ndistsons.eq.1) then
            invtmp = rinvsons(1)
            tmpStatsft = tmpStatsft * invtmp
         endif
        
         if (iKeq.eq.1) then
           if (ndistsons.eq.nfath) then
            tmpKeqft = tmpKeqft * rinvsons 
           else if (ndistsons.eq.1) then
            invtmp =  rinvsons(1)
            tmpKeqft = tmpKeqft * invtmp
           endif
         endif

         if (myrank.eq.master) then
            if (n.le.iConsStressSz) then
              stsBar(:,n)=tfact*tmpStatsft+(one-tfact)*stsBar(:,n)
            endif
            if (iKeq.eq.1) then
              stsBarKeq(:,n)=tfact*tmpKeqft+(one-tfact)*stsBarKeq(:,n)
            endif
         endif

       enddo


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
              
!      write (irstou,*) nfath, lstep
      if((itwmod.gt.0) .or. (irscale.ge.0)
     &   .or. (ispanAvg.eq.1)) then
           do i=1,nfath      
             if (iConsStress.eq.0) then 
               write (irstou,*) stsBar(i,1),stsBar(i,2),stsBar(i,3),
     &                          stsBar(i,4),stsBar(i,5),stsBar(i,6)
             else
               write (irstou,*) stsBar(i,1),stsBar(i,2),stsBar(i,3),
     &                          stsBar(i,4),stsBar(i,5),stsBar(i,6),
     &                          stsBar(i,7),stsBar(i,8),stsBar(i,9)

             endif
           enddo
      endif
      close (irstou)

      return

996   call error ('stsBar  ','opening ', irstou)

      end subroutine wstsBar

      subroutine wstsBarKeq(ifail)

      use spanStats
      include "common.h"
  
      character*25 fname1,fmt1
      character*5 cname
      integer itmp, irstou

      itmp = 1
      if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
      write (fmt1,"('(''stsbarKeq.'',i',i1,',1x)')") itmp
      write (fname1,fmt1) lstep
      fname1 = trim(fname1) // cname(myrank+1)
c     
      open (unit=irstou, file=fname1, status='unknown',
     &     err=996)
              
!      write (irstou,*) nfath, lstep
      if((ispanAvg.eq.1).and.(iKeq.eq.1)) then
        do i=1,nfath            
          write (irstou,*) stsBarKeq(i,1),stsBarKeq(i,2),stsBarKeq(i,3),
     &                     stsBarKeq(i,4),stsBarKeq(i,5),stsBarKeq(i,6),
     &                     stsBarKeq(i,7),stsBarKeq(i,8),stsBarKeq(i,9),
     &                     stsBarKeq(i,10)
        enddo
      endif
      close (irstou)

      return

996   call error ('stsBarKeq  ','opening ', irstou)

      end subroutine wstsBarKeq
