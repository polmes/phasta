c---------------------------------------------------------------------
c     
c     drvftools.f : Bundle of Fortran driver routines for ftools.f
c     
c     Each routine is to be called by les**.c
c     
c---------------------------------------------------------------------
c     
c----------------
c     drvLesPrepDiag
c----------------
c     
      subroutine drvlesPrepDiag ( flowDiag, ilwork,
     &                            iBC,      BC,      iper,
     &                            rowp,     colm,    
     &                            lhs16,    lhsP)   !keeping it in interface for a while
c     
      use pointer_data
      use pvsQbi
      use convolImpFlow !brings in the current part of convol coef for imp BC
      use convolRCRFlow !brings in the current part of convol coef for RCR BC

      include "common.h"
      include "mpif.h"
c     
      real*8 lhs16(16,nnz_tot)
      dimension flowDiag(nshg,4), ilwork(nlwork)
      dimension iBC(nshg), iper(nshg), BC(nshg,ndofBC)
      real*8 lhsP(4,nnz_tot)
      integer rowp(nnz_tot),  colm(nshg+1)
      integer	n,	k
c
      integer sparseloc
      rstart=TMRC()
c
c     
c.... Clear the flowdiag
c
      if((flmpl.eq.1).or.(ipord.gt.1)) then
         do n = 1, nshg
            k = sparseloc( rowp(colm(n)), colm(n+1)-colm(n), n )
     &       + colm(n)-1
c     
            flowdiag(n,1) = lhs16(1,k) !lhsK(1,k)
            flowdiag(n,2) = lhs16(6,k) !lhsK(5,k)
            flowdiag(n,3) = lhs16(11,k) !lhsK(9,k)
c     
            flowdiag(n,4) = lhsP(4,k)
         enddo
      else
        flowDiag = zero
        do n = 1, nshg  ! rowsum put on the diagonal instead of diag entry
           do k=colm(n),colm(n+1)-1

c
              flowdiag(n,1) = flowdiag(n,1) + abs(lhs16(1,k)) 
c     &                          + lhs16(2,k) + lhs16(3,k)
              flowdiag(n,2) = flowdiag(n,2) + abs(lhs16(6,k)) 
c     &                          + lhs16(5,k) + lhs16(7,k)
              flowdiag(n,3) = flowdiag(n,3) + abs(lhs16(11,k)) 
c     &                          + lhsK(9,k) + lhs16(10,k)
c
              flowdiag(n,4) = flowdiag(n,4) + abs(lhsP(4,k)) 
           enddo
           flowdiag(n,:)=flowdiag(n,:)*pt33
        enddo
      endif
      if(ipvsq.ge.3) then ! for first cut only do diagonal extraction
 ! this is not yet correct for multi procs I suspect if partition
 ! boundary cuts a p=QR face
         tfact=alfi * gami * Delt(1)
         do n=1,nshg
            if(numResistSrfs.gt.zero) then
               do k = 1,numResistSrfs
                  if (nsrflistResist(k).eq.ndsurf(n)) then
                     irankCoupled=k      
                     flowdiag(n,1:3) = flowdiag(n,1:3)
     &               + tfact*ValueListResist(irankCoupled)*
     &               NABI(n,:)*NABI(n,:)
                     exit
                  endif
               enddo
            elseif(numImpSrfs.gt.zero) then
               do k = 1,numImpSrfs
                  if (nsrflistImp(k).eq.ndsurf(n)) then
                     irankCoupled=k      
                     flowdiag(n,1:3) = flowdiag(n,1:3)
     &               + tfact*ImpConvCoef(ntimeptpT+2,irankCoupled)*
     &               NABI(n,:)*NABI(n,:)
                     exit
                  endif
               enddo
            elseif(numRCRSrfs.gt.zero) then
               do k = 1,numRCRSrfs
                  if (nsrflistRCR(k).eq.ndsurf(n)) then
                     irankCoupled=k      
                     flowdiag(n,1:3) = flowdiag(n,1:3)
     &               + tfact*RCRConvCoef(lstep+2,irankCoupled)* !check lstep+2 if restart from t.ne.0
     &               NABI(n,:)*NABI(n,:)
                     exit
                  endif
               enddo
            endif
         enddo
      endif
c     

c
      if(iabc==1)    !are there any axisym bc's
     &      call rotabc(flowdiag, iBC, 'in ')
c

c     
c.... communicate : add the slaves part to the master's part of flowDiag
c     
        if (numpe > 1) then 
           call commu (flowDiag, ilwork, nflow, 'in ') 
        endif
c
c.... satisfy the boundary conditions on the diagonal
c
        call bc3diag(iBC, BC,  flowDiag)
c
c     
c.... on processor periodicity was not taken care of in the setting of the 
c     boundary conditions on the matrix.  Take care of it now.
c
        call bc3per(iBC,  flowDiag, iper, ilwork, 4)
c
c... slaves and masters have the correct values
c
c     
c.... Calculate square root
c     
        do i = 1, nshg
           do j = 1, nflow
              if (flowDiag(i,j).ne.0) 
     &             flowDiag(i,j) = 1. / sqrt(abs(flowDiag(i,j)))
           enddo
        enddo
        rspmvD=rspmvD+TMRC()-rstart
        ispmvD=ispmvD+1
c     
        return
        end

c     
c-------------
c     drvsclrDiag
c-------------
c     
      subroutine drvsclrDiag ( sclrDiag, ilwork, iBC, BC, iper, 
     &                         rowp,     colm,   lhsS )
c     
      use pointer_data
      include "common.h"
      include "mpif.h"
c     
      integer  ilwork(nlwork),    iBC(nshg),     iper(nshg),
     &         rowp(nnz_tot),    colm(nshg+1)

      real*8   sclrDiag(nshg),    lhsS(nnz_tot), BC(nshg,ndofBC)
      integer sparseloc

      sclrDiag = zero
      do n = 1, nshg
         k = sparseloc( rowp(colm(n)), colm(n+1)-colm(n), n ) 
     &                               + colm(n)-1
c
         sclrDiag(n) = lhsS(k)
      enddo
c     
c.... communicate : add the slaves part to the master's part of sclrDiag
c     
        if (numpe > 1) then 
           call commu (sclrDiag, ilwork, 1, 'in ') 
        endif
c
c.... satisfy the boundary conditions on the diagonal
c
        call bc3SclrDiag(iBC,  sclrDiag)
c
c     
c.... on processor periodicity was not taken care of in the setting of the 
c     boundary conditions on the matrix.  Take care of it now.
c
        call bc3per(iBC,  sclrDiag, iper, ilwork, 1)
c
c... slaves and masters have the correct values
c
c     
c.... Calculate square root
c     
        do i = 1, nshg
           if (sclrDiag(i).ne.0) then
              sclrDiag(i) = 1. / sqrt(abs(sclrDiag(i)))
           endif
        enddo
c     
      return
      end

C============================================================================
C
C "fLesSparseApG":
C
C============================================================================
        subroutine fLesSparseApG(	col,	row,	pLhs,	
     &					p,	q,	nNodes,
     &                                  nnz_tot_hide )
c
c.... Data declaration
c
!        implicit none
      include "common.h"
        integer	nNodes, nnz_tot
        integer	col(nNodes+1),	row(nnz_tot)
        real*8	pLhs(4,nnz_tot),	p(nNodes),	q(nNodes,3)
c
        real*8	pisave
        integer	i,	j,	k
        rstart=TMRC()
c
c.... clear the vector
c
        do i = 1, nNodes
            q(i,1) = 0
            q(i,2) = 0
            q(i,3) = 0
        enddo
c
c.... Do an AP product
c
        do i = 1, nNodes
c
            pisave = p(i)
cdir$ ivdep
            do k = col(i), col(i+1)-1
              j = row(k) 
c
              q(j,1) = q(j,1) - pLhs(1,k) * pisave
              q(j,2) = q(j,2) - pLhs(2,k) * pisave
              q(j,3) = q(j,3) - pLhs(3,k) * pisave
            enddo
        enddo
c
c.... end
c
      rspmvG=rspmvG+TMRC()-rstart
      ispmvG=ispmvG+1
        return
        end

C============================================================================
C
C "fLesSparseApKG":   OPTION to TRANSPOSED P vector prior to call
C
C============================================================================

      subroutine fLesSparseApKG(col,	row,	lhs16,	pLhs,
     1                           p,	q,	nNodes,
     2                                  nnz_tot_hide ) 
c
c.... Data declaration
c
c      implicit none
      use pvsQbi
      include "common.h"
      integer	nNodes
      integer	col(nNodes+1),	row(nnz_tot)
      real*8 	p(nNodes,4),  q(nNodes,3)
      real*8    q3(3,nNodes), p3(3,nNodes)
      real*8    q4(4,nNodes), p4(4,nNodes)
      real*8    pLhs(4,nnz_tot)  ! old way passed in matrices
      real*8    lhs9(9,nnz_tot), lhs16(16,nnz_tot) ! needed for 3x3 options
c
      real*8	tmp1,	tmp2,	tmp3,	pisave
      integer	i,	j,	k
!
! This routine performs K.p_m + G.p_c
! The K.p_m is a 3x3 . 3x1 (standard). The second part is G.p_c is a 3x1 . 1x1
! Original code did not store G so  it gets pulled from PLHS but note
! PLHS in NOT  GoverC (4x1) rather it is {-G^T C}  (1x4).  What we
! see below is an on the fly negation and transpose (note j inplace summ) to
! accomplish + G p_c.  These options allow testing if this is more or less efficient 
! than directly computing and using the full matrix.
!
      rstart=TMRC()
      iwork=ieqswork/10
#ifndef USE_MKL
        iwork=5  ! safety in input selected MKL but code was not compiled for it 
#endif
! chosen: 0 as above,  1 as above^T, 2 use MKL for the K.p_m then OS G.p_c
! 3 MKL on 4x4, 4 use 4x4 ^T, 5 use 4x4  no ^T, 6, 0 using 4x4
      if(iwork.eq.0) then  ! {old way
        lhs9(1:3,:)=lhs16(1:3,:) ! left out of timer as we could store
        lhs9(4:6,:)=lhs16(5:7,:)
        lhs9(7:9,:)=lhs16(9:11,:)
        rdelta=TMRC()
        do i = 1, nNodes
          q(i,1) = 0
          q(i,2) = 0
          q(i,3) = 0
        enddo
c
c.... Do an AP product
c
        do i = 1, nNodes
c
            tmp1 = 0
            tmp2 = 0
            tmp3 = 0
            pisave   = p(i,4)
cdir$ ivdep
            do k = col(i), col(i+1)-1
                j = row(k) 
                tmp1 = tmp1
     1               + lhs9(1,k) * p(j,1)
     2               + lhs9(4,k) * p(j,2)
     3               + lhs9(7,k) * p(j,3)
                tmp2 = tmp2
     1               + lhs9(2,k) * p(j,1)
     2               + lhs9(5,k) * p(j,2)
     3               + lhs9(8,k) * p(j,3)
                tmp3 = tmp3
     1               + lhs9(3,k) * p(j,1)
     2               + lhs9(6,k) * p(j,2)
     3               + lhs9(9,k) * p(j,3)
                q(j,1) = q(j,1) - pLhs(1,k) * pisave
                q(j,2) = q(j,2) - pLhs(2,k) * pisave
                q(j,3) = q(j,3) - pLhs(3,k) * pisave
            enddo
            q(i,1) = q(i,1) + tmp1
            q(i,2) = q(i,2) + tmp2
            q(i,3) = q(i,3) + tmp3
        enddo
        rspmvphasta=rspmvphasta + TMRC()-rdelta
        ispmvphasta=ispmvphasta + 1
      endif !} original code in fast index on dof-HOLDER

      if(iwork.eq.1) then ! original transposed {
        do i = 1, nNodes
          p3(1,i)=p(i,1)
          p3(2,i)=p(i,2)
          p3(3,i)=p(i,3)
        enddo
        lhs9(1:3,:)=lhs16(1:3,:) ! left out of timer as we could store
        lhs9(4:6,:)=lhs16(5:7,:)
        lhs9(7:9,:)=lhs16(9:11,:)
c
c.... Do an AP product
c
        rdelta= TMRC()
        do i = 1, nNodes
          q3(1,i) = 0
          q3(2,i) = 0
          q3(3,i) = 0
        enddo
        do i = 1, nNodes
c
          tmp1 = 0
          tmp2 = 0
          tmp3 = 0
          pisave   = p(i,4)
cdir$ ivdep
          do k = col(i), col(i+1)-1
            j = row(k) 
            tmp1 = tmp1
     1           + lhs9(1,k) * p3(1,j)
     2           + lhs9(4,k) * p3(2,j)
     3           + lhs9(7,k) * p3(3,j)
            tmp2 = tmp2
     1           + lhs9(2,k) * p3(1,j)
     2           + lhs9(5,k) * p3(2,j)
     3           + lhs9(8,k) * p3(3,j)
            tmp3 = tmp3
     1           + lhs9(3,k) * p3(1,j)
     2           + lhs9(6,k) * p3(2,j)
     3           + lhs9(9,k) * p3(3,j)
            q3(1,j) = q3(1,j) - pLhs(1,k) * pisave
            q3(2,j) = q3(2,j) - pLhs(2,k) * pisave
            q3(3,j) = q3(3,j) - pLhs(3,k) * pisave
          enddo
! cannot accumulate into transpose due to sparse transposed G work
          q3(1,i) = q3(1,i) + tmp1
          q3(2,i) = q3(2,i) + tmp2
          q3(3,i) = q3(3,i) + tmp3
        enddo
        rspmvphasta=rspmvphasta + TMRC()-rdelta
        ispmvphasta=ispmvphasta + 1
        do i =1, nNodes ! transpose back from  dof_var first
          q(i,1)=q3(1,i)
          q(i,2)=q3(2,i)
          q(i,3)=q3(3,i)
        enddo
      endif !iwork 1 }
      if(iwork.eq.2) then ! { mkl sparse  3x3
        do i = 1, nNodes
          p3(1,i)=p(i,1)
          p3(2,i)=p(i,2)
          p3(3,i)=p(i,3)
        enddo
        lhs9(1:3,:)=lhs16(1:3,:)
        lhs9(4:6,:)=lhs16(5:7,:)
        lhs9(7:9,:)=lhs16(9:11,:)
        rdelta=TMRC()
#ifdef HAVE_MKL
        call mkl_dbsrgemv('N', nNodes, 3, lhs9, col, row, p3, q3)
#endif
        rspmvmkl=rspmvmkl + TMRC()-rdelta
        ispmvmkl=ispmvmkl + 1
        rdelta=TMRC()
        do i = 1, nNodes
c
          tmp1 = 0
          tmp2 = 0
          tmp3 = 0
          pisave   = p(i,4)
cdir$ ivdep
          do k = col(i), col(i+1)-1
            j = row(k) 
            q3(1,j) = q3(1,j) - pLhs(1,k) * pisave
            q3(2,j) = q3(2,j) - pLhs(2,k) * pisave
            q3(3,j) = q3(3,j) - pLhs(3,k) * pisave
          enddo
! cannot accumulate into transpose due to sparse transposed G work
          q3(1,i) = q3(1,i) + tmp1
          q3(2,i) = q3(2,i) + tmp2
          q3(3,i) = q3(3,i) + tmp3
        enddo
        rspmvphasta=rspmvphasta + TMRC()-rdelta
        ispmvphasta=ispmvphasta + 1
        do i =1, nNodes ! transpose back from  dof_var first
          q(i,1)=q3(1,i)
          q(i,2)=q3(2,i)
          q(i,3)=q3(3,i)
        enddo
!SLOWER ALT c
!SLOWER ALT c.... Do an AP product
!SLOWER ALT c
!SLOWER ALT         rdelta=TMRC()
!SLOWER ALT         do i = 1, nNodes
!SLOWER ALT           tmp1=0
!SLOWER ALT           tmp2=0
!SLOWER ALT           tmp3=0
!SLOWER ALT c
!SLOWER ALT cdir$ ivdep
!SLOWER ALT           do k = col(i), col(i+1)-1
!SLOWER ALT             j = row(k) 
!SLOWER ALT             tmp1 = tmp1 + lhs16(13,k) * p(j,4)
!SLOWER ALT             tmp2 = tmp2 + lhs16(14,k) * p(j,4)
!SLOWER ALT             tmp3 = tmp3 + lhs16(15,k) * p(j,4)
!SLOWER ALT           enddo
!SLOWER ALT ! accummulate into transpose to avoid transpose later
!SLOWER ALT           q(i,1) = q3(1,i) + tmp1
!SLOWER ALT           q(i,2) = q3(2,i) + tmp2
!SLOWER ALT           q(i,3) = q3(3,i) + tmp3
!SLOWER ALT         enddo
!SLOWER ALT         rspmvphasta=rspmvphasta + TMRC()-rdelta
      endif !} mkl sparse 3x3
      if(iwork.eq.3) then ! { mkl  with 4x4 matrix
        do i =1, nNodes  !mkl requires dof_var first
          p4(1,i)=p(i,1)
          p4(2,i)=p(i,2)
          p4(3,i)=p(i,3)
          p4(4,i)=p(i,4)
        enddo
        rdelta=TMRC()
#ifdef HAVE_MKL
        call mkl_dbsrgemv('N', nNodes, 4, lhs16, col, row, p4, q4)
#endif
        rspmvmkl=rspmvmkl + TMRC()-rdelta
        ispmvmkl=ispmvmkl + 1
        do i =1, nNodes ! transpose back mkl's dof_var first
          q(i,1)=q4(1,i)
          q(i,2)=q4(2,i)
          q(i,3)=q4(3,i)
        enddo
      endif !} end mkl 4x4
      if(iwork.eq.4) then !{ 4x4 transposed p
        do i =1, nNodes
          p4(1,i)=p(i,1)
          p4(2,i)=p(i,2)
          p4(3,i)=p(i,3)
          p4(4,i)=p(i,4)
        enddo
        rdelta= TMRC()
c
c.... Do an AP product
c
        do i = 1, nNodes
c
          tmp1 = 0
          tmp2 = 0
          tmp3 = 0
cdir$ ivdep
          do k = col(i), col(i+1)-1
            j = row(k) 
            tmp1 = tmp1
     1           + lhs16( 1,k) * p4(1,j)
     2           + lhs16( 5,k) * p4(2,j)
     3           + lhs16( 9,k) * p4(3,j)
     4           + lhs16(13,k) * p4(4,j)
            tmp2 = tmp2
     1           + lhs16( 2,k) * p4(1,j)
     2           + lhs16( 6,k) * p4(2,j)
     3           + lhs16(10,k) * p4(3,j)
     3           + lhs16(14,k) * p4(4,j)
            tmp3 = tmp3
     1           + lhs16( 3,k) * p4(1,j)
     2           + lhs16( 7,k) * p4(2,j)
     3           + lhs16(11,k) * p4(3,j)
     3           + lhs16(15,k) * p4(4,j)
          enddo
          q(i,1) = tmp1
          q(i,2) = tmp2
          q(i,3) = tmp3
        enddo
        rspmvphasta=rspmvphasta + TMRC()-rdelta
        ispmvphasta=ispmphasta + 1
      endif ! } 4x4 transposed 
      if(iwork.eq.5) then  ! { 4x4 used to do 3x4 with original p
        rdelta= TMRC()
c
c.... Do an AP product
c
        do i = 1, nNodes
c
          tmp1 = 0
          tmp2 = 0
          tmp3 = 0
cdir$ ivdep
          do k = col(i), col(i+1)-1
            j = row(k) 
            tmp1 = tmp1
     1           + lhs16( 1,k) * p(j,1)
     2           + lhs16( 5,k) * p(j,2)
     3           + lhs16( 9,k) * p(j,3)
     4           + lhs16(13,k) * p(j,4)
            tmp2 = tmp2
     1           + lhs16( 2,k) * p(j,1)
     2           + lhs16( 6,k) * p(j,2)
     3           + lhs16(10,k) * p(j,3)
     3           + lhs16(14,k) * p(j,4)
            tmp3 = tmp3
     1           + lhs16( 3,k) * p(j,1)
     2           + lhs16( 7,k) * p(j,2)
     3           + lhs16(11,k) * p(j,3)
     3           + lhs16(15,k) * p(j,4)
          enddo
          q(i,1) = tmp1
          q(i,2) = tmp2
          q(i,3) = tmp3
        enddo
        rspmvphasta=rspmvphasta + TMRC()-rdelta
        ispmvphasta=ispmphasta + 1
      endif ! } original p  with 4x4
      if(iwork.eq.6) then  ! {old way (using 4x3 matrix from a 4x4 data structure
        rdelta=TMRC()
        do i = 1, nNodes
          q(i,1) = 0
          q(i,2) = 0
          q(i,3) = 0
        enddo
c
c.... Do an AP product
c
        do i = 1, nNodes
c
            tmp1 = 0
            tmp2 = 0
            tmp3 = 0
            pisave   = p(i,4)
cdir$ ivdep
            do k = col(i), col(i+1)-1
                j = row(k) 
                tmp1 = tmp1
     1               + lhs16( 1,k) * p(j,1)
     2               + lhs16( 5,k) * p(j,2)
     3               + lhs16( 9,k) * p(j,3)
                tmp2 = tmp2
     1               + lhs16( 2,k) * p(j,1)
     2               + lhs16( 6,k) * p(j,2)
     3               + lhs16(10,k) * p(j,3)
                tmp3 = tmp3
     1               + lhs16( 3,k) * p(j,1)
     2               + lhs16( 7,k) * p(j,2)
     3               + lhs16(11,k) * p(j,3)
                q(j,1) = q(j,1) - lhs16( 4,k) * pisave
                q(j,2) = q(j,2) - lhs16( 8,k) * pisave
                q(j,3) = q(j,3) - lhs16(12,k) * pisave
            enddo
            q(i,1) = q(i,1) + tmp1
            q(i,2) = q(i,2) + tmp2
            q(i,3) = q(i,3) + tmp3
        enddo
        rspmvphasta=rspmvphasta + TMRC()-rdelta
        ispmvphasta=ispmphasta + 1
      endif !} original code, 4x4,  in fast index on dof-HOLDER

      if(ipvsq.ge.2) then
        tfact=alfi * gami * Delt(1)
        call ElmpvsQ(q,p,tfact)
      endif
c
c.... end
c
      rspmvKG=rspmvKG+TMRC()-rstart
      ispmvKG=ispmvKG+1
      return
      end


C============================================================================
C
C "fLesSparseApNGt":
C
C============================================================================

        subroutine fLesSparseApNGt(	col,	row,	pLhs,	
     1					p,	q,	nNodes,
     2                                  nnz_tot_hide   )
c
c.... Data declaration
c
        include "common.h"

        integer	col(nNodes+1),	row(nnz_tot)
        real*8	pLhs(4,nnz_tot),	p(nNodes,3),	q(nNodes)
c
        real*8	tmp
        integer	i,	j,	k
        rstart=TMRC()
c
c.... Do an AP product
c
        do i = nNodes, 1, -1
c
            tmp = 0
cdir$ ivdep
            do k = col(i), col(i+1)-1
        	j = row(k)
c
        	tmp = tmp
     1		    + pLhs(1,k) * p(j,1)
     2		    + pLhs(2,k) * p(j,2)
     3		    + pLhs(3,k) * p(j,3)
            enddo
            q(i) = tmp
        enddo
c
c.... end
c
        rspmvNGt=rspmvNGt+TMRC()-rstart
        ispmvNGt=ispmvNGt+1
        return
        end

C============================================================================
C
C "fLesSparseApNGtC":
C
C============================================================================

        subroutine fLesSparseApNGtC(	col,	row,	pLhs,	
     1					p,	q,	nNodes,
     2                                  nnz_tot_hide )
c
c.... Data declaration
c
!        implicit none
      include "common.h"
        integer	nNodes, nnz_tot
        integer	col(nNodes+1),	row(nnz_tot)
        real*8	pLhs(4,nnz_tot),	p(nNodes,4),	q(nNodes)
c
        real*8	tmp
        integer	i,	j,	k
        rstart=TMRC()
c
c.... Do an AP product
c
        do i = nNodes, 1, -1
c
            tmp = 0
cdir$ ivdep
            do k = col(i), col(i+1)-1
        	j = row(k)
c
        	tmp = tmp
     1		    + pLhs(1,k) * p(j,1)
     2		    + pLhs(2,k) * p(j,2)
     3		    + pLhs(3,k) * p(j,3)
     4		    + pLhs(4,k) * p(j,4)
            enddo
            q(i) = tmp
        enddo
c
c.... end
c
        rspmvNGtC=rspmvNGtC+TMRC()-rstart
        ispmvNGtC=ispmvNGtC+1
        return
        end

C============================================================================
C
C "fLesSparseApFull":
C
C============================================================================

      subroutine fLesSparseApFull(	col,	row,	lhs16,	pLhs,
     1					p,	q,	nNodes,
     2                                  nnz_tot_hide )
c
c.... Data declaration
c
c	implicit none
      use pvsQbi
      include "common.h"

      integer	nNodes, nnz_tot
      integer	col(nNodes+1),	row(nnz_tot)
      real*8	lhs9(9,nnz_tot),	pLhs(4,nnz_tot)
      real*8	lhs16(16,nnz_tot)
      real*8  p(nNodes,4),	q(nNodes,4)
      real*8  p4(4,nNodes),	q4(4,nNodes)
c
      real*8	tmp1,	tmp2,	tmp3,	tmp4
      integer	i,	j,	k
       
      iwork=ieqswork/10
! chosen: 0 original with matrices contracted back
! 3 MKL on 4x4, 4 use 4x4 ^T, 5 use 4x4  no ^T
!
! handle options available to KG but not for Full
!
      if((iwork.eq.1).or.(iwork.eq.2).or.(iwork.eq.6)) iwork=5 
#ifndef USE_MKL
        iwork=5  ! safety in input selected MKL but code was not compiled for it 
#endif
      if(iwork.eq.0) then !{ original alg with 3x3
        lhs9(1:3,:)=lhs16(1:3,:) ! left out of timer as we could store
        lhs9(4:6,:)=lhs16(5:7,:)
        lhs9(7:9,:)=lhs16(9:11,:)
        rstart=TMRC()
        do i = 1, nNodes
            q(i,1) = 0
            q(i,2) = 0
            q(i,3) = 0
        enddo
c
c.... Do an AP product
c
        do i = 1, nNodes
c
            tmp1 = 0
            tmp2 = 0
            tmp3 = 0
            tmp4 = 0
            pisave   = p(i,4)
cdir$ ivdep
            do k = col(i), col(i+1)-1
                j = row(k)
c
                tmp1 = tmp1
     1               + lhs9(1,k) * p(j,1)
     2               + lhs9(4,k) * p(j,2)
     3               + lhs9(7,k) * p(j,3)
                tmp2 = tmp2
     1               + lhs9(2,k) * p(j,1)
     2               + lhs9(5,k) * p(j,2)
     3               + lhs9(8,k) * p(j,3)
                tmp3 = tmp3
     1               + lhs9(3,k) * p(j,1)
     2               + lhs9(6,k) * p(j,2)
     3               + lhs9(9,k) * p(j,3)
c
                tmp4 = tmp4
     1               + pLhs(1,k) * p(j,1)
     2               + pLhs(2,k) * p(j,2)
     3               + pLhs(3,k) * p(j,3)
     4               + pLhs(4,k) * p(j,4)
c
                q(j,1) = q(j,1) - pLhs(1,k) * pisave
                q(j,2) = q(j,2) - pLhs(2,k) * pisave
                q(j,3) = q(j,3) - pLhs(3,k) * pisave
            enddo
            q(i,1) = q(i,1) + tmp1
            q(i,2) = q(i,2) + tmp2
            q(i,3) = q(i,3) + tmp3
            q(i,4) = tmp4
        enddo
        rspmvFull=rspmvFull+TMRC()-rstart
        ispmvFull=ispmvFull+1
      endif
      if(iwork.eq.3) then ! { mkl  with 4x4 matrix
        do i =1, nNodes  !mkl requires dof_var first
          p4(1,i)=p(i,1)
          p4(2,i)=p(i,2)
          p4(3,i)=p(i,3)
          p4(4,i)=p(i,4)
        enddo
        rstart=TMRC()
#ifdef HAVE_MKL
        call mkl_dbsrgemv('N', nNodes, 4, lhs16, col, row, p4, q4)
#endif
        rspmvFull=rspmvFull+TMRC()-rstart
        ispmvFull=ispmvFull+1
        do i =1, nNodes ! transpose back mkl's dof_var first
          q(i,1)=q4(1,i)
          q(i,2)=q4(2,i)
          q(i,3)=q4(3,i)
          q(i,4)=q4(4,i)
        enddo
      endif !} end mkl 4x4
      if(iwork.eq.4) then !{mkl 4x4
        do i =1, nNodes  !mkl requires dof_var first
          p4(1,i)=p(i,1)
          p4(2,i)=p(i,2)
          p4(3,i)=p(i,3)
          p4(4,i)=p(i,4)
        enddo
        rstart=TMRC()
c
c.... Do an AP product
c
        do i = 1, nNodes
c
            tmp1 = 0
            tmp2 = 0
            tmp3 = 0
            tmp4 = 0

cdir$ ivdep
            do k = col(i), col(i+1)-1
                j = row(k)
c
                tmp1 = tmp1
     1               + lhs16( 1,k) * p4(1,j)
     2               + lhs16( 5,k) * p4(2,j)
     3               + lhs16( 9,k) * p4(3,j)
     3               + lhs16(13,k) * p4(4,j)
                tmp2 = tmp2
     1               + lhs16( 2,k) * p4(1,j)
     2               + lhs16( 6,k) * p4(2,j)
     3               + lhs16(10,k) * p4(3,j)
     3               + lhs16(14,k) * p4(4,j)
                tmp3 = tmp3
     1               + lhs16( 3,k) * p4(1,j)
     2               + lhs16( 7,k) * p4(2,j)
     3               + lhs16(11,k) * p4(3,j)
     3               + lhs16(15,k) * p4(4,j)
c
                tmp4 = tmp4
     1               + lhs16( 4,k) * p4(1,j)
     2               + lhs16( 8,k) * p4(2,j)
     3               + lhs16(12,k) * p4(3,j)
     4               + lhs16(16,k) * p4(4,j)
c
            enddo
            q(i,1) = tmp1 !done inline to avoid transpose
            q(i,2) = tmp2
            q(i,3) = tmp3
            q(i,4) = tmp4
        enddo
        rspmvFull=rspmvFull+TMRC()-rstart
        ispmvFull=ispmvFull+1
!done inline        do i =1, nNodes ! transpose back mkl's dof_var first
!done inline          q(i,1)=q4(1,i)
!done inline          q(i,2)=q4(2,i)
!done inline          q(i,3)=q4(3,i)
!done inline          q(i,4)=q4(4,i)
!done inline        enddo
      endif
      if(iwork.eq.5) then !{
        rstart=TMRC()
c
c.... Do an AP product
c
        do i = 1, nNodes
c
            tmp1 = 0
            tmp2 = 0
            tmp3 = 0
            tmp4 = 0

cdir$ ivdep
            do k = col(i), col(i+1)-1
                j = row(k)
c
                tmp1 = tmp1
     1               + lhs16( 1,k) * p(j,1)
     2               + lhs16( 5,k) * p(j,2)
     3               + lhs16( 9,k) * p(j,3)
     3               + lhs16(13,k) * p(j,4)
                tmp2 = tmp2
     1               + lhs16( 2,k) * p(j,1)
     2               + lhs16( 6,k) * p(j,2)
     3               + lhs16(10,k) * p(j,3)
     3               + lhs16(14,k) * p(j,4)
                tmp3 = tmp3
     1               + lhs16( 3,k) * p(j,1)
     2               + lhs16( 7,k) * p(j,2)
     3               + lhs16(11,k) * p(j,3)
     3               + lhs16(15,k) * p(j,4)
c
                tmp4 = tmp4
     1               + lhs16( 4,k) * p(j,1)
     2               + lhs16( 8,k) * p(j,2)
     3               + lhs16(12,k) * p(j,3)
     4               + lhs16(16,k) * p(j,4)
c
            enddo
            q(i,1) = tmp1
            q(i,2) = tmp2
            q(i,3) = tmp3
            q(i,4) = tmp4
        enddo
        rspmvFull=rspmvFull+TMRC()-rstart
        ispmvFull=ispmvFull+1
      endif

      if(ipvsq.ge.2) then
        tfact=alfi * gami * Delt(1)
        call ElmpvsQ(q,p,tfact)
      endif
c
c.... end
c
        return
        end

C============================================================================
C
C "fLesSparseApSclr":
C
C============================================================================

        subroutine fLesSparseApSclr(	col,	row,	lhs,	
     1					p,	q,	nNodes,
     &                                  nnz_tot)
c
c.... Data declaration
c
        implicit none
        integer	nNodes, nnz_tot
        integer	col(nNodes+1),	row(nnz_tot)
        real*8	lhs(nnz_tot),	p(nNodes),	q(nNodes)
c
        real*8	tmp
        integer	i,	j,	k
c
c.... Do an AP product
c
        do i = nNodes, 1, -1
c
            tmp = 0
cdir$ ivdep
            do k = col(i), col(i+1)-1
        	tmp = tmp + lhs(k) * p(row(k))
            enddo
            q(i) = tmp
        enddo
c
c.... end
c
        return
        end
C============================================================================
        subroutine commOut(  global,  ilwork,  n, 
     &                       iper,    iBC, BC  )
        
        include "common.h"
        
        real*8  global(nshg,n), BC(nshg,ndofBC)
        integer ilwork(nlwork), iper(nshg), iBC(nshg)
c
        if ( numpe .gt. 1) then 
           call commu ( global, ilwork, n, 'out')
        endif
c
c     before doing AP product P must be made periodic
c     on processor slaves did not get updated with the 
c     commu (out) so do it here
c
        do i=1,n
           global(:,i) = global(iper(:),i)  ! iper(i)=i if non-slave so no danger
        enddo
c
c       slave has masters value, for abc we need to rotate it
c        (if this is a vector only no SCALARS)
        if((iabc==1) .and. (n.gt.1)) !are there any axisym bc's
     &     call rotabc(global, iBC,  'out')


c$$$        do j = 1,nshg
c$$$           if (btest(iBC(j),10)) then
c$$$              i = iper(j)
c$$$              res(j,:) = res(i,:) 
c$$$           endif
c$$$        enddo
        
        return 
        end

C============================================================================
        subroutine commIn(  global,  ilwork,  n, 
     &                      iper,    iBC, BC )
        
        include "common.h"
        
        real*8  global(nshg,n), BC(nshg,ndofBC)
        integer ilwork(nlwork), iper(nshg), iBC(nshg)
c
        if((iabc==1) .and. (n.gt.1)) !are there any axisym bc's
     &       call rotabc(global, iBC, 'in ')
c

        if ( numpe .gt. 1 ) then
           call commu ( global, ilwork, n, 'in ')
        endif
        	
        call bc3per ( iBC, global, iper, ilwork, n)
        
        return 
        end

