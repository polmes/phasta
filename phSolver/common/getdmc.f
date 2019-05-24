      subroutine getdmc (y,      shgl,      shp, 
     &                   iper,   ilwork,    
     &                   nsons,  ifath,     x )

      use pointer_data

      use quadfilt   ! This module gives us shglf(maxtp,nsd,maxsh,ngaussf),
c                    shpf(maxtp,maxsh,ngaussf), and Qwtf(maxtp,ngaussf). 
c                    Shpf and shglf are the shape funciotns and their 
c                    gradient evaluated using the quadrature rule desired 
c                    for computing the dmod. Qwt contains the weights of the 
c                    quad. points. 

      use eblock
      use lesArrs

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      type (LocalBlkData) blk

c
      dimension fres(nshg,24),         fwr(nshg),
     &          strnrm(nshg),
     &          xnum(nshg),           xden(nshg),
     &          xmij(nshg,6),         xlij(nshg,6),
     &          xnude(nfath,2),        xnuder(nfath,2),
     &          nsons(ndistsons),
     &          strl(numel,maxnint),           
     &          y(nshg,5), 
     &          ifath(nshg),          iper(nshg),
     &          ilwork(nlwork),!        xmudmi(numel,ngauss),
     &          x(numnp,3),
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), shp(MAXTOP,maxsh,MAXQPT)
      real*8    xi, surf, BLedge, MeshBLedge
c$$$     &          ,xnutf(nfath)  must be uncommmented for diags at bottom
c
c
c   setup the weights for time averaging of cdelsq
c
      if (irunTave.eq.0) then
        denom=max(1.0d0*(lstep),one)
        if(dtavei.lt.0) then
          wcur=one/denom
        else
          wcur=dtavei
        endif  
        whist=1.0-wcur
      else ! if irunTave=1 then doing running time average
        denom = max(1,lstep-irunTaveSt)
        wcur = one/denom
        whist = one-wcur  
      endif
c
c  hack in an interesting velocity field (uncomment to test dmod)
c
c      y(:,5) = 1.0  
c      y(:,2) = 2.0*x(:,1) - 3*x(:,2) 
c      y(:,3) = 3.0*x(:,1) + 4.0*x(:,2)
c      y(:,4) = 4.0*x(:,1) + x(:,2) + x(:,3) 
c      y(:,1) = Rgas * y(:,5) ! Necessary to make model suitable suitable
                             ! for the incompressible case.
c

     
      fres = zero
c
c       y(:,5) = 1.0  
c       y(:,1) = 2.0*x(:,1) - 3*x(:,2) 
c       y(:,2) = 3.0*x(:,1) + 4.0*x(:,2)
c       y(:,3) = 4.0*x(:,1) + x(:,2) + x(:,3) 
c       y(:,4) = 1.0    ! Necessary to make model suitable suitable
                             ! for the incompressible case.
c


      intrul=intg(1,itseq)
      intind=intpt(intrul)

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1
        blk%n   = lcblk(5,iblk) ! no. of vertices per element
        blk%s   = lcblk(10,iblk)
        blk%e   = lcblk(1,iblk+1) - iel 
        blk%l = lcblk(3,iblk)
        blk%g = nint(blk%l)
        blk%o = lcblk(4,iblk)
        blk%i = lcblk(1,iblk)

        ngauss = nint(lcsyst)
        ngaussf = nintf(lcsyst)
        
        call asithf (blk, y, x, strl(iel:inum,:), mien(iblk)%p, fres, 
     &               shglf(lcsyst,:,1:nshl,:),
     &               shpf(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf))

      enddo
c
 
      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1
        blk%n   = lcblk(5,iblk) ! no. of vertices per element
        blk%s   = lcblk(10,iblk)
        blk%e   = lcblk(1,iblk+1) - iel 
        blk%l = lcblk(3,iblk)
        blk%g = nint(blk%l)
        blk%o = lcblk(4,iblk)
        blk%i = lcblk(1,iblk)
        
        ngauss = nint(lcsyst)
        ngaussf = nintf(lcsyst)

        if (ngaussf .ne. ngauss) then
        call getstrl (blk, y, x,      mien(iblk)%p,  
     &               strl(iel:inum,:), shgl(lcsyst,:,1:nshl,:),
     &               shp(lcsyst,1:nshl,:))
        endif

      enddo
c
c
C must fix for abc and dynamic model
c      if(iabc==1)   !are there any axisym bc's
c     &      call rotabc(res, iBC,  'in ')
c
      if(numpe>1) then
         call commu (fres(:,1:4), ilwork, 4, 'in ')
         call commu (fres(:,5:8), ilwork, 4, 'in ')
         call commu (fres(:,9:12), ilwork, 4, 'in ')
         call commu (fres(:,13:16), ilwork, 4, 'in ')
         call commu (fres(:,17:20), ilwork, 4, 'in ')
         call commu (fres(:,21:24), ilwork, 4, 'in ')
      endif
!proc-masters = proc-masters + proc-slaves
c 
c account for periodicity in filtered variables
c
      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(i,:) = fres(i,:) + fres(j,:) ! masters = masters + slaves
        endif
      enddo
      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(j,:) = fres(i,:)  !slaves get copy of the complete master
        endif
      enddo

      if (nfath .eq. 0) then
         if(numpe>1)   call commu (fres, ilwork, 24, 'out')
c again we will need a  call
c      if(iabc==1)   !are there any axisym bc's
c     &   rotabc(y, iBC,  'out')
      endif

      fres(:,23) = one / fres(:,23)
      do j = 1,22
        fres(:,j) = fres(:,j) * fres(:,23)  !"solve"
      enddo
c     fres(:,24) = fres(:,24) * fres(:,23)
c
c.....at this point fres is really all of our filtered quantities
c     at the nodes
c

      strnrm = sqrt( 
     &  two * (fres(:,10)**2 + fres(:,11)**2 + fres(:,12)**2)
     &  + four * ( fres(:,13)**2 + fres(:,14)**2 + fres(:,15)**2 ) )

      fwr = fwr1 * fres(:,22) * strnrm



c... There is a difference in xmij(:,1:3) between the compressible 
c... and incompressible forms of the model.


      if(matflg(1,1).eq.0) then ! compressible

         xmij(:,1) = -fwr
     &             * pt33 * (two*fres(:,10) - fres(:,11) - fres(:,12))
     &             + pt33 * (two*fres(:,16) - fres(:,17) - fres(:,18))
         xmij(:,2) = -fwr
     &             * pt33 * (two*fres(:,11) - fres(:,10) - fres(:,12))
     &             + pt33 * (two*fres(:,17) - fres(:,16) - fres(:,18))
         xmij(:,3) = -fwr
     &             * pt33 * (two*fres(:,12) - fres(:,10) - fres(:,11))
     &             + pt33 * (two*fres(:,18) - fres(:,16) - fres(:,17))

      else

         xmij(:,1) = -fwr
     &             * fres(:,10) + fres(:,16)
         xmij(:,2) = -fwr
     &             * fres(:,11) + fres(:,17) 
         xmij(:,3) = -fwr
     &             * fres(:,12) + fres(:,18) 

      endif      

      xmij(:,4) = -fwr * fres(:,13) + fres(:,19)
      xmij(:,5) = -fwr * fres(:,14) + fres(:,20)
      xmij(:,6) = -fwr * fres(:,15) + fres(:,21)

      fres(:,22) = one / fres(:,22)

      xlij(:,1) = fres(:,4) - fres(:,1) * fres(:,1) * fres(:,22)
      xlij(:,2) = fres(:,5) - fres(:,2) * fres(:,2) * fres(:,22)
      xlij(:,3) = fres(:,6) - fres(:,3) * fres(:,3) * fres(:,22)
      xlij(:,4) = fres(:,7) - fres(:,1) * fres(:,2) * fres(:,22)
      xlij(:,5) = fres(:,8) - fres(:,1) * fres(:,3) * fres(:,22)
      xlij(:,6) = fres(:,9) - fres(:,2) * fres(:,3) * fres(:,22)

      xnum =        xlij(:,1) * xmij(:,1) + xlij(:,2) * xmij(:,2) 
     &                                    + xlij(:,3) * xmij(:,3)
     &     + two * (xlij(:,4) * xmij(:,4) + xlij(:,5) * xmij(:,5)
     &                                    + xlij(:,6) * xmij(:,6))
      xden =        xmij(:,1) * xmij(:,1) + xmij(:,2) * xmij(:,2) 
     &                                    + xmij(:,3) * xmij(:,3)
     &     + two * (xmij(:,4) * xmij(:,4) + xmij(:,5) * xmij(:,5)
     &                                    + xmij(:,6) * xmij(:,6))
      xden = two * xden

c  zero on processor periodic nodes so that they will not be added twice
      do j = 1,numnp
          i = iper(j)
          if (i .ne. j) then
            xnum(j) = zero
            xden(j) = zero
          endif
      enddo

      if (numpe.gt.1 .and. nsons(1).gt.1) then

         numtask = ilwork(1)
         itkbeg = 1
       
c zero the nodes that are "solved" on the other processors  
         do itask = 1, numtask

            iacc   = ilwork (itkbeg + 2)
            numseg = ilwork (itkbeg + 4)

            if (iacc .eq. 0) then
               do is = 1,numseg
                  isgbeg = ilwork (itkbeg + 3 + 2*is)
                  lenseg = ilwork (itkbeg + 4 + 2*is)
                  isgend = isgbeg + lenseg - 1
                  xnum(isgbeg:isgend) = zero
                  xden(isgbeg:isgend) = zero
               enddo
            endif
            
            itkbeg = itkbeg + 4 + 2*numseg
            
         enddo
         
      endif
c
c Description of arrays.   Each processor has an array of length equal
c to the total number of fathers times 2 xnude(nfathers,2). One to collect 
c the numerator and one to collect the denominator.  There is also an array
c of length nshg on each processor which tells the father number of each
c on processor node, ifath(nnshg).  Finally, there is an arry of length
c nfathers to tell the total (on all processors combined) number of sons
c for each father. 
c
c  Now loop over nodes and accumlate the numerator and the denominator
c  to the father nodes.  Only on processor addition at this point.
c  Note that serrogate fathers are collect some for the case where some
c  sons are on another processor
c
      ihomog=1
      if(maxval(nsons).eq.1) ihomog=0
      if(ihomog.eq.1) then
         xnude = zero
         do i = 1,nshg
            xnude(ifath(i),1) = xnude(ifath(i),1) + xnum(i)
            xnude(ifath(i),2) = xnude(ifath(i),2) + xden(i)
         enddo
      else                      ! no homogeneity assumed
         xnude(:,1)=xnum(:)
         xnude(:,2)=xden(:)
      endif

c
c Now  the true fathers and serrogates combine results and update
c each other.
c
c ONLY DO THIS IF  1) multiprocessors AND 2) homogenous directions exist 
  
      if(numpe .gt. 1 .and. ihomog.eq.1 )then
         call drvAllreduce(xnude, xnuder,2*nfath)
c
c  xnude is the sum of the sons for each father on this processor
c
c  xnuder is the sum of the sons for each father on all processor combined
c  (the same as if we had not partitioned the mesh for each processor)
c
c   For each father we have precomputed the number of sons (including
c   the sons off processor). 
c
c   Now divide by number of sons to get the average (not really necessary
c   for dynamic model since ratio will cancel nsons at each father)
c
c         xnuder(:,1) = xnuder(:,1) ! / nsons(:)
c         xnuder(:,2) = xnuder(:,2) ! / nsons(:)
c
c  the next line is c \Delta^2
c
            if (irunTave.eq.0) then 
              numNden(:,1) = whist*numNden(:,1)+wcur*xnuder(ifath(:),1)
              numNden(:,2) = whist*numNden(:,2)+wcur*xnuder(ifath(:),2)
              cdelsq(:,1) = numNden(:,1) / (numNden(:,2) + 1.d-09)
            else
              cdelsq(:,2) = whist*cdelsq(:,2)+wcur*xnuder(ifath(:),1)
              cdelsq(:,3) = whist*cdelsq(:,3)+wcur*xnuder(ifath(:),2)
              cdelsq(:,1) = cdelsq(:,2) / (cdelsq(:,3) + 1.d-09)
            endif
c           Perform clipping of cdelsq if desired
            if (iclipCdelsq.eq.1) then
             do ii=1,nshg
               xi = x(ii,1)
               surf = 0.0777240d0*exp(-(xi/0.178308)**2)
               BLedge = a0 + a1*cos(xi*w) + b1*sin(xi*w) + 
     &                  a2*cos(2.0*xi*w) + b2*sin(2.0*xi*w) + a3*cos(3.0*xi*w) + 
     &                  b3*sin(3.0*xi*w) + a4*cos(4.0*xi*w) + b4*sin(4.0*xi*w) + 
     &                  a5*cos(5.0*xi*w) + b5*sin(5.0*xi*w) + a6*cos(6.0*xi*w) +
     &                  b6*sin(6.0*xi*w) + a7*cos(7.0*xi*w) + b7*sin(7.0*xi*w) + 
     &                  a8*cos(8.0*xi*w) + b8*sin(8.0*xi*w)
               MeshBLedge = surf+1.40d0*(BLedge-surf)
               if (x(ii,2).gt.MeshBLedge) cdelsq(ii,1) = zero 
             enddo
            endif
c  note that we have whist and wcur in here to allow for both time
c  averaging to be used in conjunction with spatial homogenous averaging

      else
c     
c     the next line is c \Delta^2, not nu_T but we want to save the
c     memory
c     
c$$$            write(540+myrank,555) (xnude(j+500,2),j=1,5)
            if (irunTave.eq.0) then 
              numNden(:,1) = whist*numNden(:,1)+wcur*xnude(ifath(:),1)
              numNden(:,2) = whist*numNden(:,2)+wcur*xnude(ifath(:),2)
              cdelsq(:,1) = numNden(:,1) / (numNden(:,2) + 1.d-09)
            else
              cdelsq(:,2) = whist*cdelsq(:,2)+wcur*xnude(ifath(:),1)
              cdelsq(:,3) = whist*cdelsq(:,3)+wcur*xnude(ifath(:),2)
              cdelsq(:,1) = cdelsq(:,2) / (cdelsq(:,3) + 1.d-09)
            endif
c           Perform clipping of cdelsq if desired
            if (iclipCdelsq.eq.1) then
             do ii=1,nshg
               xi = x(ii,1)
               surf = 0.0777240d0*exp(-(xi/0.178308)**2)
               BLedge = a0 + a1*cos(xi*w) + b1*sin(xi*w) + 
     &                  a2*cos(2.0*xi*w) + b2*sin(2.0*xi*w) + a3*cos(3.0*xi*w) + 
     &                  b3*sin(3.0*xi*w) + a4*cos(4.0*xi*w) + b4*sin(4.0*xi*w) + 
     &                  a5*cos(5.0*xi*w) + b5*sin(5.0*xi*w) + a6*cos(6.0*xi*w) +
     &                  b6*sin(6.0*xi*w) + a7*cos(7.0*xi*w) + b7*sin(7.0*xi*w) + 
     &                  a8*cos(8.0*xi*w) + b8*sin(8.0*xi*w)
               MeshBLedge = surf+1.40d0*(BLedge-surf)
               if (x(ii,2).gt.MeshBLedge) cdelsq(ii,1) = zero 
             enddo
            endif
            
      endif
 555  format(5(2x,e14.7))

c     Write some stats
      tmp1 =  MINVAL(cdelsq(:,1))
      tmp2 =  MAXVAL(cdelsq(:,1))
      if(numpe>1) then
         call MPI_REDUCE (tmp1, tmp3, 1,MPI_DOUBLE_PRECISION,
     &        MPI_MIN, master, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE (tmp2, tmp4, 1, MPI_DOUBLE_PRECISION,
     &        MPI_MAX, master, MPI_COMM_WORLD, ierr)
         tmp1=tmp3
         tmp2=tmp4
      endif
      if (myrank .EQ. master) then !print CDelta^2 range
         write(34,*)lstep+1,tmp1,tmp2
         call flush(34)
         ! fort.34 contains the global min and max of cdelsq for each time step
      endif
c      if (myrank .eq. master) then
c         write(*,*)'xnut=',sum(cdelsq)/nshg
c         !write(*,*) 'cdelsq=', cdelsq(1),cdelsq(2)        
c      endif


c
c ... Another loop over element block to compute the eddy viscosity 
      do iblk = 1,nelblk
         lcsyst = lcblk(3,iblk)
         iel  = lcblk(1,iblk)
         npro = lcblk(1,iblk+1) - iel
         lelCat = lcblk(2,iblk)
         inum  = iel + npro - 1
         blk%n   = lcblk(5,iblk) ! no. of vertices per element
         blk%s   = lcblk(10,iblk)
         blk%e   = lcblk(1,iblk+1) - iel 
         blk%l = lcblk(3,iblk)
         blk%g = nint(blk%l)
         blk%o = lcblk(4,iblk)
         blk%i = lcblk(1,iblk)
         
         ngauss = nint(lcsyst)

         ! strl is actually the norm of the strain |S|
         ! xmudmi is the eddy viscosity nu_T=cdelsq*strl
         call scatnu (blk, mien(iblk)%p, strl(iel:inum,:), 
     &        mxmudmi(iblk)%p,shp(lcsyst,1:nshl,:))
      enddo

c     Write some eddy viscosity stats
!      tmp1 =  MINVAL(xmudmi)
!      tmp2 =  MAXVAL(xmudmi)
!      if(numpe>1) then
!        call MPI_REDUCE (tmp1, tmp3, 1, MPI_DOUBLE_PRECISION,
!     &                 MPI_MIN, master, MPI_COMM_WORLD, ierr)
!        call MPI_REDUCE (tmp2, tmp4, 1, MPI_DOUBLE_PRECISION,
!     &                 MPI_MAX, master, MPI_COMM_WORLD, ierr)
!        tmp1=tmp3
!        tmp2=tmp4
!      endif
!      if (myrank .EQ. master) then
!        write(35,*) lstep+1,tmp1,tmp2
!        call flush(35)
!        ! fort.35 contains global max and min of eddy viscosity for each time step
!      endif

c$$$c
c$$$c  if flag set, write a restart file with info (reuse xmij's memory)
c$$$c
c$$$      if(irs.eq.11) then
c$$$         lstep=999
c$$$         xmij(:,1)=xnum(:)
c$$$         xmij(:,2)=xden(:)
c$$$         xmij(:,3)=cdelsq(:)
c$$$         xmij(:,5)=xlij(:,4)    !leave M_{12} in 4 and put L_{12} here
c$$$         call restar('out ',xmij,xlij) !also dump all of L_{ij} in ac
c$$$         stop
c$$$      endif
c$$$      if(lhs.eq.1) then
c$$$         lstepsafe=lstep
c$$$         lstep=999
c$$$
c$$$         xmij(:,1)=xnum(:)
c$$$         xmij(:,2)=xden(:)
c$$$         xmij(:,3)=cdelsq(:)
c$$$         xmij(:,4)=xmij(:,6)    !put M_{23} in 4 
c$$$         xmij(:,5)=xlij(:,6)    !put L_{23} here
c$$$         call restar('out ',xmij,xlij) !also dump all of L_{ij} in ac
c$$$         lstep=lstepsafe
c$$$      endif
c$$$c
c$$$c  local clipping moved to scatnu with the creation of mxmudmi pointers
c$$$c
c$$$c$$$      rmu=datmat(1,2,1)
c$$$c$$$      xmudmi=min(xmudmi,1000.0*rmu) !don't let it get larger than 1000 mu
c$$$c$$$      xmudmi=max(xmudmi, -rmu) ! don't let (xmudmi + mu) < 0
c$$$c      stop !uncomment to test dmod
c$$$c
c$$$
c$$$
c$$$c  write out the nodal values of xnut (estimate since we don't calc strain
c$$$c  there and must use the filtered strain).
c$$$c
c$$$
c$$$      if ((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)) then
c$$$c
c$$$c  collect the average strain into xnude(2)
c$$$c
c$$$         xnude(:,2) = zero
c$$$         do i = 1,numnp
c$$$            xnude(ifath(i),2) = xnude(ifath(i),2) + strnrm(i)
c$$$         enddo
c$$$
c$$$         if(numpe .gt. 1 .and. ihomog.eq.1 )then
c$$$            call drvAllreduce(xnude(:,2), xnuder(:,2),nfath)
c$$$         else
c$$$            xnuder=xnude
c$$$         endif
c$$$c     
c$$$c          nut= cdelsq    * |S|
c$$$c 
c$$$         xnutf=xnuder(:,1)*xnuder(:,2)/nsons(:) ! off by one
c$$$c
c$$$c  collect the x and y coords into xnude
c$$$c
c$$$         xnude = zero
c$$$         do i = 1,numnp
c$$$            xnude(ifath(i),1) = xnude(ifath(i),1) + x(i,3)
c$$$            xnude(ifath(i),2) = xnude(ifath(i),2) + x(i,2)
c$$$         enddo
c$$$
c$$$         if(numpe .gt. 1 .and. nsons(1).ne.1 )then
c$$$            call drvAllreduce(xnude, xnuder,2*nfath)
c$$$            xnuder(:,1)=xnuder(:,1)/nsons(:)
c$$$            xnuder(:,2)=xnuder(:,2)/nsons(:)
c$$$         else
c$$$            xnuder=xnude
c$$$         endif
c$$$c
c$$$c  xnude is the sum of the sons for each father on this processor
c$$$c
c$$$         if((myrank.eq.master)) then
c$$$            do i=1,nfath      ! cdelsq   * |S|
c$$$               write(444,*) xnuder(i,1),xnuder(i,2),xnutf(i)
c$$$            enddo
c$$$            call flush(444)
c$$$         endif
c$$$      endif

      return
      end
