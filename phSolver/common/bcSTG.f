!-----------------------------------------------------------------------
!
!  This module handles synthetic turbulence generation at inflows and
!   later at LES switchover planes for ADES
!
!-----------------------------------------------------------------------
      module STG_BC

       real*8, allocatable :: dVect(:,:), sigVect(:,:)
       real*8, allocatable :: STGrnd(:,:), BCrestart(:,:) 
       real*8, allocatable :: phiVect(:),hx1(:),hx2(:),hy1(:),hy2(:)
       real*8, allocatable :: hz1(:),hz2(:),hMax(:)
       real*8, allocatable :: uBar(:,:), dUBardY(:), duBardYADJ(:)
       real*8, allocatable :: STGInflow(:,:)
       integer,allocatable :: n_hxl(:),n_hxh(:),n_hyl(:),n_hyh(:)
       integer,allocatable :: n_hzl(:),n_hzh(:)
       integer,allocatable :: stgSurf(:),mapBack(:)!,d2Wall(:) 
      
       integer STGmeth, nKWave,nNSurf
       real*8 alphaWave,alphaGR 
       real*8, allocatable :: le(:),lt(:), kCut(:)
       real*8 kMax,kMax_all,lCut,leMax,leMax_all,kMin
       integer iSTGInread, nPoints

!Flow Parameters
        real*8 reTau,uTau,deltaBLSTG,uNaught,kappa,b,nu
!$$$      integer itvn


        end

!----------------------------------------------------------------------
!
!    Create Initial STG random vectors
!   Commented print *, statements are used for Random vector debug/varification
!-----------------------------------------------------------------------
        subroutine initSTG(x)
        use timedata
        use STG_BC
        use pvsQbi
        use turbSA
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
        real*8, allocatable :: rTemp(:),tdump(:)
        real*8, allocatable :: tempdump(:,:)
        double precision, allocatable ::writeTOFile(:,:),thisData(:,:)
        double precision, allocatable ::writeTOFile_a(:,:)
        double precision, allocatable ::nNSurfs(:,:), nNSurfs_a(:,:)
        double precision, allocatable :: tempAllo(:,:)
        double precision :: temp
        real*8, dimension(3) :: tHat,bHat 
        real*8 :: x(numnp,nsd)
        real*8  :: theta
        real*8 :: hxS1
        real*8 :: hxS2
        real*8 :: hyS1
        real*8 :: hyS2
        real*8 :: hzS1
        real*8 :: hzS2
        real*8 :: norm

        logical exlog
        integer,allocatable :: sta(:),holdSurf(:),markHold(:)
        integer :: stopflag,pp,ii,jj,kk,prs,newtype


        leMax=zero
        kMax=zero
        deltaBLSTG=STGDelBL
        alphaGR=STGMeshGrow    ! growth rate of boundary layer mesh
! Initialize Fourier Modes for STG        
        !begin by finding nKWave
        !Read in STGinflow 

c       If a file with the baseline profile at the inflow exists, read from it and set a flag for later
        inquire(file="STGInflow.dat",exist=exlog)
        if(exlog) then ! the STG baseline profile is set in a file so read it
            open (unit=654,file="STGInflow.dat",status="old")
            read(654,*) nPoints ! first line of file MUST contain number of rows of data
c           This file contains in the columns d2wal, ubar(1:3), Reynolds stresses, eddy viscosity,k,w
c           Where Re stresses are R11,R22,R33,R12,R13,R23
            allocate(STGInflow(nPoints,12))
            do i=1,nPoints
              read(654,*) (STGInflow(i,j),j=1,12) ! read the data
            enddo
            iSTGInread=1
            if (myrank.eq.master) then
               write(*,*) "Reading from STGInflow.dat"
            endif
        else ! the file does not exist so set flag to zero and call findUBarandDeriv later 
            iSTGInread=0
        endif
        !Find the spacing
        STGmeth=0;
        alphaGR=STGMeshGrow    ! growth rate of boundary layer mesh
        call findSTG_mesh_hs(x)
        !loop over STG Plane find lcut->kCut->kMax
        allocate(le(nNSurf),lt(nNSurf),kCut(nNSurf))
        do n=1,nNSurf
            lCut=2.0*min(max(hy1(n),hz1(n),0.3*hMax(n))
     &             +0.1*d2Wall(stgSurf(n)),hMax(n))
            kCut(n)=2.0*atan2(0.0,-1.0)/lCut
            if(1.5*kCut(n).gt.kMax)then
                kMax=1.5*kCut(n)
            endif
            !Find le(:) and leMax
            
                y=d2Wall(stgSurf(n))
                call getLt(x,n) !lt(n)=3000.00
                if(2*y.le.3.0*lt(n))then
                    le(n)=2*y
                else
                    le(n)=3.0*lt(n)
                endif
                if(le(n).gt.leMax.and.y.le.deltaBLSTG)then
                    leMax=le(n)
                endif
        enddo
        !communicate kMax about proccessors
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(kMax,kMax_all,1,MPI_REAL8,MPI_MAX,
     &   MPI_COMM_WORLD,ierr)
        kMax=kMax_all
       !communicate leMax about Proccessors
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(leMax,leMax_all,1,MPI_REAL8,MPI_MAX,
     &   MPI_COMM_WORLD,ierr)
        leMax=leMax_all
        !Find kMin
        kmin= atan2(0.0,-1.0)/leMax
        !read in alpha for mode growth
        alphaWave=STGModeGrow
        !Find nKWave 
        nKWave=ceiling((log(kMax/kMin)/log(1+alphaWave))+1)
        if(nKWave.lt.0)then
            print *, "This Mesh is too course on the STG inflow plane 
     &       and PHASTA will now exit"
            call exit       
        endif
        !find other variables
        uNaught=STGUo
        kappa=0.41
        b=5.2
        nu=datmat(1,2,1) !1.81e-5  ! this should come from material properties
        reTau=400.0
        uTau=nu*reTau/deltaBLSTG
        !find nKWave
        
        allocate(rTemp(nKWave*5))
        allocate (dVect(nKWave,3),sigVect(nKWave,3),phiVect(nKWave))

c       Compute the random variables only on the first step
        if (lstep.eq.iSTGStart) then
          call random_seed()
          call random_number(rTemp)
          do id=1,nKWave
              dVect(id,1)=rTemp((id-1)*5+1)*2-1
              dVect(id,2)=rTemp((id-1)*5+2)*2-1
              dVect(id,3)=rTemp((id-1)*5+3)*2-1
              norm=sqrt(dVect(id,1)**2+dVect(id,2)**2+dVect(id,3)**2)
              dVect(id,1:3)=dVect(id,1:3)/norm
              
              phiVect(id)=rTemp((id-1)*5+4)*(2*atan2(0.0,-1.0))
              theta=rTemp((id-1)*5+5)*(2*atan2(0.0,-1.0))            
              if (maxVal(dVect(id,1:3)).eq.dVect(id,1))then
                  tHat(1)=-(dVect(id,2)+dVect(id,3))/dVect(id,1)
                  tHat(2)=1
                  tHat(3)=1
              elseif (maxVal(dVect(id,1:3)).eq.dVect(id,2))then
                  tHat(1)=1
                  tHat(2)=-(dVect(id,1)+dVect(id,3))/dVect(id,2)
                  tHat(3)=1
              else 
                  tHat(1)=1
                  tHat(2)=1
                  tHat(3)=-(dVect(id,2)+dVect(id,1))/dVect(id,3)
              endif
              norm=sqrt(tHat(1)**2+tHat(2)**2+tHat(3)**2)
              tHat(1:3)=tHat(1:3)/norm

              bHat(1)=dVect(id,2)*tHat(3)-dVect(id,3)*tHat(2)
              bHat(2)=dVect(id,3)*tHat(1)-dVect(id,1)*tHat(3)
              bHat(3)=dVect(id,1)*tHat(2)-dVect(id,2)*tHat(1)  
              norm=sqrt(bHat(1)**2+bHat(2)**2+bHat(3)**2)
              bHat(1:3)=bHat(1:3)/norm
            
              sigVect(id,1:3)=tHat*cos(theta)+bHat*sin(theta)
  
          enddo
c          Put all needed random numbers in a single array to be written to restart
           allocate(STGrnd(nKWave,7))
           STGrnd(:,1:3)=dVect
           STGrnd(:,4)=phiVect
           STGrnd(:,5:7)=sigVect

c       If not the first time step of STG, the random variables were read from the restart files
        else
           dVect=STGrnd(:,1:3)
           phiVect=STGrnd(:,4)
           sigVect=STGrnd(:,5:7)

        endif

        allocate(uBar(nNsurf,nsd),duBardy(nNsurf),duBardyADJ(nNsurf)) 


        inquire(file='xyzts.dat',exist=exts)
! Check to see if time statistics are wanted and where they are wanted
        if(iSTGspec.eq.1.and.(.not.exts.or.lstep.eq.iSTGStart))then
            if(myRank.eq.master)then
            write(*,*) "Finding Probe Points and Writing to xyzts.dat"
            endif
            allocate(tdump(4))
            do i=1,4
                tdump(i)=0.0
            enddo
            !check to see if xyzts.dat exists and if it does save the points it has in it
           inquire(file='xyzts_add.dat',exist=exts) 
           if(exts.and.myrank.eq.master)then
               open(unit=626,file='xyzts_add.dat',status='old')
               read(626,*) tdump(1),tdump(2),tdump(3),tdump(4)
                allocate(tempdump(int(tdump(1)),3))
               do jj=1,int(tdump(1))
                   read(626,*) tempdump(jj,1),tempdump(jj,2),tempdump(jj,3)
               enddo
               close(626,status='DELETE')
            endif
            inquire(file='xyzts.dat',exist=exts)
            if(exts)then
               open(unit=626,file='xyzts.dat',status='old')
               close(626,status='DELETE')
            endif
            !accumulate the STG surface nodes on master, broadcast them to all 
            
            allocate(sta(MPI_STATUS_SIZE))
            call MPI_COMM_SIZE(MPI_COMM_WORLD,prs,ierr)
            if(myRank.eq.master)then 
                !make an array that will hold all of the amounts of surface points for each process
                allocate(holdSurf(prs))
            endif
            !Halt all proccesses
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            !Gather nNSurf into holdSurf in the master process
            call MPI_GATHER(nNSurf,1,MPI_INT,holdSurf,1,MPI_INT,master
     &          ,MPI_COMM_WORLD,ierr)
            !add up lengths and allocate nNSurfs_a
            nNSurfs_tot=0
            if(myRank.eq.master)then
                do ii=1,prs
                    nNSurfs_tot=nNSurfs_tot+holdSurf(ii)
                enddo
                !Assign the final Gathered array of surface positions
                allocate(nNSurfs_a(nNSurfs_tot,3))
            endif
            !Create a MPI vector type (real8,real8,real8) 
            call MPI_TYPE_CONTIGUOUS(3,MPI_DOUBLE_PRECISION,
     &          newtype,ierr)
            call MPI_TYPE_COMMIT(newtype,ierr)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            !Send from each non-master proccess the position vectors of the nNSurf nodes
            if(myRank.ne.master.and.nNSurf.ne.0) then
                !create an array sending pos vectors of nNSurf nodes
                call MPI_SEND(x(stgSurf,:),nNSurf,
     &          newtype,master,5253,MPI_COMM_WORLD,ierr)
            endif
            
            !Recieve from each proccess to master
            if(myRank.eq.master)then
                cou=1
                i=1
                do while(i.le.prs)
                    if(i-1.eq.master)then
                        nNSurfs_a(cou:cou+nNSurf-1,:)=x(stgSurf,:)
                        cou=cou+nNSurf
                    elseif(holdSurf(i).ne.0)then
                        allocate(thisData(holdSurf(i),3))
                        thisData=0.0
                        call MPI_RECV(thisData,holdSurf(i),
     &                       newtype,i-1,5253,
     &                      MPI_COMM_WORLD,sta,ierr)
                        nNSurfs_a(cou:cou+holdSurf(i)-1,1:3)=thisData
                        deallocate(thisData)
                        cou=cou+holdSurf(i)
                    endif
                    i=i+1
                enddo
            endif
            !Brodcast nNsurfs_a to all proccesses
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call MPI_BCAST(nNSurfs_tot,1,MPI_INTEGER,master,
     &          MPI_COMM_WORLD,ierr)
            if(myRank.ne.master)then
                allocate(nNSurfs_a(nNSurfs_tot,3))
            endif
            call MPI_BCAST(nNSurfs_a,nNSurfs_tot,newtype,master,
     &          MPI_COMM_WORLD,ierr)

            !FIND POINTS ON this PROCCESS TO WRITE TO XYZTS.DAT

              next=1
              do i=1,nshg!local accumulation
                  if(d2Wall(i).gt.1e-15.and.
     &                d2Wall(i).lt.deltaBLSTG.and.
     &                size(writeToFile_a,DIM=1).le.iSTGEngPNTS)then
                      mark=0
                      NNcheck=1
                      do while(NNcheck.le.nNSurfs_tot.and.mark.ne.1)
                          jj=1
                          do while(jj.le.size(STGdes).and.mark.ne.1)!loop over all input downstream distances
                              if(abs(nNSurfs_a(NNcheck,1)+STGdes(jj)
     &                            -x(i,1)).lt.STGdesol)then
                                  mark=1
                                  !reAssignment of writeToFile_a
                                  if(allocated(writeToFile_a))then
                                      next=next+1
                                      deallocate(writeToFile_a)
                                  endif
                                  allocate(writeToFile_a(next,3))
                                  if(allocated(tempAllo))then
                                      writeToFile_a(1:next-1,:)=tempAllo
                                  endif
                                  writeToFile_a(next,:)=x(i,:)
                                  if(allocated(tempAllo))then
                                      deallocate(tempAllo)
                                  endif
                                  allocate(tempAllo(next,3))
                                  tempAllo=writeToFile_a
                              endif
                              jj=jj+1
                          enddo
                          NNcheck=NNcheck+1
                      enddo
                  endif
              enddo

            stopFlag=0
                if(myRank.eq.master)then
                    do ip=0,prs-1
                        if(ip.eq.master)then
                            StopFlag=stopFlag+size(writeToFile_a,DIM=1)
                        else
                            NNCheck=0
                            call MPI_RECV(NNCheck,1,MPI_INTEGER,ip,
     &                          5111,MPI_COMM_WORLD,sta,ierr)
                            if(recvTot.ne.-1)then
                                stopFlag=stopFlag+NNCheck
                            endif
                        endif
                        kk=kk+1    
                    enddo
                else
                    if(.not.allocated(writeToFile_a))then
                        call MPI_SEND(-1,1,
     &                  MPI_INTEGER,master,5111,MPI_COMM_WORLD,ierr)
                    else
                        NNCheck=size(writeToFile_a,DIM=1)
                        call MPI_SEND(NNCheck,1,
     &                  MPI_INTEGER,master,5111,MPI_COMM_WORLD,ierr)
                    endif
                endif

              call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            jj=0
            do ip=0,prs-1!loop over proccesses number
                if(myrank.eq.master)then
                    write(*,*) "Writing Probe Points to xyzts.dat
     &              from proccess: ",ip
                endif
              !write to xyzts.dat one by one
              if(ip.ne.myRank)then
                  call MPI_BCAST(jj,1,MPI_INTEGER,ip,MPI_COMM_WORLD,
     &            ierr)
              else

            !write to the xyzts.dat file for ip proccess
                open (unit=626,file="xyzts.dat",position="append")
                if(ip.eq.0)then
                    if(stopFlag.gt.iSTGEngPNTS)then
                        write(626,*) iSTGEngPNTS,1,1,1,1
                    else
                        write(626,*) stopFlag,1,1,1,1
                    endif
                endif
                if(.not.allocated(writeTOFile_a))then
                    call MPI_BCAST(jj,1,MPI_INTEGER,ip,MPI_COMM_WORLD,
     &              ierr)
                    close(626)
                else
                    ii=1
                    do while(ii.le.size(writeTOFile_a,dim=1)
     &                  .and.jj.lt.iSTGEngPNTS )
                        write(626,'(F18.15,F18.15,F18.15)')
     &  writeTOFile_a(ii,1), writeTOFile_a(ii,2), writeTOFile_a(ii,3)
                        ii=ii+1
                        jj=jj+1
                    enddo
                    close(626)
                    call MPI_BCAST(jj,1,MPI_INTEGER,ip,MPI_COMM_WORLD,
     &              ierr)
                endif
              endif
            enddo


            !xyzts.dat file must exist before other proccesses hit ln 687 "c...DUMP TIME SERIES"
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

               
                   
        endif 
        end
!----------------------------------------------------------------------
!
!     Apply STG at current time step (modeled after BCprofile not bctint.f)
!
!-----------------------------------------------------------------------
      subroutine applySTG(t, BC,x)
      use STG_BC
      use turbSA ! d2wall(1:numnp) gives distance to wall for each node
      use pvsQbi ! ndsurf(1:nshg) comes from this module and it has "marked" the inflow nodes for this BC
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
      real*8 ::  x(numnp,nsd)
      real*8 ::  BC(nshg,ndofBC)
      
      real*8 :: ke
      real*8 :: midTemp,hTemp,lTemp,uBar_Tol,ka,qN,totqN,rPrime1,rPrime2 
      real*8 :: rPrime3,dNdotRPrime,vPrime1,vPrime2,vPrime3
      real*8 :: turbKin,lAta
      real*8 :: fcut,leta,keta,feta
      real*8, dimension(3,3) :: Cho
      real*8 :: eng(nNSurf,nKWave),totEng(nNSurf)
      real*8 :: kN(nKWave)
      real*8 :: uPrime(nNSurf,3)
      real*8 :: reyS(6), turbVisc, dist, sumBC, tmp
      integer :: noTurbMod
      

       noTurbMod=iRANS+iLES ! 0 only if No-Model is selected (DNS)
       engTot=zero
       uBar_Tol=0.01
        !This statement clears then recreates two files pertaining to 
        !STG inflow plane coordinates and current U output
c        call openFiles()

       uBar=zero
   
       leta=(nu**3/STGeps)**0.25 ! Kolmogorov length scale
       keta=2.0*atan2(0.0,-1.0)/leta ! Kolmogorov wave number
        do n=1,nNSurf
            totEng(n)=0.0
            !This is for the energySpectrum for the n Point
             do k=1,nKWave
                 kN(k)=kmin*((1+alphaWave)**(k-1))
            enddo
            do k=1,nKWave
                  lAta=atan2(0.0,-1.0)
                  ke=2.0*atan2(0.0,-1.0)/le(n)
                  eng(n,k)=(kN(k)/ke)**4.0/((1.0+2.4*(kN(k)/ke)**2.0)**(17.0/6.0))
!                  eng(n,kN)=eng(n,kN)*exp(-((12*kN*lAta/(2*atan2(0.0,-1.0)))**2))
                  feta=exp(-(12*kN(k)/keta)**2) ! damping function
                  fcut=exp(-(((4.0*max(kN(k)-0.9*kCut(n),zero))/kCut(n))**3.0))
!                  write(*,*) d2Wall(stgSurf(n)), lcut,  kCut
                  eng(n,k)=eng(n,k)*fcut*feta
                  if(k.eq.1)then
                      totEng(n)=totEng(n)+(eng(n,k)*kN(k))
                  else
                      totEng(n)=totEng(n)+(eng(n,k)*(kN(k)-kN(k-1)))
                  endif  
!                  if(stgSurf(n).eq.96504.or.stgSurf(n).eq.314549)then       !Midplane energyspecra at n=314549. energy spectra close to wall at n=96504. 
!                    if(istep.eq.0)    then
!                         write(121,*) kN(k)
!                         write(120,*) eng(n,k)
!                    endif
!                  endif
              enddo
            
        enddo 
        do n=1,nNSurf
            if(iLES.ne.0.or.noTurbMod.eq.0) then
                call findUBarandDeriv(x,n,reyS,turbVisc)
                dist=d2Wall(stgSurf(n))
                if(abs(dist).gt.1e-15.and.
     &              dist.lt.1.2*deltaBLSTG)then
                !Find Reynold's Stress/cholesky decomp
                    if (reyS(1).gt.zero) then
                       Cho(1,1)=sqrt(reyS(1))
                    else
                       write(*,*) 'WARNING: reyS(1) <= 0 at y= ',dist
                       Cho(1,1)=sqrt(1.0e-3)
                    endif
                    Cho(1,2)=0
                    Cho(1,3)=0
                    Cho(2,1)=reyS(4)/Cho(1,1)
                    tmp=reyS(2)-((Cho(2,1))**2)
                    if (tmp.gt.zero) then
                       Cho(2,2)=sqrt(tmp)
                    else
                       write(*,*) 'WARNING: reyS(2)>=Cho(2,1)^2 at y= ',
     &                             dist
                       !Cho(2,2)=sqrt(1e-3)
                       Cho(2,2)=-one*1.150*Cho(2,1)
                    endif
                    Cho(2,3)=0
                    Cho(3,1)=reyS(5)/Cho(1,1)
                    Cho(3,2)=(reyS(6)-Cho(2,1)*Cho(3,1))/Cho(2,2)
                    tmp=reyS(3)-(Cho(3,1)**2)-(Cho(3,2)**2)
                    if (tmp.gt.zero) then
                       Cho(3,3)=sqrt(tmp)
                    else
                       write(*,*) 'WARNING: problem with Cho(3,3) at
     &                              y=',dist
                       Cho(3,3)=sqrt(1.0e-3)
                    endif
                    !assign qN,r',and v'
                    vPrime1=0.0
                    vPrime2=0.0
                    vPrime3=0.0
                    totqN=0.0
                    do k=1,nKWave
                        if(totEng(n).eq.0)then
                            qN=0.0
                        elseif(k.eq.1)then
                            qN=(eng(n,k)*kN(k))/totEng(n) 
                            totqN=qN
                        else
                            qN=(eng(n,k)*(kN(k)-kN(k-1)))/totEng(n)
                            totqN=totqN+qN
                        endif
                        rPrime1=(2.0*atan2(0.0,-1.0)*(x(stgSurf(n),1)-uNaught*t))/(kN(k)*leMax)
                        rPrime2=x(stgSurf(n),2)
                        rPrime3=x(stgSurf(n),3)
                        dNdotRPrime=dVect(k,1)*rPrime1+dVect(k,2)*rPrime2+dVect(k,3)*rPrime3
                        vPrime1=vPrime1+sqrt(qN)*(sigVect(k,1)*cos(kN(k)*dNdotRPrime+phiVect(k)))
                        vPrime2=vPrime2+sqrt(qN)*(sigVect(k,2)*cos(kN(k)*dNdotRPrime+phiVect(k)))
                        vPrime3=vPrime3+sqrt(qN)*(sigVect(k,3)*cos(kN(k)*dNdotRPrime+phiVect(k)))
                    enddo
                    vPrime1=2.0*sqrt(1.5)*vPrime1
                    vPrime2=2.0*sqrt(1.5)*vPrime2
                    vPrime3=2.0*sqrt(1.5)*vPrime3

                    uPrime(n,1)=vPrime1*Cho(1,1)+vPrime2*Cho(1,2)+vPrime3*Cho(1,3)
                    uPrime(n,2)=vPrime1*Cho(2,1)+vPrime2*Cho(2,2)+vPrime3*Cho(2,3)
                    uPrime(n,3)=vPrime1*Cho(3,1)+vPrime2*Cho(3,2)+vPrime3*Cho(3,3)
                else ! node is outside BL
                    uPrime(n,1)=0.0
                    uPrime(n,2)=0.0
                    uPrime(n,3)=0.0
                endif
                sumBC=BC(stgSurf(n),3)+BC(stgSurf(n),4)+BC(stgSurf(n),5)
                if(sumBC.ne.zero) then ! protect no slip BC from getting overwritten
                   BC(stgSurf(n),3)=uBar(n,1)+uPrime(n,1)            
                   BC(stgSurf(n),4)=uBar(n,2)+uPrime(n,2)            
                   BC(stgSurf(n),5)=uBar(n,3)+uPrime(n,3)   
                endif         
            else ! iLES.ne.0
!   if the model is not scale-resolving (RANS) then don't include u' but include eddy visc
                call findUBarandDeriv(x,n,reyS,turbVisc)
                sumBC=BC(stgSurf(n),3)+BC(stgSurf(n),4)+BC(stgSurf(n),5)
                if(sumBC.ne.zero) then ! protect no slip BC from getting overwritten
                   BC(stgSurf(n),3)=uBar(n,1)            
                   BC(stgSurf(n),4)=uBar(n,2)
                   BC(stgSurf(n),5)=uBar(n,3)  
                   BC(stgSurf(n),7)=turbVisc
                endif
            endif ! end if iLES ne 0              
!                if(istep.eq.0) write(122,*) BC(stgSurf(n),3),BC(stgSurf(n),4), BC(stgSurf(n),5) 
        enddo 
!End U calculation over the plane at this Time step

!Close varification Files
!        call closeFiles()


!        Loop to check if just put NaNs or Infs in the BC array
!         do n=1,nNSurf
!             do j=3,5
!               if(BC(stgsurf(n),j).ne.BC(stgsurf(n),j)) then
!                 write (*,*) 'Found a NaN in BC array, clm ', j
!               endif
!               if(abs(BC(stgsurf(n),j)).gt.1e3) then
!                 write(*,*) 'Found a big number in BC array, clm ',j
!               endif
!             enddo
!         enddo
 
       end


    
        subroutine openFiles()
       open(unit=117,file='coords4Matlab.dat',status='unknown',iostat=ierr)
       if (ierr.eq.0)then
            close(unit=117,status='delete')
         endif
         open(unit=117,file='coords4Matlab.dat')
 
        open(unit=118,file='currentUbar.dat',status='unknown',iostat=ierr)
        if(ierr.eq.0)then
            close(unit=118,status='delete')
         endif
         open(unit=118,file='currentUbar.dat')

         open(unit=119,file='duBardY.dat',status='unknown',iostat=ierr)
         if(ierr.eq.0)then
             close(unit=119,status='delete')
        endif
        open(unit=119,file='duBardY.dat')

        open(unit=120,file='energySpec.dat',status='unknown',iostat=ierr)
        if(ierr.eq.0)then
            close(unit=120,status='delete')
        endif
        open(unit=120,file='energySpec.dat')
         open(unit=121,file='kN.dat',status='unknown',iostat=ierr)
         if(ierr.eq.0)then
             close(unit=121,status='delete')
         endif
         open(unit=121,file='kN.dat')
          open(unit=122,file='uTot.dat',status='unknown',iostat=ierr)
          if(ierr.eq.0)then
              close(unit=122,status='delete')
          endif
          open(unit=122,file='uTot.dat')

        end

        subroutine closeFiles()
        close(unit=117)
        close(unit=118)
        close(unit=119)
        close(unit=120)
        close(unit=121)
        close(unit=122)
        end



      subroutine findSTG_mesh_hs(x)

      use turbsa
      use pvsQbi
      use STG_BC
      use pointer_data
      include "common.h"
      include "mpif.h"
      integer, allocatable :: ienb(:)
      real*8 x(numnp,3)
      real*8, allocatable :: d2wl(:),xl(:),yl(:),zl(:)
      integer, allocatable :: above(:),below(:)
      real*8 :: hxS1
      real*8 :: hxS2
      real*8 :: hyS1
      real*8 :: hyS2
      real*8 :: hzS1
      real*8 :: hzS2
    
!--------------------------------------------------------------------------------------------------------------------!
!
!                                           Finding the Spacing via elements  
!
!--------------------------------------------------------------------------------------------------------------------! 
      if(STGmeth.eq.0)then
!This is for the mapping that are still required by both methodologies so the the apply STG may loop over the survace rather than the whole mesh
        nNSurf=0
        do n=1,nshg
            if(ndsurf(n).eq.iSTGsurfID)then
                nNSurf=nNSurf+1
                
            endif
        enddo
        allocate(stgSurf(nNSurf))   ! takes in a surface node number and returns the global node number
        allocate(mapBack(numnp))   ! takes in a global node number and returns its surface node number
        nNSurf=0
        do n=1,nshg
            if(ndsurf(n).eq.iSTGsurfID)then
                !fill stgSurf with the node numbers
                nNSurf=nNSurf+1
                stgSurf(nNSurf)=n
                !Find map from numnp node number to stgSurf node number
                mapBack(n)=nNSurf
            endif
        enddo 
         allocate(hx1(nNsurf),hx2(nNsurf),hy1(nNsurf),hy2(nNsurf),hz1(nNsurf),hz2(nNsurf),hMax(nNsurf),
     & above(nNsurf),below(nNsurf))
!  we will be taking mins of hx1 and maxs of hx2 so initialize to extremes
      hx1=1e8
      hy1=1e8
      hz1=1e8
      hx2=-1e8
      hy2=-1e8
      hz2=-1e8
      above=0
      below=0
      do iblk=1,nelblb ! loop element block
         npro=lcblkb(1,iblk+1)-lcblkb(1,iblk)
         nenbl=lcblkb(6,iblk)
         nshl=lcblkb(9,iblk)
         allocate(ienb(nshl))
         allocate(xl(nshl))
         allocate(yl(nshl))
         allocate(zl(nshl))
         allocate(d2wl(nshl))
         do i=1,npro   ! loop element
            iBCB2=miBCB(iblk)%p(i,2)
            if(iBCB2.eq.iSTGSurfID) then ! this is a boundary element with an STG inflow condition
               ienb(1:nshl)=mienb(iblk)%p(i,1:nshl)
               xl(1:nshl)=x(ienb(1:nshl),1)
               yl(1:nshl)=x(ienb(1:nshl),2)
               zl(1:nshl)=x(ienb(1:nshl),3)
               xmax=maxval(xl)
               ymax=maxval(yl)
               zmax=maxval(zl)
               xmin=minval(xl)
               ymin=minval(yl)
               zmin=minval(zl)
               d2wl(1:nshl)=d2wall(ienb(1:nshl))
               d2wallmax=maxval(d2wl)
               d2wallmin=minval(d2wl)
               hx=xmax-xmin
               hy=ymax-ymin
               hz=zmax-zmin
               do j=1,nenbl  ! loop elemental vertex
                 nn=ienb(j)
                 if(ndsurf(nn).eq.iSTGSurfID) then
                  nns=mapback(nn)
                  hx1(nns)=min(hx1(nns),hx)                 
                  hy1(nns)=min(hy1(nns),hy)                 
                  hz1(nns)=min(hz1(nns),hz)                 
                  hx2(nns)=max(hx2(nns),hx)                 
                  hy2(nns)=max(hy2(nns),hy)                 
                  hz2(nns)=max(hz2(nns),hz) 
                  hMax(nns)=max(hx1(nns),max(hy1(nns),hz1(nns)))        
                  hmino2=min(hx,min(hy,hz))/2
                  if(d2wl(j).lt.(d2wallmax-hmino2)) then
                      above(nns) = 1
                  endif
                  if(d2wl(j).gt.(d2wallmin+hmino2)) then
                      below(nns)= 1
                  endif
                 endif
               enddo  ! end loop elemental vertex
           endif !this was an element with a face on the STG plane and we skip others
         enddo ! end loop element
         deallocate(ienb)
         deallocate(d2wl)
         deallocate(xl)
         deallocate(yl)
         deallocate(zl)
      enddo ! end loop element block 
!
! note the above only accounts for elements on this processor. Given the non-global nature of the h arrays, it would be 
! complicated to communicate among the ranks that have nodes on the STG plane to make each inflow node correct on each 
! rank. For now, we will avoid this complexity and assume that, at least within the boundary layer where these length 
! scales are needed, the mesh has been stretched with a geometric growth rate that we will take as a known input 
! parameter. In this way, we can identify and correct for parallel as follows
      
      do i=1,nNSurf
         if(below(i).eq.0) then  ! this surface node did not have elements below it to get a good hy min so scale hymax
            hy1(i)=hy2(i)/alphaGR
         endif
         if(above(i).eq.0) then  ! this surface node did not have elements above it to get a good hy max so scale hymax
            hy2(i)=hy1(i)*alphaGR
         endif
      enddo

! NOTE this assumes for now that BL grows in y direction.  At some point, this needs to be generalizd to 
! normal direction and "in plane" directions.  Sketching that, I suspect that we would want to define hst, hsp, and hn.  
! The STG plane would be defined with a point and a normal to the plane with an assumption that hst aligns with the
! normal to that plane, hn would be a local normal to the surface ( wall normal which could perhaps be found from taking the gradient
! of d2wall which can be computed for each element), and hsp would be found by crossing the plane normal with the local
! wall normal.

      return
!--------------------------------------------------------------------------------------------------------------------!
!
!                                           Finding the Spacing via search (nshg^2 inefficant) 
!
!--------------------------------------------------------------------------------------------------------------------!      
      elseif(STGmeth.eq.1)then
      

!Find the length of the array that holds global nodes on the STG inflow
        nNSurf=0
        do n=1,nshg
            if(ndsurf(n).eq.iSTGsurfID)then
                nNSurf=nNSurf+1
                
            endif
        enddo
        allocate(stgSurf(nNSurf))   ! takes in a surface node number and returns the global node number
        allocate(mapBack(numnp))   ! takes in a global node number and returns its surface node number
!Find the hx(1,2),hy(1,2),hz(1,2),hMax for each node on the STG inflow
!Also n_hxl,n_hxh,n_hyl,n_hyh,n_hzl,n_hzh
!Current code requires rectangular mesh
        hxS1=minval(x(:,1))-1
        hxS2=maxval(x(:,1))+1
        hyS1=minval(x(:,2))-1
        hyS2=maxval(x(:,2))+1
        hzS1=minval(x(:,3))-1
        hzS2=maxval(x(:,3))+1

        nNSurf=0
        do n=1,nshg
            if(ndsurf(n).eq.iSTGsurfID)then
                !fill stgSurf with the node numbers
                nNSurf=nNSurf+1
                stgSurf(nNSurf)=n
                !Find map from numnp node number to stgSurf node number
                mapBack(n)=nNSurf
                !if h_i is unset the defalt is node where it's at
                n_hxl(n)=n
                n_hxh(n)=n
                n_hyl(n)=n
                n_hyh(n)=n
                n_hzl(n)=n
                n_hzh(n)=n

                !initialize to outside of mesh domain
                tempHx1=hxS1
                tempHx2=hxS2
                tempHy1=hyS1
                tempHy2=hyS2
                tempHz1=hzS1
                tempHz2=hzS2
                do id=1,nshg
                    !hx conditional
                    if(n.ne.id.and.
     &              x(n,1).ne.x(id,1).and.
     &              abs(x(n,2)-x(id,2)).lt.1e-15.and.
     &              abs(x(n,3)-x(id,3)).lt.1e-15)then
                        if(x(n,1).ge.x(id,1).and.x(id,1).ge.tempHx1)then
                            tempHx1=x(id,1)
                            n_hxl(n)=id
                        endif
                        if(x(n,1).le.x(id,1).and.x(id,1).le.tempHx2)then
                            tempHx2=x(id,1)
                            n_hxh(n)=id
                        endif
                    endif
                    !hy conditional
                     if(n.ne.id.and.
     &              abs(x(n,1)-x(id,1)).lt.1e-15.and.
     &              x(n,2).ne.x(id,2).and.
     &              abs(x(n,3)-x(id,3)).lt.1e-15)then
                         if(x(n,2).ge.x(id,2).and.x(id,2).ge.tempHy1)then
                             tempHy1=x(id,2)
                             n_hyl(n)=id
                         endif
                         if(x(n,2).le.x(id,2).and.x(id,2).le.tempHy2)then
                             tempHy2=x(id,2)
                             n_hyh(n)=id
                         endif
                     endif
                    !hz conditional
                     if(n.ne.id.and.
     &              abs(x(n,1)-x(id,1)).lt.1e-15.and.
     &              abs(x(n,2)-x(id,2)).lt.1e-15.and.
     &              x(n,3).ne.x(id,3))then
                         if(x(n,3).ge.x(id,3).and.x(id,3).ge.tempHz1)then
                             tempHz1=x(id,3)
                             n_hzl(n)=id
                         endif
                         if(x(n,3).le.x(id,3).and.x(id,3).le.tempHz2)then
                             tempHz2=x(id,3)
                             n_hzh(n)=id
                         endif
                     endif
                enddo
                !Catch if the node is on the edge of the plane
                if(tempHx1.eq.hxS1)then
                    tempHx1=tempHx1+1
                endif
                 if(tempHx2.eq.hxS2)then
                     tempHx2=tempHx2-1
                 endif
                 if(tempHy1.eq.hyS1)then
                     tempHy1=tempHy1+1
                 endif
                 if(tempHy2.eq.hyS2)then
                     tempHy2=tempHy2-1
                 endif
                 if(tempHz1.eq.hzS1)then
                     tempHz1=tempHz1+1
                 endif
                 if(tempHz2.eq.hzS2)then
                     tempHz2=tempHz2+1
                 endif

                hx1(nNsurf)=abs(x(n,1)-tempHx1)
                hx2(nNsurf)=abs(x(n,1)-tempHx2)
                hy1(nNsurf)=abs(x(n,2)-tempHy1)
                hy2(nNsurf)=abs(x(n,2)-tempHy2)
                hz1(nNsurf)=abs(x(n,3)-tempHz1)
                hz2(nNsurf)=abs(x(n,3)-tempHz2)

                tempHx1=min(hx2(nNsurf),hx1(nNsurf))
                tempHy1=min(hy2(nNsurf),hy1(nNsurf))
                tempHz1=min(hz2(nNsurf),hz1(nNsurf))
                if(tempHx1.ge.tempHy1)then
                    if(tempHx1.ge.tempHz1) then
                        hMax(nNsurf)=tempHx1
                    else
                        hMax(nNsurf)=tempHz1
                    endif
                else
                    if(tempHy1.ge.tempHz1)then
                        hMax(nNsurf)=tempHy1
                    else
                        hMax(nNsurf)=tempHz1
                    endif
                endif
            endif
        enddo
        endif
      end 

      subroutine findUBarandDeriv(x,n,reyS,turbVisc)
      use STG_BC
      use turbSA ! d2wall(1:numnp) gives distance to wall for each node
      use pvsQbi ! ndsurf(1:nshg) comes from this module and it has "marked" the inflow nodes for this BC
      include "common.h"
      real*8 ::  x(numnp,nsd)
      real*8 :: reyS(6), lm, turbVisc
 
      if(STGmeth.eq.0) then ! this is the parallel way
      
           if (iSTGInread.eq.0) then ! this is the case where no profile from a file exists
              utauonu=utau/nu
              utausqonu=utau*utauonu
              rkapinv=one/kappa

              y=d2wall(stgSurf(n))
              if(y.lt.1e-8) then
                   ubar(n,1)=0
                   duBardY(n)=utausqonu
              else
                 yplus=utauonu*y               
                 ulothw=utausqonu*y
                 uloglaw=utau*(B+rkapinv*log(yplus))
c              Curve fit from RANS of Channel flow at Ret=400
                 ufit = 0.70367*log(4.73081*y)-0.66824*log(3.36056*y)
                 if(ulothw.lt.ufit) then
                    uBar(n,1)=ulothw
                    duBardY(n)=utausqonu
                 else
c                   uBar(n)=utau*(B+rkapinv*log(yplus))
c                   duBardY(n)=utau*rkapinv/y
                    uBar(n,1)=ufit
                    duBardY(n)=0.70367/y-0.66824/y
                 endif
              endif
              
              lm=kappa*y*(1.0-exp(-((y*reTau)/(deltaBLSTG*26.0))))
              turbVisc=abs(duBardY(n))*lm**2.0
              spre=sqrt(2.0)*abs(duBardY(n))
              turbKin=turbVisc*spre/0.3
              reyS(1)=2.0*turbKin/3.0
              reyS(4)=turbVisc*abs(duBardY(n))
              reyS(5)=0.0
              reyS(2)=2.0*turbKin/3.0
              reyS(6)=0.0
              reyS(3)=2.0*turbKin/3.0 

           elseif (iSTGInread.eq.1) then ! this is the case where the inflow data is read from file
c             Now need to interpolate from the file data for the current node point
              y=d2wall(stgSurf(n))
              iupper=0
              do j=2,nPoints
                 if(STGInflow(j,1).gt.y) then !bound found
                    xi=(y-STGInflow(j-1,1))/
     &                 (STGInflow(j,1)-STGInflow(j-1,1))
                    iupper=j
                    exit
                 endif
              enddo
              if(iupper.eq.0) then ! node is higher than interpolating stack so use the top of interpolating stack
                 iupper=nPoints
                 xi=1.0
              endif
              uBar(n,1:3)=(xi*STGInflow(iupper,2:4)
     &                +(one-xi)*STGInflow(iupper-1,2:4))
              reyS=(xi*STGInflow(iupper,5:10)
     &                +(one-xi)*STGInflow(iupper-1,5:10))
              turbVisc=(xi*STGInflow(iupper,11)
     &                +(one-xi)*STGInflow(iupper-1,11))

!             This is for the channel where the shear stress changes sign across the centerline
              if(iSTGChan.eq.1) then
                if (x(stgSurf(n),2).gt.deltaBLSTG) then
                   reyS(4)=-one*reyS(4)
                endif
              endif

           endif ! whether inflow profile is read from file



      elseif(STGmeth.eq.1)then !Original way (not parallel)
                lTemp=0
                hTemp=50
                do while(hTemp-lTemp.ge.uBar_Tol)
                    midTemp=(hTemp+lTemp)/2
                    if(0.ge.(midTemp+Exp(-kappa*b)*
     &              (Exp(kappa*midTemp)-1-(kappa*midTemp)-(((kappa*midTemp)**2)/2)-(((kappa*midTemp)**3)/6))
     &              -reTau*d2Wall(stgSurf(n))/deltaBLSTG))then
                        lTemp=midTemp
                    elseif(0.eq.(midTemp+Exp(-kappa*b)*
     &              (Exp(kappa*midTemp)-1-(kappa*midTemp)-(((kappa*midTemp)**2)/2)-(((kappa*midTemp)**3)/6))
     &              -reTau*d2Wall(stgSurf(n))/deltaBLSTG))then
                        lTemp=midTemp
                        hTemp=midTemp
                    else
                        hTemp=midTemp
                    endif
                enddo
           !turn u+ to uBar for this node (1st midTemp=u+ then midTemp=ubar)
                uBar(n,1)=midTemp*uTau
                if(uBar(n,1).ge.uNaught) then
                    uNaught=uBar(n,1)
                endif

                duBardy(n)=(uBar(mapBack(n_hyh(stgSurf(n))),1)-uBar(mapBack(n_hyl(stgSurf(n))),1))/
     &           (x(n_hyh(stgSurf(n)),2)-x(n_hyl(stgSurf(n)),2))
                !for lm hypothesis duBardy should always be increasing away from a wall
                !and so the centerline for channel flow must be fixed. Use duBardyADJ at a node away
                !for this mesh
!                if(istep.eq.0)         write(119,*) duBardy
                    if(x(n,2)-d2Wall(n).eq.0)then
                        duBardyADJ(n)=
     &                  (uBar(mapBack(n_hyl(n_hyl(stgSurf(n)))),1)
     &                  -uBar(n,1))/
     &                  (x(n_hyl(n_hyl(stgSurf(n))),2)-
     &                  x(stgSurf(n),2))
                    else
                        duBardyADJ(n)=
     &                  (uBar(mapBack(n_hyh(n_hyh(stgSurf(n)))),1)-
     &                  uBar(n,1))/
     &                  (x(n_hyh(n_hyh(stgSurf(n))),2)-
     &                  x(stgSurf(n),2))
                    endif
                if(abs(duBardy(n)).lt.1e-2)then
                    duBardy(n)=duBardyADJ(n)
                endif
             endif !end of original way abandonded for parallel
                   ! if(istep.eq.0)      write(118,*) uBar(n),0,0
        
      end


      subroutine getLt(x,n)
      use STG_BC
      use turbSA ! d2wall(1:numnp) gives distance to wall for each node
      use pvsQbi ! ndsurf(1:nshg) comes from this module and it has "marked" the inflow nodes for this BC
      include "common.h"
      
      real*8 ::  x(numnp,nsd)

c      Now need to interpolate from the file data for the current node point
       y=d2wall(stgSurf(n))
       iupper=0
       do j=2,nPoints
           if(STGInflow(j,1).gt.y) then !bound found
               xi=(y-STGInflow(j-1,1))/
     &             (STGInflow(j,1)-STGInflow(j-1,1))
               iupper=j
               exit
           endif
       enddo
       if(iupper.eq.0) then ! node is higher than interpolating stack so use the top of interpolating stack
            iupper=nPoints
            xi=1.0
       endif
       lt(n)=(xi*STGInflow(iupper,12)
     &          +(one-xi)*STGInflow(iupper-1,12))

      end


