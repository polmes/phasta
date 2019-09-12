!-----------------------------------------------------------------------
!
!  This module handles synthetic turbulence generation at inflows and
!   later at LES switchover planes for ADES
!
!-----------------------------------------------------------------------
      module STG_BC

       real*8, allocatable :: dVect(:,:), sigVect(:,:)
       real*8, allocatable :: STGrnd(:,:), BCrestart(:,:) 
       real*8, allocatable :: phiVect(:), phiVectTWO(:)
       real*8, allocatable :: hx1(:),hx2(:),hy1(:),hy2(:)
       real*8, allocatable :: hz1(:),hz2(:),hMax(:)
       real*8, allocatable :: uBar(:,:), dUBardY(:), duBardYADJ(:)
       real*8, allocatable :: STGInflow(:,:)
       integer,allocatable :: n_hxl(:),n_hxh(:),n_hyl(:),n_hyh(:)
       integer,allocatable :: n_hzl(:),n_hzh(:)
       integer,allocatable :: stgSurf(:),mapBack(:) 
      
       integer STGmeth, nKWave,nNSurf
       real*8 alphaWave,alphaGR 
       real*8, allocatable :: le(:),lt(:), kCut(:), keta(:)
       real*8 kMax,kMax_all,lCut,leMax,leMax_all,kMin
       integer iSTGInread, nPoints

!Flow Parameters
        real*8 deltaBLSTG,uNaught,nu
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
        real*8, allocatable :: rTemp(:)
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
        real*8 eps
        real*8 leta, y

        logical exlog
        integer,allocatable :: sta(:),holdSurf(:),markHold(:)
        integer :: stopflag,pp,ii,jj,kk,prs,newtype


        leMax=zero
        kMax=zero
        deltaBLSTG=STGDelBL
        alphaGR=STGMeshGrow    ! growth rate of boundary layer mesh
        nu = datmat(1,2,1) 
! Initialize Fourier Modes for STG        
        !begin by finding nKWave
        !Read in STGinflow 

c       If a file with the baseline profile at the inflow exists, read from it and set a flag for later
        inquire(file="STGInflow.dat",exist=exlog)
        if(exlog) then ! the STG baseline profile is set in a file so read it
            open (unit=654,file="STGInflow.dat",status="old")
            read(654,*) nPoints ! first line of file MUST contain number of rows of data
c           This file contains in the columns d2wal, ubar(1:3), Reynolds stresses, scalar 1-2,lt,epsilon
c           Where Re stresses are R11,R22,R33,R12,R13,R23
            allocate(STGInflow(nPoints,14))
            do i=1,nPoints
              read(654,*) (STGInflow(i,j),j=1,14) ! read the data
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
        allocate(le(nNSurf),lt(nNSurf),kCut(nNSurf),keta(nNSurf))
        do n=1,nNSurf
            lCut=2.0*min(max(hy1(n),hz1(n),0.3*hMax(n))
     &             +0.1*d2Wall(stgSurf(n)),hMax(n))
            kCut(n)=2.0*atan2(0.0,-1.0)/lCut
            call getEps(x,n,eps)
            leta = (nu**3/eps)**0.25
            keta(n) = 2.0*atan2(0.0,-1.0)/leta
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
!        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(kMax,kMax_all,1,MPI_REAL8,MPI_MAX,
     &   MPI_COMM_WORLD,ierr)
        kMax=kMax_all
       !communicate leMax about Proccessors
!        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
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
        
        allocate (rTemp(nKWave*5))
        allocate (dVect(nKWave,3),sigVect(nKWave,3))
        allocate (phiVect(nKWave),phiVectTWO(nKWave))

c       Compute the random variables only on the first step
        if (lstep.eq.iSTGStart) then
          if (myrank.eq.master) write(*,*) 
     &                           'Computing the STG random numbers'
          call random_seed()
          call random_number(rTemp)
          do id=1,nKWave
c               Latest version from John
                theta=acos(-2.0*rTemp((id-1)*5+1)+1.0)
                phiVect(id)=rTemp((id-1)*5+2)*(2.0*atan2(0.0,-1.0))
                phiVectTWO(id)=rTemp((id-1)*5+3)*(2.0*atan2(0.0,-1.0))
                dVect(id,1)=cos(theta)
                dVect(id,2)=sin(theta)*sin(phiVect(id))
                dVect(id,3)=sin(theta)*cos(phiVect(id))
                tHat(1)=cos(theta +atan2(0.0,-1.0)/2)
                tHat(2)=sin(theta+atan2(0.0,-1.0)/2)*sin(phiVect(id))
                tHat(3)=sin(theta+atan2(0.0,-1.0)/2)*cos(phiVect(id))
                bHat(1)=dVect(id,2)*tHat(3)-dVect(id,3)*tHat(2)
                bHat(2)=dVect(id,3)*tHat(1)-dVect(id,1)*tHat(3)
                bHat(3)=dVect(id,1)*tHat(2)-dVect(id,2)*tHat(1)  
                sigVect(id,1:3)=cos(phiVectTWO(id))*bHat(1:3)+
     &                          sin(phiVectTWO(id))*tHat(1:3)
c               Previous version used on flat plate and bump DNS 
!                theta=acos(-2.0*rTemp((id-1)*5+1)+1.0)
!                phiVect(id)=rTemp((id-1)*5+4)*(2.0*atan2(0.0,-1.0))
!                dVect(id,1)=cos(theta)
!                dVect(id,2)=sin(theta)*sin(phiVect(id))
!                dVect(id,3)=sin(theta)*cos(phiVect(id))              
!                theta=acos(-2.0*rTemp((id-1)*5+2)+1.0)
!                phiVect(id)=rTemp((id-1)*5+3)*(2.0*atan2(0.0,-1.0))
!                tHat(1)=cos(theta)
!                tHat(2)=sin(theta)*sin(phiVect(id))
!                tHat(3)=sin(theta)*cos(phiVect(id))
!                bHat(1)=dVect(id,2)*tHat(3)-dVect(id,3)*tHat(2)
!                bHat(2)=dVect(id,3)*tHat(1)-dVect(id,1)*tHat(3)
!                bHat(3)=dVect(id,1)*tHat(2)-dVect(id,2)*tHat(1)  
!                norm=sqrt(bHat(1)**2+bHat(2)**2+bHat(3)**2)
!                sigVect(id,1:3)=bHat(1:3)/norm
          enddo
c          Put all needed random numbers in a single array to be written to restart
           if(allocated(STGrnd)) then 
              STGrnd(:,1:3)=dVect
              STGrnd(:,4)=phiVect
              STGrnd(:,5:7)=sigVect
           else
              allocate(STGrnd(nKWave,7))
              STGrnd(:,1:3)=dVect
              STGrnd(:,4)=phiVect
              STGrnd(:,5:7)=sigVect
           endif

c       If not the first time step of STG, the random variables were read from the restart files
        else
           dVect=STGrnd(:,1:3)
           phiVect=STGrnd(:,4)
           sigVect=STGrnd(:,5:7)

        endif

        allocate(uBar(nNsurf,nsd),duBardy(nNsurf),duBardyADJ(nNsurf)) 


        end
!----------------------------------------------------------------------
!
!     Apply STG at current time step (modeled after BCprofile not bctint.f)
!
!-----------------------------------------------------------------------
      subroutine applySTG(t, BC,x,GradV)
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
      real*8 :: rPrime3,dNdotRPrime,vPrime1,vPrime2,vPrime3,rPrime1Den
      real*8 :: turbKin,lAta
      real*8 :: fcut,feta
      real*8, dimension(3,3) :: Cho
      real*8 :: eng(nNSurf,nKWave),totEng(nNSurf)
      real*8 :: kN(nKWave)
      real*8 :: uPrime(nNSurf,3)
      real*8 :: reyS(6), sclr1, sclr2, dist, sumBC, tmp, eps
      real*8 :: ss1, ss2, ss, ssq, arg1, arg2, tmp1, tmp2, nut,
     &          GradV(nshg,nsdsq), yp
      integer :: inode
      integer :: noTurbMod
      

       noTurbMod=iRANS+iLES ! 0 only if No-Model is selected (DNS)
       engTot=zero
       uBar_Tol=0.01
        !This statement clears then recreates two files pertaining to 
        !STG inflow plane coordinates and current U output
c        call openFiles()

       uBar=zero
   
       !leta=(nu**3/STGeps)**0.25 ! Kolmogorov length scale
       !keta=2.0*atan2(0.0,-1.0)/leta ! Kolmogorov wave number
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
                  feta=exp(-(12*kN(k)/keta(n))**2) ! damping function
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
                call findUBarandDeriv(x,n,reyS,sclr1,sclr2)
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
                        rPrime1Den = kN(k)*min(leMax,10.0*2*atan2(0.0,-1.0)/kN(k))
                        rPrime1=(2.0*atan2(0.0,-1.0)*(x(stgSurf(n),1)-uNaught*t))/rPrime1Den
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

                    if (iRANS.eq.-5.and.iSTGSSTMeth.eq.2) then
                      inode = stgSurf(n)
                      ss1 = GradV(inode,1) ** 2 + GradV(inode,2) ** 2
     &                     + GradV(inode,3) ** 2 + GradV(inode,4) ** 2
     &                     + GradV(inode,5) ** 2 + GradV(inode,6) ** 2
     &                     + GradV(inode,7) ** 2 + GradV(inode,8) ** 2
     &                     + GradV(inode,9) ** 2
                      ss2 = GradV(inode,1) ** 2 + GradV(inode,5) ** 2
     &                     + GradV(inode,9) ** 2 + 2*(
     &                       GradV(inode,2)*GradV(inode,4)
     &                       + GradV(inode,3)*GradV(inode,7)
     &                       + GradV(inode,6)*GradV(inode,8)) 
                      ss = 0.5 * (ss1 + ss2)
                      ssq = sqrt(two*ss)
                      arg1 = 0.410*dist
                      arg2 = 0.3250*hMax(n)
                      tmp1 = min(arg1**2,arg2**2)
                      yp = dist*STGSSTuTau/nu
                      tmp2 = one-exp(-(yp/25)**3)
                      nut = tmp1*tmp2*ssq
                      if (istep.eq.0) then
                         sclr1 = BC(stgSurf(n),7)
                      else
                         sclr1 = nut*sclr2*0.1
                      endif
                    endif
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
                   if (iRANS.eq.-1) BC(stgSurf(n),7)=sclr1
                   if (iRANS.eq.-5) then
                     if (iSTGSSTMeth.eq.1) then
                        BC(stgSurf(n),7)=STGSSTnut*nu*sclr2
                     elseif (iSTGSSTMeth.eq.2) then
                        BC(stgSurf(n),7)=sclr1
                     endif
                     BC(stgSurf(n),8)=sclr2
                   endif
                endif         
            else ! iLES.ne.0
!   if the model is not scale-resolving (RANS) then don't include u' but include eddy visc
                call findUBarandDeriv(x,n,reyS,sclr1,sclr2)
                sumBC=BC(stgSurf(n),3)+BC(stgSurf(n),4)+BC(stgSurf(n),5)
                if(sumBC.ne.zero) then ! protect no slip BC from getting overwritten
                   BC(stgSurf(n),3)=uBar(n,1)            
                   BC(stgSurf(n),4)=uBar(n,2)
                   BC(stgSurf(n),5)=uBar(n,3)  
                   BC(stgSurf(n),7)=sclr1
                   if (iRANS.eq.-5) BC(stgSurf(n),8)=sclr2
                endif
            endif ! end if iLES ne 0              
!                if(istep.eq.0) write(122,*) BC(stgSurf(n),3),BC(stgSurf(n),4), BC(stgSurf(n),5) 
        enddo 
!End U calculation over the plane at this Time step

c        Loop to check if just put NaNs or Infs in the BC array
c         do n=1,nNSurf
c             do j=3,7
c               if(BC(stgsurf(n),j).ne.BC(stgsurf(n),j)) then
c                 write (*,*) 'Found a NaN in BC array, clm ', j
c               endif
c               if(abs(BC(stgsurf(n),j)).gt.1e3) then
c                 write(*,*) 'Found a big number in BC array, clm ',j
c               endif
c             enddo
c         enddo
 
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
      end 

      subroutine findUBarandDeriv(x,n,reyS,sclr1,sclr2)
      use STG_BC
      use turbSA ! d2wall(1:numnp) gives distance to wall for each node
      use pvsQbi ! ndsurf(1:nshg) comes from this module and it has "marked" the inflow nodes for this BC
      include "common.h"
      real*8 ::  x(numnp,nsd)
      real*8 :: reyS(6), lm, sclr1, sclr2
 
      
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
              sclr1=abs(duBardY(n))*lm**2.0
              spre=sqrt(2.0)*abs(duBardY(n))
              turbKin=sclr1*spre/0.3
              reyS(1)=2.0*turbKin/3.0
              reyS(4)=sclr1*abs(duBardY(n))
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
              sclr1=(xi*STGInflow(iupper,11)
     &                +(one-xi)*STGInflow(iupper-1,11))
              sclr2=(xi*STGInflow(iupper,12)
     &                +(one-xi)*STGInflow(iupper-1,12))

!             This is for the channel where the shear stress changes sign across the centerline
              if(iSTGChan.eq.1) then
                if (x(stgSurf(n),2).gt.deltaBLSTG) then
                   reyS(4)=-one*reyS(4)
                endif
              endif

           endif ! whether inflow profile is read from file
       
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
       lt(n)=(xi*STGInflow(iupper,13)
     &          +(one-xi)*STGInflow(iupper-1,13))

      end


      subroutine getEps(x,n,eps)
      use STG_BC
      use turbSA ! d2wall(1:numnp) gives distance to wall for each node
      use pvsQbi ! ndsurf(1:nshg) comes from this module and it has "marked" the inflow nodes for this BC
      include "common.h"
      
      real*8 ::  x(numnp,nsd), eps
      integer n

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
       eps=(xi*STGInflow(iupper,14)
     &          +(one-xi)*STGInflow(iupper-1,14))

      end


