        subroutine proces
c
c----------------------------------------------------------------------
c
c This subroutine generates the problem data and calls the solution
c  driver.
c
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        use readarrays          ! used to access x, iper, ilwork
        use turbsa          ! used to access d2wall
        use dtnmod
        use periodicity
        use pvsQbi
        use mlirans
        use spanStats
        include "common.h"
        include "mpif.h"
c
c arrays in the following 2 lines are now dimensioned in readnblk
c        dimension x(numnp,nsd)
c        dimension iper(nshg), ilwork(nlwork)
c
        dimension y(nshg,ndof),
     &            iBC(nshg),
     &            BC(nshg,ndofBC),
     &            ac(nshg,ndof)
c
c.... shape function declarations
c     
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
c  stuff for dynamic model s.w.avg and wall model
c
c        dimension velbar(nfath,nflow)
c
c stuff to interpolate profiles at inlet
c
        real*8, allocatable :: bcinterp(:,:)
        real*8, allocatable :: IC3d(:,:)
        real*8, allocatable :: ICSTG(:,:)
        real*8 maxdwl,maxx,minx,dwl,z
        integer numSTG
        integer interp_mask(ndof)
        integer Map2ICfromBC(ndof)
        character*256 fname2
        logical exlog

        !Duct
        real*8 c0, c1, c2, c3, x1, x2
        integer nn

c        if ((irscale .ge. 0) .and. (myrank.eq.master)) then
c           call setSPEBC(numnp, point2nfath, nsonmax)
c        endif
c     
c.... generate the geometry and boundary conditions data
c
        call gendat (y,              ac,             point2x,
     &               iBC,            BC,
     &               point2iper,     point2ilwork,   shp,
     &               shgl,           shpb,           shglb,
     &               point2ifath,    point2nsons )
        call setper(nshg)
        call perprep(iBC,point2iper,nshg)
        if (iLES/10 .eq. 1) then
        call keeplhsG ! Allocating space for the mass (Gram) matrix to be
                      ! used for projection filtering and reconstruction
                      ! of the strain rate tensor.

        call setrls   ! Allocating space for the resolved Leonard stresses
c                         See bardmc.f 
        endif
c
c.... time averaged statistics
c
        if (ioform .eq. 2) then
           call initStats(point2x, iBC, point2iper, point2ilwork) 
        endif
c
c.... RANS turbulence model
c
        if (iRANS .lt. 0.or.iSTG.eq.1) then
           call initTurb( point2x )
        endif
c
c.... p vs. Q boundary
c
           call initNABI( point2x, shpb )
c     
c.... check for execution mode
c
        if (iexec .eq. 0) then
           lstep = 0
           call restar ('out ',  y  ,ac)
           return
        endif
c
c.... initialize AutoSponge
c
        if(matflg(5,1).ge.4) then ! cool case (sponge)
           call initSponge( y,point2x) 
        endif
c
c....  initialize the STG inflow arrays
c
        if(iSTG.eq.1) then 
           call initSTG(point2x)
        endif

c
c
c.... adjust BC's to interpolate from file
c
        
        inquire(file="inlet.dat",exist=exlog)
        if(exlog) then

           !display to screen if BC's will be adjusted from file
           if(myrank == 0)
     &         write(*,*)"adjusting BCs or ICs from inlet.dat file"

           !open inlet.dat file for reading
           open (unit=654,file="inlet.dat",status="old")

           !read inlet.dat meta data
           read(654,*) ninterp,ICset,iBCset,isrfidmatch,(interp_mask(j),j=1,ndof)
           !ninterp = number of data points to interpolate between
           !ICset = 0,1 depending on whether IC need to be set
           !iBCset = 0,1 depending on whether BCs need to be set
           !isrfidmatch = surf ID where BCs need to be set
           !interp_mask(j) = 0,1 whether p,T,u,v,w,scalar needs to be set


           if(ICset.eq.1) then  ! we need a map that takes IC# and returns IC #
             Map2ICfromBC(1)=4
             Map2ICfromBC(2)=5
             Map2ICfromBC(3)=1
             Map2ICfromBC(4)=2
             Map2ICfromBC(5)=3
             do j=6,ndof  ! scalars have an identity map
               Map2ICfromBC(j)=j
             enddo
           endif
            
 
           allocate(bcinterp(ninterp,ndof+1))

           !loop through all lines in file to read in dwall,p,T,u,v,w,scalar
           !note: values must be in increasing distance from wall order
           do i=1,ninterp
              read(654,*) (bcinterp(i,j),j=1,ndof+1)
           enddo
           do i=1,nshg  ! only correct for linears at this time
              iBcon = 0
              if((ndsurf(i).eq.isrfidmatch).and.(iBCset.eq.1)) iBCon=1
              if((IBCon.eq.1).or.(ICset.eq.1)) then ! need to compute the bounding points and interpolator, xi
                 iupper=0
                 do j=2,ninterp
                    if(bcinterp(j,1).gt.d2wall(i)) then !bound found
                       xi=(d2wall(i)-bcinterp(j-1,1))/
     &                    (bcinterp(j,1)-bcinterp(j-1,1))
                       iupper=j
                       exit
                    endif
                 enddo
                 if(iupper.eq.0) then ! node is higher than interpolating stack
! so use the top of interpolating stack
                    iupper=ninterp
                    xi=1.0
                 endif
                 if(iBCon.eq.1) then
                   do j=1,nflow
                     if(interp_mask(j).eq.1) then 
                        BC(i,j)=(xi*bcinterp(iupper,j+1)
     &                    +(one-xi)*bcinterp(iupper-1,j+1))
                     endif
                   enddo
                   if(solheat.eq.1) 
     &                  BC(i,2)=(i*bcinterp(iupper,3)
     &                    +(one-xi)*bcinterp(iupper-1,3))
                   do j=1,nsclr
                     if(interp_mask(j+5).eq.1) then 
                        BC(i,j+6)=(xi*bcinterp(iupper,j+6)
     &                    +(one-xi)*bcinterp(iupper-1,j+6))
                     endif
                   enddo
                 endif ! IBCon
                 if(ICset.eq.1) then  ! note at this time no mask so all IC's are getting set
! we could add a separate mask for IC or potentially use the same one (permuted through the map)
! but it is really not that hard to get profiles for all needed variables. 
                   do j=1,ndof   ! BC index
                        ic=Map2ICfromBC(j) !map to IC index
                        y(i,ic)=(xi*bcinterp(iupper,j+1)
     &                    +(one-xi)*bcinterp(iupper-1,j+1))
                   enddo
                 endif ! ICset

              endif ! either IBC active or icset active
           enddo
        endif  ! inlet.dat existed
c$$$$$$$$$$$$$$$$$$$$

!======================================================================
!Modifications for Duct. Need to encapsulate in a function call. 
        !specify an initial eddy viscosity ramp
        if(isetEVramp .gt. 0) then
          if(myrank .eq. 0) then
            write(*,*) "Setting eddy viscosity ramp with:" 
            write(*,*) "  - ramp X min = ", EVrampXmin
            write(*,*) "  - ramp X max = ", EVrampXmax
            write(*,*) "  - EV min = ", EVrampMin
            write(*,*) "  - EV max = ", EVrampMax
          endif
             
          x1 = EVrampXmin  !stuff in a shorter variable name to
          x2 = EVrampXmax  !make the formulas cleaner
          !Newton Divided differences to generate a polynomial of
          !the form P(x) = c0 + x*(c1 + x*(c2 + (x - (x2 - x1))*c3))
          !satisfying P(x1) = EVrampMin, P(x2) = EVrampMax,
          ! P'(x1) = 0, and P'(x2) = 0
          
          c0 = EVrampMin
          c1 = 0            !zero derivative
          c2 = (EVrampMax - EVrampMin)/(x2 - x1)
          c3 = 0            !zero derivative
          c3 = (c3 - c2)/(x2 - x1)
          c2 = (c2 - c1)/(x2 - x1)
          c3 = (c3 - c2)/(x2 - x1)
          
          do nn = 1,nshg
            if(y(nn,6) .eq. 0) cycle  !don't change wall boundary conditions, should be iTurbWall == 1
              
            if(point2x(nn,1) .gt. EVrampXmax) then !downstream of the ramp
              y(nn,6) = EVrampMax
            elseif(point2x(nn,1) .gt. EVrampXmin) then !and x(:,1) <= EVrampXmax
             
              !P(x) = c0 + x*(c1 + x*(c2 + (x - (x2 - x1))*c3)) 
              !     = c0 + x*(c1 + x*(c2 - (x2 - x1)*c3 + x*c3
              y(nn,6) = c0                 + (point2x(nn,1) - x1)*(
     &                  c1                 + (point2x(nn,1) - x1)*(
     &                 (c2 - (x2 - x1)*c3) + (point2x(nn,1) - x1)*c3))
            else
              y(nn,6) = EVrampMin
            endif
          enddo
        endif
!End modifications for Duct
!======================================================================
c.... Load a 2D field as an initial condition
        if(i2DIC == 1)then
           call set_wmles_ic(point2x, y)
        endif
c.... Load a 3D field as an initial condition
        if(i3DIC == 1)then
         if(myrank.eq.master) write(*,*) 
     &                       'Setting 3D IC from interpolated solution'
         maxdwl = maxval(d2wall)
         if (maxdwl.lt.0.04) then
           fname2 = "DNSIC/solLES."
           if (myrank.lt.9) then
             write (fname2, "(A13,I1)") trim(fname2), myrank+1
           else if (myrank.lt.99) then 
             write (fname2, "(A13,I2)") trim(fname2), myrank+1
           else if (myrank.lt.999) then
             write (fname2, "(A13,I3)") trim(fname2), myrank+1
           else if (myrank.lt.9999) then
             write (fname2, "(A13,I4)") trim(fname2), myrank+1
           else if (myrank.lt.99999) then
             write (fname2, "(A13,I5)") trim(fname2), myrank+1
           else if (myrank.lt.999999) then
             write (fname2, "(A13,I6)") trim(fname2), myrank+1
           endif           
           inquire(file=fname2,exist=exlog)
           if(exlog) then
              open (unit=456,file=fname2,status="old")
              allocate(IC3d(numnp,7))
              do i=1,numnp
                read(456,*) (IC3d(i,j),j=1,7)
              enddo
              close(456)
           else
              write(*,*) 'Did not find file ',fname2
           endif
           do i=1,numnp
              dx = abs(point2x(i,1)-IC3d(i,1))
              dy = abs(point2x(i,2)-IC3d(i,2))
              dz = abs(point2x(i,3)-IC3d(i,3))
              if (dx.lt.1.0e-5.and.dy.lt.1.0e-5.and.dz.lt.1.0e-5) then
                 y(i,1) = IC3d(i,5) ! u
                 y(i,2) = IC3d(i,6) ! v
                 y(i,3) = IC3d(i,7) ! w
                 !y(i,4) = IC3d(i,4) ! p
              else
                 write(*,*) 'PROBLEM: nodes do not match'
              endif
            enddo
            deallocate(IC3d)
         endif
        endif
c.... Load a 3D field as an initial condition
        if(iSTGIC == 1)then
         if(myrank.eq.master) write(*,*)
     &                       'Setting 3D IC from saved STG plane'
         maxdwl = maxval(d2wall)
         minx = minval(point2x(:,1))
         maxx = maxval(point2x(:,1))
         if ((minx.gt.-0.12)
     &       .and.(maxx.lt.0.1)) then
           inquire(file='STGplane.dat',exist=exlog)
           if(exlog) then
            open (unit=456,file='STGplane.dat',status="old")
            read(456,*) numSTG
            allocate(ICSTG(numSTG,5))
            do i=1,numSTG
              read(456,*) (ICSTG(i,j),j=1,5)
            enddo
            close(456)
            ! Note the ICSTG array form the file needs to be ordered by
            ! increasing dwal and increasing z
            do i=1,nshg
              dwl = d2wall(i)
              if (dwl.le.1.0e-3) then
              z = point2x(i,3)
              do j=1,numSTG
                if (ICSTG(j,1).ge.dwl) then ! found dwal
                  if (ICSTG(j,2).ge.z) then
                     y(i,1:3) = y(i,1:3)+ICSTG(j,3:5)
                     exit
                  endif
                endif
              enddo
              endif
            enddo
           deallocate(ICSTG)
           else
             write(*,*) 'Did not read STGplane.dat'
           endif
         endif
        endif
!======================================================================
c
c
c.... call the semi-discrete predictor multi-corrector iterative driver
c
        call itrdrv (y,              ac,             
     &               uold,           point2x,
     &               iBC,            BC,
     &               point2iper,     point2ilwork,   shp,
     &               shgl,           shpb,           shglb,
     &               point2ifath,    point2nsons,
     &               point2cdelsq ) 
c
c.... return
c
c
c.... stop CPU-timer
c
c
c.... close echo file
c

      if (numpe > 1) call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(myrank.eq.0)  then
          write(*,*) 'process - before closing iecho'
      endif

        close (iecho)

      if (numpe > 1) call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if(myrank.eq.0)  then
          write(*,*) 'process - after closing iecho'
      endif


c
c.... end of the program
c
        deallocate(point2iper)
        if(numpe.gt.1) then
          call Dctypes(point2ilwork(1))
        endif
        deallocate(point2ilwork)
        deallocate(point2x)
        deallocate(point2nsons)
        deallocate(point2ifath)
        if (allocated(velbar)) deallocate(velbar)
        if (allocated(stsBar)) deallocate(stsBar)
        if (allocated(stsBarKeq)) deallocate(stsBarKeq)
        if (allocated(velf)) deallocate(velf)
        if (allocated(velftG)) deallocate(velftG)
        if (allocated(locifath)) deallocate(locifath)
        if (allocated(ifathG)) deallocate(ifathG)
        if (allocated(rcounts)) deallocate(rcounts)
        if (allocated(displs)) deallocate(displs)
        if (allocated(velf)) deallocate(tmpStatsf)
        if (allocated(velftG)) deallocate(tmpStatsftG)
        if (allocated(velf)) deallocate(tmpKeqf)
        if (allocated(velftG)) deallocate(tmpKeqftG)
        deallocate(point2cdelsq)
        deallocate(uold)
        deallocate(wnrm)
        deallocate(otwn)
        call finalizeDtN
        call clearper
        call finalizeNABI

        return
        end


