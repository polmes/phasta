c-----------------------------------------------------------------------
        subroutine BCprofileScale(vbc_prof,BC,yold)
        use pvsQbi
        include "common.h"
        real*8 vbc_prof(nshg,3)
        real*8 BC(nshg,ndofBC),yold(nshg,ndof)
        real*8 offphase
        integer factor

        r_amp =rampmdot(1,1)
        r_freq=rampmdot(2,1)
!  Usual sinusoidal in time syn jet is given in the next line.  It assumes that you are NOT changing the time step from previous runs since it computes the 
!  current time in the sin function to be (lstep+1)* Dt where lstep is a running total of all of the runs up to now.  This will be "ok" under the following
!  assumption:  lstep*dt_current = m * T where m is an integer and T is the period of the jet.  That is because this will evaluate to 2*m*Pi and the sin
!  that is zero

        if(abs(rampmdot(1,2)) .lt. 1.0e-5) then
           r_time_factor = r_amp*sin(two*pi*r_freq*(lstep+1)*Delt(1)) ! BC set for next time step - all phases are computed
        else
           r_time_factor = rampmdot(1,2)*r_amp   ! read parameter scales Vmax
        endif

        

        icount = 0
        do kk=1,nshg
          if(ndsurf(kk).eq.18) then ! this means diaphragm for the Cube Test case

            offphase = 1.d0
! 
! the following chunk of code will zero out the amplitude  for (iduty-1) cycles of the jet before 
! using the time varying amplitude computed above
!
            iduty=idnint(rampmdot(2,2))
            if(iduty.gt.1) then
                nperiods=r_freq*(lstep+1)*Delt(1)  ! compute period number of the current step. NOTE lstep step number across all runs
                imod=mod(nperiods,iduty)           ! will be the remeinder of nperiods/iduty
                if(imod.gt.0) offphase=0           ! set to zero except for the period with no remainder
            endif

            BC(kk,3)=r_time_factor*vbc_prof(kk,1)*offphase
            BC(kk,4)=r_time_factor*vbc_prof(kk,2)*offphase
            BC(kk,5)=r_time_factor*vbc_prof(kk,3)*offphase

            icount = icount + 1

          endif
        enddo

!        if(istep.eq.0 .and. icount.ne.0)
!     &     write(*,*) 'BCprofile count',myrank,icount
    
        return
        end 

c--------------------------------------------------------------
        subroutine BCprofileInit(vbc_prof,x)
      
        use pvsQbi
        include "common.h"
        real*8 vbc_prof(nshg,3), x(numnp,nsd)
        real*8 rcenter(3),rnorml(3)
        do kk=1,nshg
c.............Factors below are negative for desired blowing direction
           if(ndsurf(kk).eq.18) then
             x1=-0.025d0
             x2=0.025d0
             z1=-0.01d0
             z2=0.01d0
             vbc_prof(kk,1)=0.d0
             vbc_prof(kk,2)=-1.6*1E7*(x(kk,1)-x1)*(x(kk,1)-x2)*
     &                              (x(kk,3)-z1)*(x(kk,3)-z2)
             vbc_prof(kk,3)=0.d0 
           else
              vbc_prof(kk,:)=zero
           endif

        enddo
        close(789)
987     format(i6,3(2x,e14.7))

        return
        end
