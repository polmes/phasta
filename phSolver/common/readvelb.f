        subroutine readvelb (q, ifath, ifail)
c
c        use quadfilt
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
c
c        character*5 cname
c        character*255 fname1,fnamer,fmt1
c
        dimension q(nfath,nflow), ifath(numnp)   
        real*8    read1(nfath,nflow)
c        integer     irstin, lstep, itmp, iarray(10)


        q=0
        ifail=1
        read1=0
        call read_velbarS(myrank,lstep,ifail,read1)
        if(ifail.eq.0) then
          do i=1,nflow
           q(:,i)=read1(:,i)
          enddo
C..MODIFICATION FOR INCORPORATING TEMP. FIELD
CG          q(:,5)=q(:,2)
        endif

        return

        end

