      subroutine hello 
      use omp_lib
      include "common.h"
      integer idshared(4)
#ifdef HAVE_OMP
!$OMP  parallel do 
!$OMP& private (ith,id)
!$OMP& shared (idshared)
      do ith = 1, 4
          id = omp_get_thread_num ( )
!$OMP critical
          write (*,*) 'HELLO from rank, process, ith ', myrank,id,ith
!$OMP end critical
          idshared(ith)=id
      enddo
!     do ith = 1, 4
!         write (*,*) 'HELLO from rank, process, ith ', myrank,idshared(ith),ith
!     enddo
#endif
      write (*,*) 'Back to MPI ', myrank
      return
      end
