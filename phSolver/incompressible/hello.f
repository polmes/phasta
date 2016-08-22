      subroutine hello 
      use omp_lib
      include "common.h"
!$OMP  parallel do 
!$OMP& private (ith,id)
      do ith = 1, 4
          id = omp_get_thread_num ( )
          write (*,*) 'HELLO from rank, process, ith ', myrank,id,ith
      enddo
      write (*,*) 'Back to MPI ', myrank
      return
      end
