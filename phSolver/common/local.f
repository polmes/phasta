!
! new threaded blocks only node for local and local x
!
      subroutine localT (blk,global, rlocal, ientmp, n, code)
c
c----------------------------------------------------------------------
c
c This subroutine performs a vector gather/scatter operation.
c
c input:
c  global (nshg,n)             : global array
c  rlocal (bsz,blk%s,n)         : local array
c  ien    (blk%e,blk%s)      : nodal connectivity
c  n                            : number of d.o.f.'s to be copied
c  code                         : the transfer code
c                                  .eq. 'gather  ', from global to local
c                                  .eq. 'scatter ', add  local to global 
c                                  .eq. 'globaliz', from local to global
c
c
c Zdenek Johan, Winter 1992.
c----------------------------------------------------------------------
c
      use eblock
      include "common.h"
      type (LocalBlkData) blk


        dimension global(nshg,n),           rlocal(bsz,blk%s,n),
     &            ien(blk%e,blk%s),           ientmp(blk%e,blk%s)
c
        character*8 code
        
c
c.... cubic basis has negatives in ien
c
        if (blk%o > 2) then
           ien = abs(ientmp)
        else
           ien = ientmp
        endif
c
c.... ------------------------>  'localization  '  <--------------------
c
        if (code .eq. 'gather  ') then
c
c.... gather the data
c

          do j = 1, n
            do i = 1, blk%s
              rlocal(1:blk%e,i,j) = global(ien(1:blk%e,i),j)
            enddo
          enddo


c
c.... transfer count
c
          gbytes = gbytes + n*blk%s*blk%e
c
c.... return
c
          return
        endif
c
c.... ------------------------->  'assembling '  <----------------------
c
        if (code .eq. 'scatter ') then
c
c.... scatter the data (possible collisions)
c
          do j = 1, n
            do i = 1, blk%s
              do nel = 1,blk%e
                global(ien(nel,i),j) = global(ien(nel,i),j) 
     &                               + rlocal(nel,i,j)
              enddo
            enddo
          enddo

c
c.... transfer and flop counts
c
          sbytes = sbytes + n*blk%s*blk%e
          flops  = flops  + n*blk%s*blk%e
c
c.... return
c
          return
        endif
c
c.... ------------------------->  'globalizing '  <----------------------
c
        if (code .eq. 'globaliz') then
c
c.... scatter the data (possible collisions)
c
          do j = 1, n
            do i = 1, blk%s
              do nel = 1,blk%e
                global(ien(nel,i),j) = rlocal(nel,i,j)
              enddo
            enddo
          enddo
c
c.... return
c
          return
        endif
c
c.... --------------------------->  error  <---------------------------
c
        call error ('local   ', code, 0)
c
c.... end
c
        end
c
!
! new threaded blocks only node for local above
!
      subroutine local (blk,global, rlocal, ientmp, n, code)
c
c----------------------------------------------------------------------
c
c This subroutine performs a vector gather/scatter operation.
c
c input:
c  global (nshg,n)             : global array
c  rlocal (blk%e,blk%s,n)         : local array
c  ien    (blk%e,blk%s)      : nodal connectivity
c  n                            : number of d.o.f.'s to be copied
c  code                         : the transfer code
c                                  .eq. 'gather  ', from global to local
c                                  .eq. 'scatter ', add  local to global 
c                                  .eq. 'globaliz', from local to global
c
c
c Zdenek Johan, Winter 1992.
c----------------------------------------------------------------------
c
      use eblock
      include "common.h"
      type (LocalBlkData) blk


        dimension global(nshg,n),           rlocal(blk%e,blk%s,n),
     &            ien(blk%e,blk%s),           ientmp(blk%e,blk%s)
c
        character*8 code
        
c
c.... cubic basis has negatives in ien
c
        if (blk%o > 2) then
           ien = abs(ientmp)
        else
           ien = ientmp
        endif
c
c.... ------------------------>  'localization  '  <--------------------
c
        if (code .eq. 'gather  ') then
c
c.... gather the data
c

          do j = 1, n
            do i = 1, blk%s
              rlocal(1:blk%e,i,j) = global(ien(1:blk%e,i),j)
            enddo
          enddo


c
c.... transfer count
c
          gbytes = gbytes + n*blk%s*blk%e
c
c.... return
c
          return
        endif
c
c.... ------------------------->  'assembling '  <----------------------
c
        if (code .eq. 'scatter ') then
c
c.... scatter the data (possible collisions)
c
          do j = 1, n
            do i = 1, blk%s
              do nel = 1,blk%e
                global(ien(nel,i),j) = global(ien(nel,i),j) 
     &                               + rlocal(nel,i,j)
              enddo
            enddo
          enddo

c
c.... transfer and flop counts
c
          sbytes = sbytes + n*blk%s*blk%e
          flops  = flops  + n*blk%s*blk%e
c
c.... return
c
          return
        endif
c
c.... ------------------------->  'globalizing '  <----------------------
c
        if (code .eq. 'globaliz') then
c
c.... scatter the data (possible collisions)
c
          do j = 1, n
            do i = 1, blk%s
              do nel = 1,blk%e
                global(ien(nel,i),j) = rlocal(nel,i,j)
              enddo
            enddo
          enddo
c
c.... return
c
          return
        endif
c
c.... --------------------------->  error  <---------------------------
c
        call error ('local   ', code, 0)
c
c.... end
c
        end
c

c
        subroutine localx (blk,global, rlocal, ien, n, code)
c
c----------------------------------------------------------------------
c
c This subroutine performs a vector gather/scatter operation for the
c nodal coordinates array.
c
c input:
c  global (numnp,n)             : global array
c  rlocal (blk%e,blk%n,n)         : local array
c  ien    (blk%e,blk%s)      : nodal connectivity
c  n                            : number of d.o.f.'s to be copied
c  code                         : the transfer code
c                                  .eq. 'gather  ', from global to local
c                                  .eq. 'scatter ', add  local to global 
c
c
c Zdenek Johan, Winter 1992.
c----------------------------------------------------------------------
c
      use eblock
      include "common.h"
      type (LocalBlkData) blk


        dimension global(numnp,n),           rlocal(blk%e,blk%n,n),
     &            ien(blk%e,blk%s)
c
        character*8 code
c
c.... ------------------------>  'localization  '  <--------------------
c
        if (code .eq. 'gather  ') then
c
c.... gather the data
c
          do j = 1, n
            do i = 1, blk%n
              rlocal(1:blk%e,i,j) = global(ien(1:blk%e,i),j)
            enddo
          enddo


c
c.... transfer count
c
          gbytes = gbytes + n*blk%n*blk%e
c
c.... return
c
          return
        endif
c
c.... ------------------------->  'assembling '  <----------------------
c
        if (code .eq. 'scatter ') then
c
c.... scatter the data (possible collisions)
c

          do j = 1, n
            do i = 1, blk%n
              do nel = 1,blk%e
                global(ien(nel,i),j) = global(ien(nel,i),j) 
     &                               + rlocal(nel,i,j)
              enddo
            enddo
          enddo


c
c.... transfer and flop counts
c
          sbytes = sbytes + n*blk%s*blk%e
          flops  = flops  + n*blk%s*blk%e
c
c.... return
c
          return
        endif
c
c.... --------------------------->  error  <---------------------------
c
        call error ('local   ', code, 0)
c
c.... end
c
        end
c

        subroutine localSumT (blk,global, rlocal, ientmp, nHits, n)
c
c----------------------------------------------------------------------
c
c  sum the data from the local array to the global degrees of
c  freedom and keep track of the number of locals contributing
c  to each global dof. This may be used to find the average.
c
c----------------------------------------------------------------------
c
      use eblock
      include "common.h"
      type (LocalBlkData) blk

        dimension global(nshg,n),           rlocal(bsz,blk%s,n),
     &            ien(blk%e,blk%s),           ientmp(blk%e,blk%s),
     &            nHits(nshg)
c
c.... cubic basis has negatives in ien
c
        if (ipord > 2) then
           ien = abs(ientmp)
        else
           ien = ientmp
        endif
c
c.... ------------------------->  'assembling '  <----------------------
c
        do j = 1, n
           do i = 1, blk%s
              do nel = 1,blk%e
                 idg = ien(nel,i)
                 global(idg,j) = global(idg,j) + rlocal(nel,i,j)
              enddo
           enddo
        enddo
        do i = 1, blk%s
           do nel = 1,blk%e
              idg = ien(nel,i)
              nHits(idg) = nHits(idg) + 1
           enddo
        enddo
c
c.... end
c
        end
        subroutine localSum (blk,global, rlocal, ientmp, nHits, n)
c
c----------------------------------------------------------------------
c
c  sum the data from the local array to the global degrees of
c  freedom and keep track of the number of locals contributing
c  to each global dof. This may be used to find the average.
c
c----------------------------------------------------------------------
c
      use eblock
      include "common.h"
      type (LocalBlkData) blk

        dimension global(nshg,n),           rlocal(blk%e,blk%s,n),
     &            ien(blk%e,blk%s),           ientmp(blk%e,blk%s),
     &            nHits(nshg)
c
c.... cubic basis has negatives in ien
c
        if (ipord > 2) then
           ien = abs(ientmp)
        else
           ien = ientmp
        endif
c
c.... ------------------------->  'assembling '  <----------------------
c
        do j = 1, n
           do i = 1, blk%s
              do nel = 1,blk%e
                 idg = ien(nel,i)
                 global(idg,j) = global(idg,j) + rlocal(nel,i,j)
              enddo
           enddo
        enddo
        do i = 1, blk%s
           do nel = 1,blk%e
              idg = ien(nel,i)
              nHits(idg) = nHits(idg) + 1
           enddo
        enddo
c
c.... end
c
        end
