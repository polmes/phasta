c48---------------------------------------------------------------------
c
c ftools.f : Bundle of Fortran routines
c
c Each routine is to be called by drvftools.f
c
c Various operations run based on les**.c
c
c---------------------------------------------------------------------
c--------------------------------
c fMtxVdimVecMult 
c Farzin's implementation
c row and column index exchanged
c--------------------------------
c
        subroutine fMtxVdimVecMult( a, b, c, na, nb, nc, m, n,rblasvvm,iblasvvm,ieqswork )
c
c.... Data declaration
c
        implicit none
        include "mymkl_vml.fi"
        integer na,     nb,     nc,     m,      n
        real*8  a(n,na),        b(n,nb),        c(n,nc)
c
        integer i,      j
        integer iwork, iblasvvm, ieqswork
        real*8 rdelta,TMRC,rblasvvm
        rdelta=TMRC()
        iwork=mod(ieqswork,10)
        if(iwork.eq.2) then
         do j=1,m
           call vdmul(n,a(1,j),b(1,j),c(1,j))
         enddo
        else if(iwork.eq.8) then
            do i = 1, n 
                do j = 1, m 
                    c(i,j) = a(i,j) * b(i,j)
                enddo
            enddo
        else
c
c.... Do the work
c
C WIP: change to F90
C
        if ( m .eq. 1 ) then
c
            do i = 1, n
                c(i,1) = a(i,1) * b(i,1)
            enddo
c
        else if ( m .eq. 2 ) then
c
            do i = 1, n
                c(i,1) = a(i,1) * b(i,1)
                c(i,2) = a(i,2) * b(i,2)
            enddo
c
        else if ( m .eq. 3 ) then
c
            do i = 1, n
                c(i,1) = a(i,1) * b(i,1)
                c(i,2) = a(i,2) * b(i,2)
                c(i,3) = a(i,3) * b(i,3)
            enddo
c
        else if ( m .eq. 4 ) then
c
            do i = 1, n
                c(i,1) = a(i,1) * b(i,1)
                c(i,2) = a(i,2) * b(i,2)
                c(i,3) = a(i,3) * b(i,3)
                c(i,4) = a(i,4) * b(i,4)
            enddo
c
        else
c
            do i = 1, n 
                do j = 1, m 
                    c(i,j) = a(i,j) * b(i,j)
                enddo
            enddo
c
        endif
      endif
      rblasvvm=rblasvvm+TMRC()-rdelta
      iblasvvm=iblasvvm+m
c
      return
      end
c
c---------- 
c flesZero
c----------
c
	subroutine flesZero ( a, m, n )
c
        implicit none
        integer  m, n, i, j
        real*8   a(n,m)
c
        do i = 1, n
          do j = 1, m
            a(i,j) = 0.0e-0
          enddo
        enddo
c
        return 
        end
c
c--------
c flesCp
c--------
c
	subroutine flesCp ( a, b, m, n )
c
        implicit none
        integer  m, n, i, j
        real*8   a(n,m), b(n,m)
c
        do i = 1, n
          do j = 1, m
            b(i,j) = a(i,j)
          enddo
        enddo
c
        return
        end
c
c-----------
c flesScale
c-----------
c
 	subroutine flesScale ( a, s, m, n )
c
        implicit none
        integer  m, n, i, j   
        real*8   a(n,m), s
c
        do i = 1, n
          do j = 1, m
            a(i,j) = a(i,j) * s
          enddo
        enddo
c
        return
        end
c
c-------------
c flesScaleCp
c-------------
c
	subroutine flesScaleCp ( a, b, s, m, n )
c
        implicit none
        integer  m, n, i, j
        real*8   a(n,m), b(n,m), s
c
        do i = 1, n
          do j = 1, m
            b(i,j) = a(i,j) * s
          enddo
        enddo
c
        return 
        end
c
c---------
c flesAdd
c---------
c
	subroutine flesAdd ( a, b, m, n )
c
        implicit none
        integer  m, n, i, j
        real*8   a(n,m), b(n,m)
c
        do i = 1, n
          do j = 1, m
            b(i,j) = b(i,j) + a(i,j)
          enddo
        enddo
c
        return
        end
c
c---------
c flesSub 
c---------
c
	subroutine flesSub ( a, b, m, n )
c
        implicit none
        integer  m, n, i, j
        real*8   a(n,m), b(n,m)
c
        do i = 1, n
          do j = 1, m
            b(i,j) = b(i,j) - a(i,j)
          enddo
        enddo
c
        return
        end
c
c----------
c flesDot1
c----------
c
 	real*8 function flesDot1 ( a, m, n )
c
        implicit none
        integer  m, n, i, j
        real*8   a(n,m)
c
	flesDot1 = 0
        do i = 1, n
          do j = 1, m
            flesDot1 = flesDot1 + a(i,j) * a(i,j)
          enddo
        enddo
c
        return 
        end
c
c----------
c flesDot2
c----------
c
	real*8 function flesDot2 ( a, b, m, n )
c
        implicit none
        integer  m, n, i, j
        real*8   a(n,m), b(n,m)
c
	flesDot2 = 0
        do i = 1, n
          do j = 1, m
            flesDot2 = flesDot2 + a(i,j) * b(i,j)
          enddo
        enddo
c
        return
        end
c
c-----------
c flesDaxpy
c-----------
c
	subroutine flesDaxpy ( x, y, a, m, n )
c
        implicit none
        include "mymkl.fi"
        integer  m, n, i, j
        real*8   x(n,m), y(n,m)
        real*8   a
c
        if(1.eq.0) then
         do j=1,m
           call daxpy(n,a,x(1,j),1,y(1,j),1)
         enddo
        else
         do i = 1, n
          do j= 1, m
            y(i,j) = y(i,j) + a * x(i,j)
          enddo
         enddo
        endif
c
        return
        end
c
c-----------
c flesDxpay
c-----------
c
	subroutine flesDxpay ( x, y, a, m, n )
c
        implicit none
        integer  m, n, i, j
        real*8   x(n,m), y(n,m)
        real*8   a
c
        do i = 1, n
          do j = 1, m
            y(i,j) = a * y(i,j) + x(i,j)
          enddo
        enddo
c
        return 
        end
c
c---------
c flesInv
c---------
c
	subroutine flesInv ( x, m, n )
c
        implicit none
        integer  m, n, i, j
        real*8   x(n,m)
c
        do i = 1, n
          do j = 1, m
            if ( x(i,j) .ne. 0 ) x(i,j) = 1. / x(i,j)
          enddo
        enddo

c
        return
        end
c
c--------------------------
c fMtxBlkDot2
c Farzin's implementation
c row and column exchanged
c--------------------------
c
        subroutine fMtxBlkDot2( x, y, c, m, n, rblasphasta,rblasmkl,
     &             iblasphasta,iblasmkl,ieqswork )

        implicit none
        include "mymkl.fi"
!        include "blas.f90"
!        include "kmkl.fi"

c
c.... Data declaration
c
        implicit none
        integer m,      n
        real*8  x(n,m), y(n),   c(m), d(m)

        real*8 alpha,beta
        real*8 rdelta, TMRC, rblasphasta,rblasmkl
c
        real*8  tmp1,   tmp2,   tmp3,   tmp4
        real*8  tmp5,   tmp6,   tmp7,   tmp8
        integer i,      j,      m1,lda,incx,incy
        integer iwork, ieqswork,     icut
        integer iblasphasta,      iblasmkl
!DIR$ ASSUME_ALIGNED x: 64, y:64, c:64

      iwork=mod(ieqswork,10)

!       0 chunk8, 1 chunk4, 2  mkl ddot, 3  mkl dgemmv 
!       4-6  3 for m< icut below, else chunk8, 7 open, 8 double loop 
        
        if((iwork.ge.3).and.(iwork.lt.7)) then ! let dgemv do some/all of the work
          icut=24*2**(iwork-4)
          if(m.gt.icut) then
            alpha=1.0
            beta=0.0
            incx=1
            incy=1
            lda=n
        
!   note matrix is expected to be a(m,n) but ours is x(n,m)
!   I assume we handle this by n<->m in arguments 2, 3
!   since m is the number of rows of A NOT A^T 
!  further complicating things A is x a 
!  their x is our y
!  their y is our output 
!  them    y=alpha*A^T x +beta*y
!  us     d=alpha*x^T y +beta*d
            rdelta=TMRC()
            call dgemv('T',n,m,alpha,x,lda,y,incx,beta,c,incy)
            rblasmkl=rblasmkl+TMRC()-rdelta
            iblasmkl=iblasmkl+m
            return
          else ! jump to chunk 8 to do the rest of the work or 
              ! change iwork  to who you want to do it
            iwork=0
            goto 1
          endif
        endif ! end of dgemv
1       continue 
        if(iwork.eq.2) then  
cdir$ ivdep
          rdelta=TMRC()
          do i=1,m
            c(i)=ddot(n,x(1,i),1,y,1)
          enddo
          rblasmkl=rblasmkl+TMRC()-rdelta
          iblasmkl=iblasmkl+m
          return
        endif
        if(iwork.eq.8) then
          rdelta=TMRC()
          c=0
cdir$ ivdep
          do j = 1, m
            do i = 1, n
              c(j) = c(j) + x(i,j) * y(i)
            enddo
          enddo
          rblasphasta=rblasphasta+TMRC()-rdelta
          iblasphasta=iblasphasta+m
          return
        endif
        if(iwork.eq.0) then
          rdelta=TMRC()

c
c.... Determine the left overs
c
        m1 = mod(m,8) + 1

c
c.... Do the small pieces
c
        goto ( 8000, 1000, 2000, 3000, 4000, 5000, 6000, 7000 ) m1
c
1000    continue
        tmp1 = 0
        do i = 1, n
            tmp1 = tmp1 + x(i,1) * y(i)
        enddo
        c(1) = tmp1
        goto 8000
c
2000    continue
        tmp1 = 0
        tmp2 = 0
        do i = 1, n
            tmp1 = tmp1 + x(i,1) * y(i)
            tmp2 = tmp2 + x(i,2) * y(i)
        enddo
        c(1) = tmp1
        c(2) = tmp2
        goto 8000
c
3000    continue
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        do i = 1, n
            tmp1 = tmp1 + x(i,1) * y(i)
            tmp2 = tmp2 + x(i,2) * y(i)
            tmp3 = tmp3 + x(i,3) * y(i)
        enddo
        c(1) = tmp1
        c(2) = tmp2
        c(3) = tmp3
        goto 8000
c
4000    continue
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        tmp4 = 0
        do i = 1, n
            tmp1 = tmp1 + x(i,1) * y(i)
            tmp2 = tmp2 + x(i,2) * y(i)
            tmp3 = tmp3 + x(i,3) * y(i)
            tmp4 = tmp4 + x(i,4) * y(i)
        enddo
        c(1) = tmp1
        c(2) = tmp2
        c(3) = tmp3
        c(4) = tmp4
        goto 8000
c
5000    continue
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        tmp4 = 0
        tmp5 = 0
        do i = 1, n
            tmp1 = tmp1 + x(i,1) * y(i)
            tmp2 = tmp2 + x(i,2) * y(i)
            tmp3 = tmp3 + x(i,3) * y(i)
            tmp4 = tmp4 + x(i,4) * y(i)
            tmp5 = tmp5 + x(i,5) * y(i)
        enddo
        c(1) = tmp1
        c(2) = tmp2
        c(3) = tmp3
        c(4) = tmp4
        c(5) = tmp5
        goto 8000
c
6000    continue
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        tmp4 = 0
        tmp5 = 0
        tmp6 = 0
        do i = 1, n
            tmp1 = tmp1 + x(i,1) * y(i)
            tmp2 = tmp2 + x(i,2) * y(i)
            tmp3 = tmp3 + x(i,3) * y(i)
            tmp4 = tmp4 + x(i,4) * y(i)
            tmp5 = tmp5 + x(i,5) * y(i)
            tmp6 = tmp6 + x(i,6) * y(i)
        enddo
        c(1) = tmp1
        c(2) = tmp2
        c(3) = tmp3
        c(4) = tmp4
        c(5) = tmp5
        c(6) = tmp6
        goto 8000
c
7000    continue
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        tmp4 = 0
        tmp5 = 0
        tmp6 = 0
        tmp7 = 0
        do i = 1, n
            tmp1 = tmp1 + x(i,1) * y(i)
            tmp2 = tmp2 + x(i,2) * y(i)
            tmp3 = tmp3 + x(i,3) * y(i)
            tmp4 = tmp4 + x(i,4) * y(i)
            tmp5 = tmp5 + x(i,5) * y(i)
            tmp6 = tmp6 + x(i,6) * y(i)
            tmp7 = tmp7 + x(i,7) * y(i)
        enddo
        c(1) = tmp1
        c(2) = tmp2
        c(3) = tmp3
        c(4) = tmp4
        c(5) = tmp5
        c(6) = tmp6
        c(7) = tmp7
        goto 8000
c
c.... Do the remaining part
c
8000    continue
c
        do j = m1, m, 8
            tmp1 = 0
            tmp2 = 0
            tmp3 = 0
            tmp4 = 0
            tmp5 = 0
            tmp6 = 0
            tmp7 = 0
            tmp8 = 0
            do i = 1, n
                tmp1 = tmp1 + x(i,j) * y(i)
                tmp2 = tmp2 + x(i,j+1) * y(i)
                tmp3 = tmp3 + x(i,j+2) * y(i)
                tmp4 = tmp4 + x(i,j+3) * y(i)
                tmp5 = tmp5 + x(i,j+4) * y(i)
                tmp6 = tmp6 + x(i,j+5) * y(i)
                tmp7 = tmp7 + x(i,j+6) * y(i)
                tmp8 = tmp8 + x(i,j+7) * y(i)
            enddo
            c(j+0) = tmp1
            c(j+1) = tmp2
            c(j+2) = tmp3
            c(j+3) = tmp4
            c(j+4) = tmp5
            c(j+5) = tmp6
            c(j+6) = tmp7
            c(j+7) = tmp8
        enddo
        rblasphasta=rblasphasta+TMRC()-rdelta
        iblasphasta=iblasphasta+m
        return
        endif !iwork=0
        if(iwork.eq.1) then ! try a max 4 version of original
          rdelta=TMRC()
c
c.... Determine the left overs
c
        m1 = mod(m,4) + 1

c
c.... Do the small pieces
c
        goto ( 4001, 1001, 2001, 3001 ) m1
c
1001    continue
        tmp1 = 0
        do i = 1, n
            tmp1 = tmp1 + x(i,1) * y(i)
        enddo
        c(1) = tmp1
        goto 4001
c
2001    continue
        tmp1 = 0
        tmp2 = 0
        do i = 1, n
            tmp1 = tmp1 + x(i,1) * y(i)
            tmp2 = tmp2 + x(i,2) * y(i)
        enddo
        c(1) = tmp1
        c(2) = tmp2
        goto 4001
c
3001    continue
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        do i = 1, n
            tmp1 = tmp1 + x(i,1) * y(i)
            tmp2 = tmp2 + x(i,2) * y(i)
            tmp3 = tmp3 + x(i,3) * y(i)
        enddo
        c(1) = tmp1
        c(2) = tmp2
        c(3) = tmp3
        goto 4001
c
c.... Do the remaining part
c
4001   continue
c
        do j = m1, m, 4
            tmp1 = 0
            tmp2 = 0
            tmp3 = 0
            tmp4 = 0
            do i = 1, n
                tmp1 = tmp1 + x(i,j+0) * y(i)
                tmp2 = tmp2 + x(i,j+1) * y(i)
                tmp3 = tmp3 + x(i,j+2) * y(i)
                tmp4 = tmp4 + x(i,j+3) * y(i)
            enddo
            c(j+0) = tmp1
            c(j+1) = tmp2
            c(j+2) = tmp3
            c(j+3) = tmp4
        enddo
        rblasphasta=rblasphasta+TMRC()-rdelta
        iblasphasta=iblasphasta+m
        endif
c
        return
        end
c
c--------------------------
c fMtxBlkDaxpy
c Farzin's implementation 
c row and column exchanged
c--------------------------
c
        subroutine fMtxBlkDaxpy( x, y, c, m, n )
c
c.... Data declaration
c
        implicit none
        integer m,      n
        real*8  x(n,m), y(n),   c(m)
c
        real*8  tmp1,   tmp2,   tmp3,   tmp4
        real*8  tmp5,   tmp6,   tmp7,   tmp8
        integer i,      j,      m1
c
c.... Determine the left overs
c
        m1 = mod(m,8) + 1
c
c.... Do the small pieces
c
        goto ( 8000, 1000, 2000, 3000, 4000, 5000, 6000, 7000 ) m1
c
1000    continue
        tmp1 = c(1)
        do i = 1, n
            y(i) = y(i)
     1           + tmp1 * x(i,1)
        enddo
        goto 8000
c
2000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        do i = 1, n
            y(i) = y(i)
     1           + tmp1 * x(i,1) + tmp2 * x(i,2)
        enddo
        goto 8000
c
3000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)

        do i = 1, n
            y(i) = y(i)
     1           + tmp1 * x(i,1) + tmp2 * x(i,2)
     2           + tmp3 * x(i,3)
        enddo
        goto 8000
c
4000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        do i = 1, n
            y(i) = y(i)
     1           + tmp1 * x(i,1) + tmp2 * x(i,2)
     2           + tmp3 * x(i,3) + tmp4 * x(i,4)
        enddo
        goto 8000
c
5000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        tmp5 = c(5)
        do i = 1, n
            y(i) = y(i)
     1           + tmp1 * x(i,1) + tmp2 * x(i,2)
     2           + tmp3 * x(i,3) + tmp4 * x(i,4)
     3           + tmp5 * x(i,5)
        enddo
        goto 8000
c
6000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        tmp5 = c(5)
        tmp6 = c(6)
        do i = 1, n
            y(i) = y(i)
     1           + tmp1 * x(i,1) + tmp2 * x(i,2)
     2           + tmp3 * x(i,3) + tmp4 * x(i,4)
     3           + tmp5 * x(i,5) + tmp6 * x(i,6)
        enddo
        goto 8000
c
7000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        tmp5 = c(5)
        tmp6 = c(6)
        tmp7 = c(7)
        do i = 1, n
            y(i) = y(i)
     1           + tmp1 * x(i,1) + tmp2 * x(i,2)
     2           + tmp3 * x(i,3) + tmp4 * x(i,4)
     3           + tmp5 * x(i,5) + tmp6 * x(i,6)
     4           + tmp7 * x(i,7)
        enddo
        goto 8000
c
c.... Do the remaining part
c
8000    continue
c
        do j = m1, m, 8
            tmp1 = c(j+0)
            tmp2 = c(j+1)
            tmp3 = c(j+2)
            tmp4 = c(j+3)
            tmp5 = c(j+4)
            tmp6 = c(j+5)
            tmp7 = c(j+6)
            tmp8 = c(j+7)
            do i = 1, n
                y(i) = y(i)
     1               + tmp1 * x(i,j+0) + tmp2 * x(i,j+1)
     2               + tmp3 * x(i,j+2) + tmp4 * x(i,j+3)
     3               + tmp5 * x(i,j+4) + tmp6 * x(i,j+5)
     4               + tmp7 * x(i,j+6) + tmp8 * x(i,j+7)
            enddo
        enddo
c
        return 
        end
c
c--------------------------
c fMtxBlkDyeax
c Farzin's implementation
c row and column exchanged
c--------------------------
c
        subroutine fMtxBlkDyeax( x, y, c, m, n )
c
c.... Data declaration
c
        implicit none
!        include "mymkl.fi"
        integer m,      n
        real*8  x(n,m), y(n),   c(m)
c
        real*8  tmp1,   tmp2,   tmp3,   tmp4
        real*8  tmp5,   tmp6,   tmp7,   tmp8
        integer i,      j,      m1
        if(1.eq.0) then
          tmp1=c(1)
          do i = 1, n
              y(i) = 
     1           + tmp1 * x(i,1)
          enddo
          do j=2,m
           tmp1=c(j)
           call daxpy(n,tmp1,x(1,j),1,y,1)
          enddo
        else
       
c
c.... Determine the left overs
c
        m1 = mod(m,8) + 1
c
c.... Do the small pieces
c
        goto ( 8000, 1000, 2000, 3000, 4000, 5000, 6000, 7000 ) m1
c
1000    continue
        tmp1 = c(1)
        do i = 1, n
            y(i) =
     1           + tmp1 * x(i,1)
        enddo
        goto 8001
c
2000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        do i = 1, n
            y(i) =
     1           + tmp1 * x(i,1) + tmp2 * x(i,2)
        enddo
        goto 8001
c
3000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        do i = 1, n
            y(i) =
     1           + tmp1 * x(i,1) + tmp2 * x(i,2)
     2           + tmp3 * x(i,3)
        enddo
        goto 8001
c
4000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        do i = 1, n
            y(i) =
     1           + tmp1 * x(i,1) + tmp2 * x(i,2)
     2           + tmp3 * x(i,3) + tmp4 * x(i,4)
        enddo
        goto 8001
c
5000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        tmp5 = c(5)
        do i = 1, n
            y(i) =
     1           + tmp1 * x(i,1) + tmp2 * x(i,2)
     2           + tmp3 * x(i,3) + tmp4 * x(i,4)
     3           + tmp5 * x(i,5)
        enddo
        goto 8001
c
6000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        tmp5 = c(5)
        tmp6 = c(6)
        do i = 1, n
            y(i) =
     1           + tmp1 * x(i,1) + tmp2 * x(i,2)
     2           + tmp3 * x(i,3) + tmp4 * x(i,4)
     3           + tmp5 * x(i,5) + tmp6 * x(i,6)
       enddo
        goto 8001
c
7000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        tmp5 = c(5)
        tmp6 = c(6)
        tmp7 = c(7)
        do i = 1, n
            y(i) =
     1           + tmp1 * x(i,1) + tmp2 * x(i,2)
     2           + tmp3 * x(i,3) + tmp4 * x(i,4)
     3           + tmp5 * x(i,5) + tmp6 * x(i,6)
     4           + tmp7 * x(i,7)
        enddo
        goto 8001
c
8000    continue
        do i = 1, n
            y(i) = 0
        enddo
        goto 8001
c
c.... Do the remaining part
c
8001    continue
c
        do j = m1, m, 8
            tmp1 = c(j+0)
            tmp2 = c(j+1)
            tmp3 = c(j+2)
            tmp4 = c(j+3)
            tmp5 = c(j+4)
            tmp6 = c(j+5)
            tmp7 = c(j+6)
            tmp8 = c(j+7)
            do i = 1, n
                y(i) = y(i)
     1               + tmp1 * x(i,j+0) + tmp2 * x(i,j+1)
     2               + tmp3 * x(i,j+2) + tmp4 * x(i,j+3)
     3               + tmp5 * x(i,j+4) + tmp6 * x(i,j+5)
     4               + tmp7 * x(i,j+6) + tmp8 * x(i,j+7)
            enddo
        enddo
        endif
c
        return 
        end
c
c--------------------------
c fMtxBlkDmaxpy
c Farzin's implementation
c row and column exchanged 
c--------------------------
c
       subroutine fMtxBlkDmaxpy( x, y, c, m, n, rblasmaxpy,iblasmaxpy,ieqswork )
c
c.... Data declaration
c
        implicit none
        include "mymkl.fi"
!        include "blas.f90"
        integer m,      n
        real*8  x(n,m), y(n),   c(m)
c
        real*8  tmp1,   tmp2,   tmp3,   tmp4
        real*8  tmp5,   tmp6,   tmp7,   tmp8
        integer i,      j,      m1
        integer iwork,iblasmaxpy, ieqswork
        real*8 a,rdelta,TMRC,rblasmaxpy
        rdelta=TMRC()
!        include "mymkl.fi"
        iwork=mod(ieqswork,10)
        if(iwork.eq.2) then
          do j=1,m
            a=-1.0d0*c(j)
            call daxpy(n,a,x(1,j),1,y,1)
          enddo
        else if(iwork.eq.8) then
          do j=1,m
           tmp1=c(j)
           do i = 1, n
             y(i) = y(i) - tmp1 * x(i,j)
           enddo
          enddo
        else
c
c.... Determine the left overs
c
        m1 = mod(m,8) + 1
c
c.... Do the small pieces
c
        goto ( 8000, 1000, 2000, 3000, 4000, 5000, 6000, 7000 ) m1
c
1000    continue
        tmp1 = c(1)
        do i = 1, n
            y(i) = y(i)
     1           - tmp1 * x(i,1)
        enddo
        goto 8000
c
2000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        do i = 1, n
            y(i) = y(i)
     1           - tmp1 * x(i,1) - tmp2 * x(i,2)
        enddo
        goto 8000
c
3000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        do i = 1, n
            y(i) = y(i)
     1           - tmp1 * x(i,1) - tmp2 * x(i,2)
     2           - tmp3 * x(i,3)
        enddo
        goto 8000
c
4000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        do i = 1, n
            y(i) = y(i)
     1           - tmp1 * x(i,1) - tmp2 * x(i,2)
     2           - tmp3 * x(i,3) - tmp4 * x(i,4)
        enddo
        goto 8000
c
5000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        tmp5 = c(5)
        do i = 1, n
            y(i) = y(i)
     1           - tmp1 * x(i,1) - tmp2 * x(i,2)
     2           - tmp3 * x(i,3) - tmp4 * x(i,4)
     3           - tmp5 * x(i,5)
        enddo
        goto 8000
c
6000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        tmp5 = c(5)
        tmp6 = c(6)
        do i = 1, n
            y(i) = y(i)
     1           - tmp1 * x(i,1) - tmp2 * x(i,2)
     2           - tmp3 * x(i,3) - tmp4 * x(i,4)
     3           - tmp5 * x(i,5) - tmp6 * x(i,6)
        enddo
        goto 8000

7000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        tmp5 = c(5)
        tmp6 = c(6)
        tmp7 = c(7)
        do i = 1, n
            y(i) = y(i)
     1           - tmp1 * x(i,1) - tmp2 * x(i,2)
     2           - tmp3 * x(i,3) - tmp4 * x(i,4)
     3           - tmp5 * x(i,5) - tmp6 * x(i,6)
     4           - tmp7 * x(i,7)
        enddo
        goto 8000
c
c.... Do the remaining part
c
8000    continue
c
        do j = m1, m, 8
            tmp1 = c(j+0)
            tmp2 = c(j+1)
            tmp3 = c(j+2)
            tmp4 = c(j+3)
            tmp5 = c(j+4)
            tmp6 = c(j+5)
            tmp7 = c(j+6)
            tmp8 = c(j+7)
            do i = 1, n
                y(i) = y(i)
     1               - tmp1 * x(i,j+0) - tmp2 * x(i,j+1)
     2               - tmp3 * x(i,j+2) - tmp4 * x(i,j+3)
     3               - tmp5 * x(i,j+4) - tmp6 * x(i,j+5)
     4               - tmp7 * x(i,j+6) - tmp8 * x(i,j+7)
            enddo
        enddo
c
        endif
        rblasmaxpy=rblasmaxpy+TMRC()-rdelta
        iblasmaxpy=iblasmaxpy+m
        return
        end
c
c--------------------------
c fMtxVdimVecCp
c Farzin's implementation
c row and column exchanged 
c--------------------------
c
        subroutine fMtxVdimVecCp( a, b, na, nb, m, n )
c
c.... Data declaration
c
        implicit none
        include "mymkl.fi"
        integer na,     nb,     m,      n
        real*8  a(n,na),        b(n,nb)
c
        integer i,      j
!        include "mymkl.fi"
        if(1.eq.0) then  ! this is slower
         do j=1,m
           call dcopy(n,a(1,j),1,b(1,j),1)
         enddo
        else
c
c.... Do the work
c
        if ( m .eq. 1 ) then

            do i = 1, n
                b(i,1) = a(i,1)
            enddo

        else if ( m .eq. 2 ) then

            do i = 1, n
                b(i,1) = a(i,1)
                b(i,2) = a(i,2)
            enddo

        else if ( m .eq. 3 ) then

            do i = 1, n
                b(i,1) = a(i,1)
                b(i,2) = a(i,2)
                b(i,3) = a(i,3)
            enddo

        else if ( m .eq. 4 ) then

            do i = 1, n
                b(i,1) = a(i,1)
                b(i,2) = a(i,2)
                b(i,3) = a(i,3)
                b(i,4) = a(i,4)
            enddo

        else

            do i = 1, n
                do j = 1, m 
                    b(i,j) = a(i,j)
                enddo
            enddo

        endif
        endif
c
        return
        end
c
c--------------------------
c fMtxVdimVecDot2
c Farzin's implementation
c row and column exchanged
c--------------------------
c
        subroutine fMtxVdimVecDot2( a, b, c, na, nb, m, n )
c
c.... Data declaration
c
        implicit none
        integer na,     nb,     m,      n
        real*8  a(n,na),        b(n,nb),        c(m)
c
        integer i,      j
c
c.... Do the work
c
        if ( m .eq. 1 ) then

            c(1) = 0
            do i = 1, n
                c(1) = c(1) + a(i,1) * b(i,1)
            enddo

        else if ( m .eq. 2 ) then

            c(1) = 0
            c(2) = 0
            do i = 1, n
                c(1) = c(1) + a(i,1) * b(i,1)
                c(2) = c(2) + a(i,2) * b(i,2)
            enddo

        else if ( m .eq. 3 ) then

            c(1) = 0
            c(2) = 0
            c(3) = 0
            do i = 1, n
                c(1) = c(1) + a(i,1) * b(i,1)
                c(2) = c(2) + a(i,2) * b(i,2)
                c(3) = c(3) + a(i,3) * b(i,3)
            enddo

        else if ( m .eq. 4 ) then

            c(1) = 0
            c(2) = 0
            c(3) = 0
            c(4) = 0
            do i = 1, n
                c(1) = c(1) + a(i,1) * b(i,1)
                c(2) = c(2) + a(i,2) * b(i,2)
                c(3) = c(3) + a(i,3) * b(i,3)
                c(4) = c(4) + a(i,4) * b(i,4)
            enddo

        else

            do j = 1, m 
                c(j) = 0
                do i = 1, n 
                    c(j) = c(j) + a(i,j) * b(i,j)
                enddo
            enddo

        endif
c
        return
        end
c
c--------------------------
c fMtxVdimVecDaxpy
c Farzin's implementation
c row and column exchanged
c--------------------------
c
        subroutine fMtxVdimVecDaxpy( a, b, c, na, nb, m, n )
c
c.... Data declaration
c
        implicit none
        integer na,     nb,     m,      n
        real*8  a(n,na),        b(n,nb),        c(m)
c
        integer i,      j
c
c.... Do the work
c
        if ( m .eq. 1 ) then

            do i = 1, n
                b(i,1) = b(i,1) + c(1) * a(i,1)
            enddo

        else if ( m .eq. 2 ) then

            do i = 1, n
                b(i,1) = b(i,1) + c(1) * a(i,1)
                b(i,2) = b(i,2) + c(2) * a(i,2)
            enddo

        else if ( m .eq. 3 ) then

            do i = 1, n
                b(i,1) = b(i,1) + c(1) * a(i,1)
                b(i,2) = b(i,2) + c(2) * a(i,2)
                b(i,3) = b(i,3) + c(3) * a(i,3)
            enddo

        else if ( m .eq. 4 ) then

            do i = 1, n
                b(i,1) = b(i,1) + c(1) * a(i,1)
                b(i,2) = b(i,2) + c(2) * a(i,2)
                b(i,3) = b(i,3) + c(3) * a(i,3)
                b(i,4) = b(i,4) + c(4) * a(i,4)
            enddo

        else

            do j = 1, m 
                do i = 1, n 
                    b(i,j) = b(i,j) + c(j) * a(i,j)
                enddo
            enddo

        endif
c
        return
        end
