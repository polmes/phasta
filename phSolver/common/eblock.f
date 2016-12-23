      module eblock
      type LocalBlkData
          integer :: n,s,e,g,l,o,b,t,i
!n=nenl, s=nshl, e=nelInThisBlock=npro, g=ngauss, l=lcyst, o=order,b=iblk, t=ith, i=iel (first element number of the block)
      end type LocalBlkData


      type LocalBlkDataB
          integer :: ni,si,nb,sb,e,g,l,o
!i=interior, b=boundary
!n=nenl, s=nshl, e=nelInThisBlock=npro, g=ngauss, l=lcyst, o=order
      end type LocalBlkDataB

! if at some point in the future you are tempted to define blk here and remove it from all of the subroutine interfaces first 
! check to be sure that the fortran 90 standard allows threads to have private copies of module data....apparantly at this time
! they do not so putting it here would kill threading

      end module
