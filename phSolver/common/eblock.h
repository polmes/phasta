      type LocalBlkData
          integer :: n,s,e,g,l,o,b,t,i
!n=nenl, s=nshl, e=nelInThisBlock=npro, g=ngauss, l=lcyst, o=order,b=iblk, t=ith, i=iel (first element number of the block)
      end type LocalBlkData


      type LocalBlkDataB
          integer :: ni,si,nb,sb,e,g,l,o
!i=interior, b=boundary
!n=nenl, s=nshl, e=nelInThisBlock=npro, g=ngauss, l=lcyst, o=order
      end type LocalBlkDataB

