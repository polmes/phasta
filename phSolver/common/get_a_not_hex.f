

      subroutine get_a_not_hex(blk,xc,anot)
    
      use eblock
      include "common.h"
      type (LocalBlkData) blk  

      dimension xc(blk%e,blk%n,nsd), anot(blk%e,blk%n,nsd)


      do i = 1, nsd

         anot(1:blk%e,1,i) = 
     &      pt125*(xc(1:blk%e,1,i)
     &        +xc(1:blk%e,2,i)+xc(1:blk%e,3,i)+xc(1:blk%e,4,i)
     &        +xc(1:blk%e,5,i)+xc(1:blk%e,6,i)+xc(1:blk%e,7,i)
     &        +xc(1:blk%e,8,i))
         
         anot(1:blk%e,2,i) = 
     &      pt125*(-xc(1:blk%e,1,i)
     &        +xc(1:blk%e,2,i)+xc(1:blk%e,3,i)-xc(1:blk%e,4,i)
     &        -xc(1:blk%e,5,i)+xc(1:blk%e,6,i)+xc(1:blk%e,7,i)
     &        -xc(1:blk%e,8,i))

         anot(1:blk%e,3,i) = 
     &      pt125*(-xc(1:blk%e,1,i)
     &        -xc(1:blk%e,2,i)+xc(1:blk%e,3,i)+xc(1:blk%e,4,i)
     &        -xc(1:blk%e,5,i)-xc(1:blk%e,6,i)+xc(1:blk%e,7,i)
     &        +xc(1:blk%e,8,i))
         
         anot(1:blk%e,4,i) = 
     &      pt125*(-xc(1:blk%e,1,i)
     &        -xc(1:blk%e,2,i)-xc(1:blk%e,3,i)-xc(1:blk%e,4,i)
     &        +xc(1:blk%e,5,i)+xc(1:blk%e,6,i)+xc(1:blk%e,7,i)
     &        +xc(1:blk%e,8,i))

         anot(1:blk%e,5,i) = 
     &      pt125*(xc(1:blk%e,1,i)
     &        -xc(1:blk%e,2,i)+xc(1:blk%e,3,i)-xc(1:blk%e,4,i)
     &        +xc(1:blk%e,5,i)-xc(1:blk%e,6,i)+xc(1:blk%e,7,i)
     &        -xc(1:blk%e,8,i))

         anot(1:blk%e,6,i) = 
     &      pt125*(xc(1:blk%e,1,i)
     &        +xc(1:blk%e,2,i)-xc(1:blk%e,3,i)-xc(1:blk%e,4,i)
     &        -xc(1:blk%e,5,i)-xc(1:blk%e,6,i)+xc(1:blk%e,7,i)
     &        +xc(1:blk%e,8,i))

         anot(1:blk%e,7,i) = 
     &      pt125*(xc(1:blk%e,1,i)
     &        -xc(1:blk%e,2,i)-xc(1:blk%e,3,i)+xc(1:blk%e,4,i)
     &        -xc(1:blk%e,5,i)+xc(1:blk%e,6,i)+xc(1:blk%e,7,i)
     &        -xc(1:blk%e,8,i))

         anot(1:blk%e,8,i) = 
     &      pt125*(-xc(1:blk%e,1,i)
     &        +xc(1:blk%e,2,i)-xc(1:blk%e,3,i)+xc(1:blk%e,4,i)
     &        +xc(1:blk%e,5,i)-xc(1:blk%e,6,i)+xc(1:blk%e,7,i)
     &        -xc(1:blk%e,8,i))

      enddo

      return
      end


      subroutine get_a_not_tet(blk,xc,anot)
        
      use eblock
      include "common.h"
      type (LocalBlkData) blk
      dimension xc(blk%e,blk%n,nsd), anot(blk%e,blk%n,nsd)


      do i = 1, nsd

         anot(1:blk%e,1,i) = xc(1:blk%e,4,i)
         anot(1:blk%e,2,i) = xc(1:blk%e,1,i)-xc(1:blk%e,4,i)
         anot(1:blk%e,3,i) = xc(1:blk%e,2,i)-xc(1:blk%e,4,i)
         anot(1:blk%e,4,i) = xc(1:blk%e,3,i)-xc(1:blk%e,4,i)
         
      enddo
      
      return
      end
      
