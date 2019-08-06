
      subroutine timeseries(blk,ycl, xl, ien, sgn)

      use timedata
      use eblock
      include "common.h"
      type (LocalBlkData) blk


      dimension shape(blk%s), ycl(bsz,blk%s,ndofl),
     &     ien(blk%e,blk%s), xl(bsz,blk%n,nsd),         
     &     sgn(blk%e,blk%s)
      real*8 al(blk%e,blk%n,nsd), 
     &     zi0(blk%e,nsd), detaij(blk%e), dzi0(blk%e,nsd),
     &     m11(blk%e), m12(blk%e), m13(blk%e), m21(blk%e), m22(blk%e),
     &     m23(blk%e), m31(blk%e), m32(blk%e), m33(blk%e), 
     &     r1(blk%e), r2(blk%e), r3(blk%e), shgradl(blk%s,nsd)
     
      real*8 xts1, xts2, xts3
      real*8 soln(ndof)
      integer e, founde

      do jj = 1, ntspts
         founde = 0
         if(statptts(jj,1).gt.0) then
            if(statptts(jj,1).eq.iblkts) then
               if(blk%l.eq.2) then ! hex
                  call shphex (blk%o, parptts(jj,:),shape(1:blk%s),
     &                 shgradl(1:blk%s,:))
               elseif(blk%l.eq.1) then
                  call shptet (blk%o, parptts(jj,:),shape(1:blk%s),
     &                 shgradl(1:blk%s,:))
               endif
               founde=statptts(jj,2)
            endif
         else
            xts1 = ptts(jj,1)
            xts2 = ptts(jj,2)
            xts3 = ptts(jj,3)

            if(blk%l.eq.2) then ! hex

               call get_a_not_hex(blk,xl,al) ! get mapping poly. coeff.
         
c...  get initial guess for Newton Iteration procedure
         
               detaij(1:blk%e) = 
     &              -al(1:blk%e,2,1)*al(1:blk%e,3,2)*al(1:blk%e,4,3)  
     &              +al(1:blk%e,2,1)*al(1:blk%e,4,2)*al(1:blk%e,3,3)
     &              + al(1:blk%e,2,2)*al(1:blk%e,3,1)*al(1:blk%e,4,3) 
     &              - al(1:blk%e,2,2)*al(1:blk%e,4,1)*al(1:blk%e,3,3) 
     &              - al(1:blk%e,2,3)*al(1:blk%e,3,1)*al(1:blk%e,4,2)
     &              +al(1:blk%e,2,3)*al(1:blk%e,4,1)*al(1:blk%e,3,2)
            
               detaij = 1./detaij
            
               zi0(1:blk%e,1) = 
     &           detaij(1:blk%e)
     &              *((al(1:blk%e,4,2)*al(1:blk%e,3,3)
     &              - al(1:blk%e,3,2)*al(1:blk%e,4,3))
     &              *(xts1-al(1:blk%e,1,1)) +
     &              (al(1:blk%e,3,1)*al(1:blk%e,4,3)
     &              - al(1:blk%e,4,1)*al(1:blk%e,3,3))
     &              *(xts2-al(1:blk%e,1,2)) +
     &              (al(1:blk%e,4,1)*al(1:blk%e,3,2)
     &              - al(1:blk%e,3,1)*al(1:blk%e,4,2))
     &              *(xts3-al(1:blk%e,1,3)))
            
            
               zi0(1:blk%e,2) = 
     &           detaij(1:blk%e)
     &              *((al(1:blk%e,2,2)*al(1:blk%e,4,3)
     &              - al(1:blk%e,4,2)*al(1:blk%e,2,3))
     &              *(xts1-al(1:blk%e,1,1)) +
     &              (al(1:blk%e,4,1)*al(1:blk%e,2,3)
     &              - al(1:blk%e,2,1)*al(1:blk%e,4,3))
     &              *(xts2-al(1:blk%e,1,2)) +
     &              (al(1:blk%e,2,1)*al(1:blk%e,4,2)
     &              - al(1:blk%e,4,1)*al(1:blk%e,2,2))
     &              *(xts3-al(1:blk%e,1,3)))
            
               zi0(1:blk%e,3) =
     &           detaij(1:blk%e)
     &              *((al(1:blk%e,3,2)*al(1:blk%e,2,3)
     &              - al(1:blk%e,2,2)*al(1:blk%e,3,3))
     &              *(xts1-al(1:blk%e,1,1)) +
     &              (al(1:blk%e,2,1)*al(1:blk%e,3,3)
     &              - al(1:blk%e,3,1)*al(1:blk%e,2,3))
     &              *(xts2-al(1:blk%e,1,2)) +
     &              (al(1:blk%e,3,1)*al(1:blk%e,2,2)
     &              - al(1:blk%e,2,1)*al(1:blk%e,3,2))
     &              *(xts3-al(1:blk%e,1,3)))
            
            
c...  iterate to convergence
            
               do it = 1, iterat
               
c...  build matrix
               
                  m11(1:blk%e)=
     &              al(1:blk%e,2,1)+al(1:blk%e,5,1)*zi0(1:blk%e,2)
     &              +al(1:blk%e,7,1)*zi0(1:blk%e,3)
     &                 +al(1:blk%e,8,1)*zi0(1:blk%e,2)*zi0(1:blk%e,3)
                  m12(1:blk%e)=
     &              al(1:blk%e,3,1)+al(1:blk%e,5,1)*zi0(1:blk%e,1)
     &                  +al(1:blk%e,6,1)*zi0(1:blk%e,3)
     &                 +al(1:blk%e,8,1)*zi0(1:blk%e,1)*zi0(1:blk%e,3)
                  m13(1:blk%e)=
     &              al(1:blk%e,4,1)+al(1:blk%e,6,1)*zi0(1:blk%e,2)
     &                  +al(1:blk%e,7,1)*zi0(1:blk%e,1)
     &                 +al(1:blk%e,8,1)*zi0(1:blk%e,1)*zi0(1:blk%e,2)

                  m21(1:blk%e)=
     &              al(1:blk%e,2,2)+al(1:blk%e,5,2)*zi0(1:blk%e,2)
     &                  +al(1:blk%e,7,2)*zi0(1:blk%e,3)
     &                 +al(1:blk%e,8,2)*zi0(1:blk%e,2)*zi0(1:blk%e,3)
                  m22(1:blk%e)=
     &              al(1:blk%e,3,2)+al(1:blk%e,5,2)*zi0(1:blk%e,1)
     &                  +al(1:blk%e,6,2)*zi0(1:blk%e,3)
     &                 +al(1:blk%e,8,2)*zi0(1:blk%e,1)*zi0(1:blk%e,3)
                  m23(1:blk%e)=
     &              al(1:blk%e,4,2)+al(1:blk%e,6,2)*zi0(1:blk%e,2)
     &                  +al(1:blk%e,7,2)*zi0(1:blk%e,1)
     &                 +al(1:blk%e,8,2)*zi0(1:blk%e,1)*zi0(1:blk%e,2)

                  m31(1:blk%e)=
     &              al(1:blk%e,2,3)+al(1:blk%e,5,3)*zi0(1:blk%e,2)
     &                  +al(1:blk%e,7,3)*zi0(1:blk%e,3)
     &                 +al(1:blk%e,8,3)*zi0(1:blk%e,2)*zi0(1:blk%e,3)
                  m32(1:blk%e)=
     &              al(1:blk%e,3,3)+al(1:blk%e,5,3)*zi0(1:blk%e,1)
     &                  +al(1:blk%e,6,3)*zi0(1:blk%e,3)
     &                 +al(1:blk%e,8,3)*zi0(1:blk%e,1)*zi0(1:blk%e,3)
                  m33(1:blk%e)=
     &              al(1:blk%e,4,3)+al(1:blk%e,6,3)*zi0(1:blk%e,2)
     &                  +al(1:blk%e,7,3)*zi0(1:blk%e,1)
     &                 +al(1:blk%e,8,3)*zi0(1:blk%e,1)*zi0(1:blk%e,2)
               
               
c...  build rhs
               
                  r1(1:blk%e)=
     &              al(1:blk%e,1,1)+al(1:blk%e,2,1)*zi0(1:blk%e,1)
     &              +al(1:blk%e,3,1)*zi0(1:blk%e,2)+
     &                 al(1:blk%e,4,1)*zi0(1:blk%e,3)
     &              +al(1:blk%e,5,1)*zi0(1:blk%e,1)*zi0(1:blk%e,2)+
     &                 al(1:blk%e,6,1)*zi0(1:blk%e,2)*zi0(1:blk%e,3)
     &              +al(1:blk%e,7,1)*
     &                 zi0(1:blk%e,1)*zi0(1:blk%e,3)
     &              +al(1:blk%e,8,1)*zi0(1:blk%e,1)*
     &                 zi0(1:blk%e,2)*zi0(1:blk%e,3) - xts1

                  r2(1:blk%e)=
     &              al(1:blk%e,1,2)+al(1:blk%e,2,2)*zi0(1:blk%e,1)
     &              +al(1:blk%e,3,2)*zi0(1:blk%e,2)+
     &                 al(1:blk%e,4,2)*zi0(1:blk%e,3)
     &              +al(1:blk%e,5,2)*zi0(1:blk%e,1)*zi0(1:blk%e,2)+
     &                 al(1:blk%e,6,2)*zi0(1:blk%e,2)*zi0(1:blk%e,3)
     &              +al(1:blk%e,7,2)*
     &                 zi0(1:blk%e,1)*zi0(1:blk%e,3)
     &              +al(1:blk%e,8,2)*zi0(1:blk%e,1)*
     &                 zi0(1:blk%e,2)*zi0(1:blk%e,3) - xts2
               
                  r3(1:blk%e)=
     &              al(1:blk%e,1,3)+al(1:blk%e,2,3)*zi0(1:blk%e,1)+
     &              al(1:blk%e,3,3)*zi0(1:blk%e,2)+
     &                 al(1:blk%e,4,3)*zi0(1:blk%e,3)
     &              +al(1:blk%e,5,3)*zi0(1:blk%e,1)*zi0(1:blk%e,2)+
     &                 al(1:blk%e,6,3)*zi0(1:blk%e,2)*zi0(1:blk%e,3)
     &              +al(1:blk%e,7,3)*
     &                 zi0(1:blk%e,1)*zi0(1:blk%e,3)
     &              +al(1:blk%e,8,3)*zi0(1:blk%e,1)*
     &                 zi0(1:blk%e,2)*zi0(1:blk%e,3) - xts3

c...  get solution

                  detaij = m11*m22*m33-m11*m23*m32-m21*m12*m33+
     &                 m21*m13*m32+m31*m12*m23-m31*m13*m22  

                  detaij = 1./detaij

                  dzi0(1:blk%e,1) = 
     &              -detaij(1:blk%e)*((m22(1:blk%e)*m33(1:blk%e)
     &              -m23(1:blk%e)*m32(1:blk%e))*
     &                 r1(1:blk%e) + (m13(1:blk%e)*m32(1:blk%e)
     &              -m12(1:blk%e)*m33(1:blk%e))*r2(1:blk%e) + 
     &                 (m12(1:blk%e)*m23(1:blk%e)
     &              -m13(1:blk%e)*m22(1:blk%e))*r3(1:blk%e))
                  dzi0(1:blk%e,2) = 
     &              -detaij(1:blk%e)*((m23(1:blk%e)*m31(1:blk%e)
     &              -m21(1:blk%e)*m32(1:blk%e))*
     &                 r1(1:blk%e) + (m11(1:blk%e)*m33(1:blk%e)
     &              -m13(1:blk%e)*m31(1:blk%e))*r2(1:blk%e) + 
     &                 (m13(1:blk%e)*m21(1:blk%e)
     &              -m11(1:blk%e)*m23(1:blk%e))*r3(1:blk%e))
                  dzi0(1:blk%e,3) =
     &              -detaij(1:blk%e)*((m21(1:blk%e)*m32(1:blk%e)
     &              -m22(1:blk%e)*m31(1:blk%e))*
     &                 r1(1:blk%e) + (m12(1:blk%e)*m31(1:blk%e)
     &              -m11(1:blk%e)*m32(1:blk%e))*r2(1:blk%e) + 
     &                 (m11(1:blk%e)*m22(1:blk%e)
     &              -m12(1:blk%e)*m21(1:blk%e))*r3(1:blk%e))

                  zi0(1:blk%e,:) = zi0(1:blk%e,:) + dzi0(1:blk%e,:)
               
               enddo
            
               do e = 1, blk%e
                  if ((abs(zi0(e,1)).lt.(one+tolpt)).and.
     &                 (abs(zi0(e,2)).lt.(one+tolpt)).and.
     &                 (abs(zi0(e,3)).lt.(one+tolpt))) then ! got the element

                     call shphex (blk%o, zi0(e,:),shape(1:blk%s),
     &                    shgradl(1:blk%s,:))
                  
                     founde=e
                     exit
                  endif
            
               enddo
            elseif (blk%l.eq.1) then !tet

               call get_a_not_tet(blk,xl,al)
                       
c
c solve for r, s, t  for each elements
c 
               do e = 1, blk%e
                  detaij(e) = al(e,2,1)*(-al(e,3,2)*al(e,4,3) + 
     &                 al(e,4,2)*al(e,3,3)) + al(e,2,2)*
     &                 (al(e,3,1)*al(e,4,3) - al(e,4,1)*
     &                 al(e,3,3)) + al(e,2,3)*(-al(e,3,1)*al(e,4,2)+
     &                 al(e,4,1)*al(e,3,2))
            
                  detaij(e) = 1./detaij(e)

                  zi0(e,1) = detaij(e)*((al(e,4,2)*al(e,3,3)
     &                 - al(e,3,2)*al(e,4,3))*(xts1-al(e,1,1)) +
     &                 (al(e,3,1)*al(e,4,3)
     &                 - al(e,4,1)*al(e,3,3))*(xts2-al(e,1,2)) +
     &                 (al(e,4,1)*al(e,3,2)
     &                 - al(e,3,1)*al(e,4,2))*(xts3-al(e,1,3)))

                  zi0(e,2) = detaij(e)*((al(e,2,2)*al(e,4,3)
     &                 - al(e,4,2)*al(e,2,3))*(xts1-al(e,1,1)) +
     &                 (al(e,4,1)*al(e,2,3)
     &                 - al(e,2,1)*al(e,4,3))*(xts2-al(e,1,2)) +
     &                 (al(e,2,1)*al(e,4,2)
     &                 - al(e,4,1)*al(e,2,2))*(xts3-al(e,1,3)))

                  zi0(e,3) = detaij(e)*((al(e,3,2)*al(e,2,3)
     &                 - al(e,2,2)*al(e,3,3))*(xts1-al(e,1,1)) +
     &                 (al(e,2,1)*al(e,3,3)
     &                 - al(e,3,1)*al(e,2,3))*(xts2-al(e,1,2)) +
     &                 (al(e,3,1)*al(e,2,2)
     &                 - al(e,2,1)*al(e,3,2))*(xts3-al(e,1,3)))
               
                  if ((zi0(e,1)+zi0(e,2)+zi0(e,3)).lt.(one+tolpt).and.
     &                 zi0(e,1).lt.(one +tolpt).and.    !should not be necessary; the limit in the tet is already defined by the previous line
     &                 zi0(e,1).gt.(zero-tolpt).and.
     &                 zi0(e,2).lt.(one +tolpt).and.    !should not be necessary 
     &                 zi0(e,2).gt.(zero-tolpt).and.
     &                 zi0(e,3).lt.(one +tolpt).and.    !should not be necessary 
     &                 zi0(e,3).gt.(zero-tolpt)) then
                     
                     call shptet (blk%o, zi0(e,:), shape(1:blk%s),
     &                    shgradl(1:blk%s,:))
                     if(jj.eq.2) then
                       write(*,*) 'second point found: rst=', zi0(e,1),zi0(e,2),zi0(e,3)
                       write(*,*) 'xyz_nd1=', xl(e,1,1),xl(e,1,2),xl(e,1,3)
                       write(*,*) 'soln_nd1=', ycl(e,1,1),ycl(e,1,2),ycl(e,1,3),ycl(e,1,4)
                     endif
                     founde=e
                     exit
                  endif
               enddo
            endif

            if(founde.ne.0) then !store values for next time steps
               statptts(jj,1)=iblkts
               statptts(jj,2)=founde
               parptts(jj,:)=zi0(founde,:)
            else
              if(iblkts.eq.nelblk) then !not found in last elm-blk of this part
                 statptts(jj,1)=nelblk+1 ! not searched for following steps
              endif
            endif
         endif

         if(founde.ne.0) then
            soln(1:ndofl) = zero

            do i = 1,blk%n
               soln(1:ndofl) = soln(1:ndofl)
     &              +ycl(founde,i,1:ndofl)*shape(i)
            enddo
            do i = 1+blk%n,blk%s
               soln(1:ndofl) = soln(1:ndofl)
     &              +ycl(founde,i,1:ndofl)*shape(i)*sgn(founde,i)
            enddo

            varts(jj,:) = soln(:)
         endif
         
      enddo

      return 
      end
      
