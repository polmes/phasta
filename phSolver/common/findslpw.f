       subroutine findslpw(x,ilwork,iper,iBC)
       use pointer_data    
       use slpw       
       include "common.h"
       include "mpif.h"
       real*8 x(numnp,nsd)
       integer iper(nshg)
       integer ilwork(nlwork)
       integer ifirstvisit(nshg)
       integer, allocatable :: ienb(:)
       real*8 elnrm(nsd),wn1(nsd),wn2(nsd),wn3(nsd),aedg(nsd),bedg(nsd) 
       real*8 mag
       integer iBC(nshg)
       integer, allocatable :: IDslpw(:)
       real*8 AA(nsd),BB(nsd),tmp
       integer iparall

       allocate(idxslpw(nshg))
       allocate(mpslpw(nshg,2))
       allocate(wlnorm(nshg,nsd,9))
             
              
      nslpwnd=0 
      idxslpw=0
      mpslpw=0
      wlnorm=0
      ifirstvisit(:)=1
      do iblk=1,nelblb   ! loop over boundary element blocks
        npro=lcblkb(1,iblk+1)-lcblkb(1,iblk)   ! number of elements in blk
        nenbl=lcblkb(6,iblk) !boundary nodes 
        nshl=lcblkb(9,iblk) ! elment nodes
        allocate(ienb(nshl))
        do i=1,npro 
          isfID=miBCB(iblk)%p(i,2) ! surfID for this element
          if(isfID.gt.709.or.isfID.lt.701)cycle ! skip out if not in slip range
          islpw=mod(isfID,700) ! single digit index for different surfIDs
          ienb(1:nshl)=mienb(iblk)%p(i,1:nshl) !extract node numbers
          do j=1,nenbl ! loop over the only the nodes on the bdy
            nn=ienb(j) ! current node number
            if(ifirstvisit(nn).eq.1)then  !array to track first visit
              ifirstvisit(nn)=0  !mark to zero means visited
              nslpwnd=nslpwnd+1  !track unique nodes treated
              idxslpw(nslpwnd)=nn ! create a map(1:nslpwnd)=global node #
            endif
!
! mslpw is a key data structure: This node's  (nn)
! mslpw(nn,2) is bit tested with the single digit surfID
! because if it is not set, the conditional is going to set it to one 
! AND it is going to add one to the mslpw(nn,1) for this node
! if this node is on a model edge it will get inside of this conditional
! once for each surfID and thus mslpw is counting the number of times
! this node computes a normal (with the help of mslpw(nn,2)
!
            if(.not.btest(mpslpw(nn,2),islpw))then
              mpslpw(nn,1)=mpslpw(nn,1)+1   
!
!set the single digit surfID bit of the second rank to 1 so that we do not! count repeats on the same for the boundary elements with the same surfID
!
              mpslpw(nn,2)=mpslpw(nn,2)+2**islpw  
            endif
!
! rather complicated way to get edge vectors depending on boundary element
! node number
!
            if(j.ne.1.and.j.ne.nenbl)then
              aedg(:)=x(ienb(j+1),:)-x(nn,:)
              bedg(:)=x(ienb(j-1),:)-x(nn,:)
            elseif(j.eq.1)then
              aedg(:)=x(ienb(j+1),:)-x(nn,:)
              bedg(:)=x(ienb(nenbl),:)-x(nn,:)
            elseif(j.eq.nenbl)then
              aedg(:)=x(ienb(1),:)-x(nn,:)
              bedg(:)=x(ienb(j-1),:)-x(nn,:)   
            endif  
! perform cross product
            elnrm(1)=aedg(2)*bedg(3)-aedg(3)*bedg(2)
            elnrm(2)=aedg(3)*bedg(1)-aedg(1)*bedg(3)
            elnrm(3)=aedg(1)*bedg(2)-aedg(2)*bedg(1)
!
! add this normal to any existing ones for this node
! note it will be re-normalized later and thus this is averaging
!
            wlnorm(nn,:,islpw)=wlnorm(nn,:,islpw)+elnrm(:)
          enddo
        enddo
        deallocate(ienb)
      enddo
      
       do k=1,9
         if(numpe.gt.1) call commu(wlnorm(:,:,k),ilwork,3,'in ') 
         call bc3per(iBC,wlnorm(:,:,k),iper,ilwork,3)
         if(numpe.gt.1) call commu(wlnorm(:,:,k),ilwork,3,'out')
       enddo
cccccccccccccccccccccccccccccccccccccccccccccccccc   
      
     
      do id=1,nslpwnd
        nn=idxslpw(id)
        nslpw=mpslpw(nn,1)
        ihis=mpslpw(nn,2)
        nslpwg=nslpw
 
        allocate(IDslpw(nslpw))
        jslpw=0 
        do ni=1,nslpw
          do islpw=1+jslpw,9
             if(btest(ihis,islpw))exit 
          enddo
          IDslpw(ni)=islpw
          jslpw=islpw 
       enddo  

       ired=0   
       do i=1,nslpwg
          if(btest(ired,i))cycle 
       do j=i+1,nslpwg
          AA(:)=wlnorm(nn,:,IDslpw(i))          
          BB(:)=wlnorm(nn,:,IDslpw(j))
          if(iparall(AA,BB).eq.1)then
             nslpw=nslpw-1
             ihis=ihis-2**IDslpw(j)
             wlnorm(nn,:,IDslpw(i))=AA(:)+BB(:) 
             ired=ired+2**j  
          else 
!
! what follows before is expected to work for reconcilliation of the normal 
! to constrain when two model surfaces with unique surfID's have normal 
! vector differences larger than angtol.  Not tested for more than 2 (e.g.,
! the corners of a cube have 3 faces but not really sure I want to try to 
! get an inviscid solution there anyway...could be an issue for wall modeled
! turbulence but, again, not sure I expect a physical prediction there 
! either).
!
! if the two normals are not parallel then we are on a model 
! edge and we want to allow flow along the edge AND along a vector that
! is the average of the two normals that are constrained on the mesh faces that share this edge
! BUT remove/constrain-to-zero flow perpendicular to those two free components.
! If we approximate the edge normal as the average of the two 
! face normals AND find the normal to that plane as the free edge vector, our desired vector
! will just be the cross of those 2 vectors
!
! make A and B unit vectors
             mag=(AA(1)**2+AA(2)**2+AA(3)**2)**0.5
             if(mag.ne.0)AA(:)=AA(:)/mag             
             mag=(BB(1)**2+BB(2)**2+BB(3)**2)**0.5
             if(mag.ne.0)BB(:)=BB(:)/mag             
! perform cross product to get the normal to n_A and n_B
            elnrm(1)=AA(2)*BB(3)-AA(3)*BB(2)
            elnrm(2)=AA(3)*BB(1)-AA(1)*BB(3)
            elnrm(3)=AA(1)*BB(2)-AA(2)*BB(1)
! build edge normal vector
             AA=AA+BB
! normalize the edge normal vector
             mag=(AA(1)**2+AA(2)**2+AA(3)**2)**0.5
             if(mag.ne.0)AA(:)=AA(:)/mag             
! normalize the edge slip vector
             mag=(elnrm(1)**2+elnrm(2)**2+elnrm(3)**2)**0.5
             if(mag.ne.0)elnrm(:)=elnrm(:)/mag             
! perform cross product to get the vector we want to constrain to zero
            wlnorm(nn,1,IDslpw(i))=AA(2)*elnrm(3)-AA(3)*elnrm(2)
            wlnorm(nn,2,IDslpw(i))=AA(3)*elnrm(1)-AA(1)*elnrm(3)
            wlnorm(nn,3,IDslpw(i))=AA(1)*elnrm(2)-AA(2)*elnrm(1)

             nslpw=nslpw-1
             ihis=ihis-2**IDslpw(j)
             ired=ired+2**j  
          endif
       enddo
       enddo            

       do k=1,9
          wn1(:)=wlnorm(nn,:,k)
          mag=(wn1(1)**2+wn1(2)**2+wn1(3)**2)**0.5
          if(mag.ne.0)wlnorm(nn,:,k)=wn1(:)/mag             
       enddo
            
      mpslpw(nn,1)=nslpw
      mpslpw(nn,2)=ihis
 
      deallocate(IDslpw) 
      enddo
    

c        call write_field(myrank,'a'//char(0),
c     &      'wnrm'//char(0),4,wlnorm(:,:,1),'d'//char(0),nshg,3,701) 
cc 7 is the number of charaters of the name of array wnrm701
c        call write_field(myrank,'a'//char(0),
c     &      'wnrm'//char(0),4,wlnorm(:,:,2),'d'//char(0),nshg,3,702)
c      stop
 
      return
      end

cccccccccccccccccccccc

      integer function iparall(A,B)
      real*8 A(3),B(3)
      real*8 angtol,tmp
 
      tmp=(A(1)**2+A(2)**2+A(3)**2)**0.5
      A(:)=A(:)/tmp
      tmp=(B(1)**2+B(2)**2+B(3)**2)**0.5
      B(:)=B(:)/tmp
      pi=3.1415926535
      angtol=10
      angtol=sin(angtol*(pi/180))
      C1=A(2)*B(3)-A(3)*B(2)
      C2=A(3)*B(1)-A(1)*B(3)   
      C3=A(1)*B(2)-A(2)*B(1)
      Cmag=(C1**2+C2**2+C3**2)**0.5
      dotp=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
      if(dotp.gt.0.and.Cmag.lt.angtol)then
         iparall=1
      else
         iparall=0
      endif

      return
      end 
