        subroutine hessian ( y,         x,     
     &                       shp,       shgl,      iBC,
     &                       shpb,      shglb,     iper,      
     &                       ilwork,    uhess,     gradu  )
      use pointer_data  ! brings in the pointers for the blocked arrays
      use eblock

      include "common.h"

c
        dimension y(nshg,ndof),      
     &            x(numnp,nsd),         iBC(nshg),           
     &            iper(nshg)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
        dimension gradu(nshg,9),     rmass(nshg),
     &            uhess(nshg,27)
c
        dimension ilwork(nlwork)

c
           gradu = zero
           rmass = zero
        
           do iblk = 1, nelblk
              iel    = lcblk(1,iblk)
              lelCat = lcblk(2,iblk)
              lcsyst = lcblk(3,iblk)
              iorder = lcblk(4,iblk)
              nenl   = lcblk(5,iblk) ! no. of vertices per element
              nshl   = lcblk(10,iblk)
              mattyp = lcblk(7,iblk)
              ndofl  = lcblk(8,iblk)
              nsymdl = lcblk(9,iblk)
              npro   = lcblk(1,iblk+1) - iel 
              ngauss = nint(lcsyst)
          blk%n   = lcblk(5,iblk) ! no. of vertices per element
          blk%s   = lcblk(10,iblk)
          blk%e   = lcblk(1,iblk+1) - iel
          blk%g = nint(lcsyst)
          blk%l = lcblk(3,iblk)
          blk%o = lcblk(4,iblk)

c     
              call velocity_gradient ( blk, y,                
     &                                 x,                       
     &                                 shp(lcsyst,1:nshl,:), 
     &                                 shgl(lcsyst,:,1:nshl,:),
     &                                 mien(iblk)%p,     
     &                                 gradu, 
     &                                 rmass )

           end do

c
           call reconstruct( rmass, gradu, iBC, iper, ilwork, 9 )       

           uhess = zero
           rmass = zero
        
           do iblk = 1, nelblk
              iel    = lcblk(1,iblk)
              lelCat = lcblk(2,iblk)
              lcsyst = lcblk(3,iblk)
              iorder = lcblk(4,iblk)
              nenl   = lcblk(5,iblk) ! no. of vertices per element
              nshl   = lcblk(10,iblk)
              mattyp = lcblk(7,iblk)
              ndofl  = lcblk(8,iblk)
              nsymdl = lcblk(9,iblk)
              npro   = lcblk(1,iblk+1) - iel 
              ngauss = nint(lcsyst)
          blk%n   = lcblk(5,iblk) ! no. of vertices per element
          blk%s   = lcblk(10,iblk)
          blk%e   = lcblk(1,iblk+1) - iel
          blk%g = nint(lcsyst)
          blk%l = lcblk(3,iblk)
          blk%o = lcblk(4,iblk)
c     

              call velocity_hessian (  blk,gradu,                
     &                                 x,                       
     &                                 shp(lcsyst,1:nshl,:), 
     &                                 shgl(lcsyst,:,1:nshl,:),
     &                                 mien(iblk)%p,     
     &                                 uhess,
     &				       rmass  )    
           end do
       

           call reconstruct( rmass, uhess, iBC, iper, ilwork, 27 )       
c
      return
      end

c-----------------------------------------------------------------------------

        subroutine velocity_gradient (blk, y,       x,       shp,     shgl, 
     &                                 ien,     gradu,   rmass    )

      use eblock
      include "common.h"
c
        dimension y(nshg,ndof),               x(numnp,nsd),            
     &            shp(nshl,ngauss),           shgl(nsd,nshl,ngauss),
     &            ien(npro,nshl),             gradu(nshg,9), 
     &            shdrv(npro,nsd,nshl),       shape( npro, nshl ),      
     &            gradul(npro,9) ,            rmass( nshg ) 
c
        dimension yl(bsz,nshl,ndof),          xl(bsz,nenl,nsd),
     &            ql(bsz,nshl,9),             dxidx(npro,nsd,nsd),
     &            WdetJ(npro),		       rmassl(bsz,nshl)
c
c
        dimension sgn(npro,nshl)
c
c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
        if(ipord.gt.1) then
        do i=1,nshl
           where ( ien(:,i) < 0 )
              sgn(:,i) = -one
           elsewhere
              sgn(:,i) = one
           endwhere
        enddo
        endif

c
c.... gather the variables
c

        call localy (blk,y,    yl,     ien,    ndof,   'gather  ')
        call localx (blk,x,    xl,     ien,    nsd,    'gather  ')
c
c.... get the element residuals 
c
        ql     = zero
        rmassl = zero

        do intp = 1, ngauss

            if ( Qwt( lcsyst, intp ) .eq. zero ) cycle

            gradul = zero
            call getshp(blk,intp, shp, shgl, sgn, shape, shdrv )
            call local_gradient(blk, yl(:,:,2:4), 3,  shdrv, xl, 
     &                           gradul , dxidx, WdetJ )

c.... assemble contribution of gradu to each element node
c     
            do i=1,nshl
                do j = 1, 9
                    ql(1:blk%e,i,j) = ql(1:blk%e,i,j)+shape(:,i)*WdetJ*gradul(:,j)
                end do
                
                rmassl(1:blk%e,i) = rmassl(1:blk%e,i) + shape(:,i)*WdetJ

             end do

        end do
c
c
        call local (blk,gradu,  ql,     ien,  9,  'scatter ')
        call local (blk,rmass,  rmassl, ien,  1,  'scatter ')
c
c.... end
c
        return
        end


c-----------------------------------------------------------------------------

        subroutine velocity_hessian ( blk,gradu,   x,     shp,   shgl, 
     &                                ien,     uhess, rmass  )

      use eblock
      include "common.h"
c
        dimension gradu(nshg,9),              x(numnp,nsd),            
     &            shp(nshl,ngauss),           shgl(nsd,nshl,ngauss),
     &            ien(npro,nshl),             uhess(nshg,27), 
     &            shdrv(npro,nsd,nshl),       shape( npro, nshl ),
     &            uhessl(bsz,27),            rmass( nshg ) 
c
        dimension gradul(npro,nshl,9),          xl(bsz,nenl,nsd),         
     &            ql(bsz,nshl,27),             dxidx(npro,nsd,nsd),    
     &            WdetJ(npro),                  rmassl(npro, nshl)
c
c
        dimension sgn(npro,nshl)
c
c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
       if (ipord.gt.1) then
        do i=1,nshl
           where ( ien(:,i) < 0 )
              sgn(:,i) = -one
           elsewhere
              sgn(:,i) = one
           endwhere
        enddo
        endif

c
c.... gather the variables
c

        call local  (blk,gradu,  gradul, ien,    9 ,   'gather  ')
        call localx (blk,x,      xl,     ien,    nsd,  'gather  ')
c
c.... get the element residuals 
c
        ql     = zero
	rmassl = zero

        do intp = 1, ngauss

            if ( Qwt( lcsyst, intp ) .eq. zero ) cycle

            uhessl = zero
            call getshp(blk,intp, shp, shgl, sgn, shape, shdrv )
            call local_gradient(blk, gradul, 9,  shdrv, xl, 
     &                           uhessl , dxidx, WdetJ )

c.... assemble contribution of gradu .,
c     
            do i=1,nshl
                do j = 1,27 
                    ql(1:blk%e,i,j)=ql(1:blk%e,i,j)+shape(:,i)*WdetJ*uhessl(:,j )
                end do

                rmassl(1:blk%e,i) = rmassl(1:blk%e,i) + shape(:,i)*WdetJ
             end do

        end do
c
c
        call local (blk,uhess,  ql,     ien,  27,     'scatter ')
        call local (blk,rmass,  rmassl, ien,   1,     'scatter ')
c
c.... end
c
        return
        end


c--------------------------------------------------------------------
      subroutine reconstruct( rmass, qres, iBC, iper, ilwork, vsize )

      include "common.h"
      
      integer vsize
      dimension rmass(nshg), qres( nshg, vsize),
     &          iBC(nshg), iper(nshg), ilwork(nlwork)
c
c
c.... compute qi for node A, i.e., qres <-- qres/rmass
c
       if (numpe > 1) then
          call commu ( qres  , ilwork,  vsize  , 'in ')
          call commu ( rmass , ilwork,  1  , 'in ')
       endif
c
c  take care of periodic boundary conditions
c
        do j= 1,nshg
          if ((btest(iBC(j),10))) then
            i = iper(j)
            rmass(i) = rmass(i) + rmass(j)
            qres(i,:) = qres(i,:) + qres(j,:)
          endif
        enddo

        do j= 1,nshg
          if ((btest(iBC(j),10))) then
            i = iper(j)
            rmass(j) = rmass(i)
            qres(j,:) = qres(i,:)
          endif
        enddo
c
c.... invert the diagonal mass matrix and find q
c
        rmass = one/rmass
       
       do i=1,vsize 
          qres(:,i) = rmass*qres(:,i)
       enddo

       if(numpe > 1) then
          call commu (qres, ilwork, vsize, 'out')    
       endif

c.... return
c    
        return
        end

c-------------------------------------------------------------------------

        subroutine local_gradient (blk, vector,   vsize, shgl,   xl, 
     &                              gradient, dxidx,   WdetJ )
c
c
      use eblock
      include "common.h"
c
c  passed arrays

        integer vsize
c
        dimension vector(npro,nshl,vsize), 
     &            shgl(npro,nsd,nshl),        xl(bsz,nenl,nsd),
     &            gradient(npro,vsize*3),     shg(npro,nshl,nsd), 
     &            dxidx(npro,nsd,nsd),        WdetJ(npro)
c
c  local arrays
c
        dimension tmp(npro),           dxdxi(npro,nsd,nsd)

c
c.... compute the deformation gradient
c
        dxdxi = zero
c
          do n = 1, nenl
            dxdxi(:,1,1) = dxdxi(:,1,1) + xl(1:blk%e,n,1) * shgl(:,1,n)
            dxdxi(:,1,2) = dxdxi(:,1,2) + xl(1:blk%e,n,1) * shgl(:,2,n)
            dxdxi(:,1,3) = dxdxi(:,1,3) + xl(1:blk%e,n,1) * shgl(:,3,n)
            dxdxi(:,2,1) = dxdxi(:,2,1) + xl(1:blk%e,n,2) * shgl(:,1,n)
            dxdxi(:,2,2) = dxdxi(:,2,2) + xl(1:blk%e,n,2) * shgl(:,2,n)
            dxdxi(:,2,3) = dxdxi(:,2,3) + xl(1:blk%e,n,2) * shgl(:,3,n)
            dxdxi(:,3,1) = dxdxi(:,3,1) + xl(1:blk%e,n,3) * shgl(:,1,n)
            dxdxi(:,3,2) = dxdxi(:,3,2) + xl(1:blk%e,n,3) * shgl(:,2,n)
            dxdxi(:,3,3) = dxdxi(:,3,3) + xl(1:blk%e,n,3) * shgl(:,3,n)
          enddo
c
c.... compute the inverse of deformation gradient
c
        dxidx(:,1,1) =   dxdxi(:,2,2) * dxdxi(:,3,3) 
     &                 - dxdxi(:,3,2) * dxdxi(:,2,3)
        dxidx(:,1,2) =   dxdxi(:,3,2) * dxdxi(:,1,3) 
     &                 - dxdxi(:,1,2) * dxdxi(:,3,3)
        dxidx(:,1,3) =   dxdxi(:,1,2) * dxdxi(:,2,3) 
     &                 - dxdxi(:,1,3) * dxdxi(:,2,2)
        tmp          = one / ( dxidx(:,1,1) * dxdxi(:,1,1) 
     &                       + dxidx(:,1,2) * dxdxi(:,2,1)  
     &                       + dxidx(:,1,3) * dxdxi(:,3,1) )
        dxidx(:,1,1) = dxidx(:,1,1) * tmp
        dxidx(:,1,2) = dxidx(:,1,2) * tmp
        dxidx(:,1,3) = dxidx(:,1,3) * tmp
        dxidx(:,2,1) = (dxdxi(:,2,3) * dxdxi(:,3,1) 
     &                - dxdxi(:,2,1) * dxdxi(:,3,3)) * tmp
        dxidx(:,2,2) = (dxdxi(:,1,1) * dxdxi(:,3,3) 
     &                - dxdxi(:,3,1) * dxdxi(:,1,3)) * tmp
        dxidx(:,2,3) = (dxdxi(:,2,1) * dxdxi(:,1,3) 
     &                - dxdxi(:,1,1) * dxdxi(:,2,3)) * tmp
        dxidx(:,3,1) = (dxdxi(:,2,1) * dxdxi(:,3,2) 
     &                - dxdxi(:,2,2) * dxdxi(:,3,1)) * tmp
        dxidx(:,3,2) = (dxdxi(:,3,1) * dxdxi(:,1,2) 
     &                - dxdxi(:,1,1) * dxdxi(:,3,2)) * tmp
        dxidx(:,3,3) = (dxdxi(:,1,1) * dxdxi(:,2,2) 
     &                - dxdxi(:,1,2) * dxdxi(:,2,1)) * tmp
c
        WdetJ = Qwt(lcsyst,intp)/ tmp

c
c.... --------------------->  Global Gradients  <-----------------------
c
        gradient = zero
c
c
        do n = 1, nshl
c
c.... compute the global gradient of shape-function
c
c            ! N_{a,x_i}= N_{a,xi_i} xi_{i,x_j}
c
          shg(:,n,1) = shgl(:,1,n) * dxidx(:,1,1) + 
     &                 shgl(:,2,n) * dxidx(:,2,1) +
     &                 shgl(:,3,n) * dxidx(:,3,1)
          shg(:,n,2) = shgl(:,1,n) * dxidx(:,1,2) + 
     &                 shgl(:,2,n) * dxidx(:,2,2) +
     &                 shgl(:,3,n) * dxidx(:,3,2) 
          shg(:,n,3) = shgl(:,1,n) * dxidx(:,1,3) + 
     &                 shgl(:,2,n) * dxidx(:,2,3) +
     &                 shgl(:,3,n) * dxidx(:,3,3) 
c
c
c  Y_{,x_i}=SUM_{a=1}^nenl (N_{a,x_i}(int) Ya)
c
          do i = 1, 3
            do j = 1, vsize
               k = (i-1)*vsize+j
               gradient(:,k) = gradient(:,k) + shg(:,n,i)*vector(:,n,j)
            end do
          end do

       end do

c
c.... return
c
       return
       end
