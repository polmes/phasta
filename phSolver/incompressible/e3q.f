        subroutine e3q (blk,yl,      dwl,     shp,     shgl,
     &                  xl,      ql,      rmassl, 
     &                  xmudmi,  sgn )
c                                                                      
c----------------------------------------------------------------------
c
c This routine computes the element contribution to the 
c diffusive flux vector and the lumped mass matrix.
c
c input: 
c  yl     (bsz,blk%s,ndof)       : Y variables
c  shp    (nen,blk%g)            : element shape-functions
c  shgl   (nsd,nen,blk%g)        : element local-grad-shape-functions
c  xl     (bsz,blk%s,nsd)        : nodal coordinates at current step
c  
c output:
c  ql     (bsz,blk%s,idflx) : element RHS diffusion residual 
c  rmassl     (bsz,blk%s)        : element lumped mass matrix
c
c----------------------------------------------------------------------
c
        use spat_var_eps   ! use spatially-varying epl_ls
c
        include "common.h"
      include "eblock.h"
      type (LocalBlkData) blk

c
        dimension yl(bsz,blk%s,ndof),     dwl(bsz,blk%n),
     &            shp(blk%s,blk%g),      shgl(nsd,blk%s,blk%g),
     &            xl(bsz,blk%n,nsd),
     &            ql(bsz,blk%s,idflx),  rmassl(bsz,blk%s),
     &            xmudmi(blk%e,blk%g)
c
c local arrays
c
        dimension g1yi(blk%e,nflow),           g2yi(blk%e,nflow),
     &            g3yi(blk%e,nflow),           shg(blk%e,blk%s,nsd),
     &            dxidx(blk%e,nsd,nsd),       WdetJ(blk%e),
     &            rmu(blk%e) 
c
        dimension qdi(blk%e,idflx),alph1(blk%e),alph2(blk%e)
c
        dimension sgn(blk%e,blk%s),          shape(blk%e,blk%s),
     &            shdrv(blk%e,nsd,blk%s),    shpsum(blk%e)

        real*8 tmp(blk%e)
c
c.... for surface tension
c     
        dimension g1yti(blk%e),          g2yti(blk%e),
     &            g3yti(blk%e)
        integer idflow
c
c.... loop through the integration points
c
        
        
        alph1 = 0.d0
        alph2 = 0.d0
        
        do intp = 1, blk%g
        if (Qwt(blk%l,intp) .eq. zero) cycle          ! precaution
c     
        call getshp(blk,intp, shp,          shgl,      sgn, 
     &              shape,        shdrv)
        
c
c.... initialize
c
        qdi = zero
c
c
c.... calculate the integration variables necessary for the
c     formation of q
c

        call e3qvar   (blk,yl,        shdrv,   
     &                 xl,           g1yi,
     &                 g2yi,      g3yi,         shg,
     &                 dxidx,     WdetJ )      
c  
        idflow = 9   ! we ALWAYS save space for tau_{ij} in q_i 
                     ! even if idiff is not greater than 1

        if(idiff >= 1) then   !so taking care of all the idiff=1,3
c
c.... compute diffusive fluxes 
c
c.... compute the viscosity
c
        call getdiff(blk,intp,dwl, yl, shape, xmudmi, xl, rmu, tmp,
     &               elem_local_size(blk%i))
c
c.... diffusive flux in x1-direction
c
        qdi(:,1) =  two * rmu *  g1yi(:,2)
        qdi(:,4) =        rmu * (g1yi(:,3) + g2yi(:,2))
        qdi(:,7) =        rmu * (g1yi(:,4) + g3yi(:,2))
c     
c.... diffusive flux in x2-direction
c
        qdi(:,2) =        rmu * (g1yi(:,3) + g2yi(:,2))
        qdi(:,5) =  two * rmu *  g2yi(:,3)
        qdi(:,8) =        rmu * (g2yi(:,4) + g3yi(:,3))
c     
c.... diffusive flux in x3-direction
c
        qdi(:,3) =        rmu * (g1yi(:,4) + g3yi(:,2))
        qdi(:,6)=         rmu * (g2yi(:,4) + g3yi(:,3))
        qdi(:,9)=  two * rmu *  g3yi(:,4)
c
c
c.... assemble contribution of qdi to ql,i.e., contribution to 
c     each element node
c     
        do i=1,blk%s
           ql(1:blk%e,i,1 ) = ql(1:blk%e,i,1 )+ shape(:,i)*WdetJ*qdi(:,1 )
           ql(1:blk%e,i,2 ) = ql(1:blk%e,i,2 )+ shape(:,i)*WdetJ*qdi(:,2 )
           ql(1:blk%e,i,3 ) = ql(1:blk%e,i,3 )+ shape(:,i)*WdetJ*qdi(:,3 )

           ql(1:blk%e,i,4 ) = ql(1:blk%e,i,4 )+ shape(:,i)*WdetJ*qdi(:,4 )
           ql(1:blk%e,i,5 ) = ql(1:blk%e,i,5 )+ shape(:,i)*WdetJ*qdi(:,5 )
           ql(1:blk%e,i,6 ) = ql(1:blk%e,i,6 )+ shape(:,i)*WdetJ*qdi(:,6 )

           ql(1:blk%e,i,7 ) = ql(1:blk%e,i,7 )+ shape(:,i)*WdetJ*qdi(:,7 )
           ql(1:blk%e,i,8 ) = ql(1:blk%e,i,8 )+ shape(:,i)*WdetJ*qdi(:,8 )
           ql(1:blk%e,i,9 ) = ql(1:blk%e,i,9 )+ shape(:,i)*WdetJ*qdi(:,9 )

        enddo
c
c.... compute and assemble the element contribution to the lumped
c     mass matrix
c
c
c.... row sum technique
c
        if ( idiff == 1 ) then
           do i=1,blk%s
              rmassl(1:blk%e,i) = rmassl(1:blk%e,i) + shape(:,i)*WdetJ
           enddo
        endif
c
c.... "special lumping technique" (Hughes p. 445)
c
        if ( idiff == 3 ) then
           shpsum = zero
           do i=1,blk%s
              shpsum = shpsum + shape(:,i)*shape(:,i)
              rmassl(1:blk%e,i)=rmassl(1:blk%e,i)+shape(:,i)*shape(:,i)*WdetJ
           enddo
           alph1 = alph1+WdetJ
           alph2 = alph2+shpsum*WdetJ
        endif
      endif                     ! end of idiff=1 .or. 3 
c
c Compute gradient of level set scalar to be used to compute
c normal vector later in qpbc & e3ivar
c
      if(isurf .eq. 1) then
c
c.... initialize
c
        g1yti   = zero
        g2yti   = zero
        g3yti   = zero
c
c.... calculate the integration variables necessary for the
c     formation of q
c
c.... compute the global gradient of Yt-variables, assuming 6th entry as 
c.... the phase indicator function 
c
c  Yt_{,x_i}=SUM_{a=1}^blk%s (N_{a,x_i}(int) Yta)
c
        do n = 1, blk%s
          g1yti(:)  = g1yti(:)  + shg(:,n,1) * yl(1:blk%e,n,6)
          g2yti(:)  = g2yti(:)  + shg(:,n,2) * yl(1:blk%e,n,6)
          g3yti(:)  = g3yti(:)  + shg(:,n,3) * yl(1:blk%e,n,6)
        enddo
c
c    computing N_{b}*N_{a,x_i)*yta*WdetJ
c
        do i=1,blk%s
           ql(1:blk%e,i,idflow+1)  = ql(1:blk%e,i,idflow+1)  
     &                       + shape(:,i)*WdetJ*g1yti
           ql(1:blk%e,i,idflow+2)  = ql(1:blk%e,i,idflow+2)  
     &                       + shape(:,i)*WdetJ*g2yti
           ql(1:blk%e,i,idflow+3)  = ql(1:blk%e,i,idflow+3)  
     &                       + shape(:,i)*WdetJ*g3yti
           rmassl(1:blk%e,i) = rmassl(1:blk%e,i) + shape(:,i)*WdetJ
        enddo
      endif  !end of the isurf  
c
c.... end of the loop over integration points
c
      enddo
c
c.... normalize the mass matrix for idiff == 3
c
      if ( idiff == 3 ) then
         do i=1,blk%s
            rmassl(1:blk%e,i) = rmassl(1:blk%e,i)*alph1/alph2
         enddo
      endif
      

c
c.... return
c
       return
       end


        subroutine e3qSclr (blk,yl,      dwl,     shp,     shgl,
     &                      xl,      ql,      rmassl, 
     &                      sgn )
c                                                                      
c----------------------------------------------------------------------
c
c This routine computes the element contribution to the 
c diffusive flux vector and the lumped mass matrix.
c
c----------------------------------------------------------------------
c
        include "common.h"
      include "eblock.h"
      type (LocalBlkData) blk

c
        dimension yl(bsz,blk%s,ndof),    dwl(bsz,blk%n),
     &            shp(blk%s,blk%g),      shgl(nsd,blk%s,blk%g),
     &            xl(bsz,blk%n,nsd),
     &            ql(bsz,blk%s,nsd),     rmassl(bsz,blk%s)
c
c local arrays
c
        dimension gradT(blk%e,nsd),     
     &            dxidx(blk%e,nsd,nsd),       WdetJ(blk%e)
c
        dimension qdi(blk%e,nsd),alph1(blk%e),alph2(blk%e)
c
        dimension sgn(blk%e,blk%s),          shape(blk%e,blk%s),
     &            shdrv(blk%e,nsd,blk%s),    shpsum(blk%e)

        real*8 diffus(blk%e)
c
c.... loop through the integration points
c
        
        
        alph1 = 0.d0
        alph2 = 0.d0
        
        do intp = 1, blk%g
        if (Qwt(blk%l,intp) .eq. zero) cycle          ! precaution
c     
        call getshp(blk,intp, shp,          shgl,      sgn, 
     &              shape,        shdrv)
        
c
c.... initialize
c
        qdi = zero
c
c
c.... calculate the integration variables necessary for the
c     formation of q 
c
        call e3qvarSclr   (blk,yl,        shdrv,   
     &                     xl,        gradT,
     &                     dxidx,     WdetJ )        

c
c.... compute diffusive flux vector at this integration point
c
        call getdiffsclr(blk,shape, dwl, yl, diffus)

c
c.... diffusive flux 
c
        qdi(:,1) =  diffus * gradT(:,1)
        qdi(:,2) =  diffus * gradT(:,2)
        qdi(:,3) =  diffus * gradT(:,3)
c
c
c.... assemble contribution of qdi to ql,i.e., contribution to 
c     each element node
c     
        do i=1,blk%s
           ql(1:blk%e,i,1 ) = ql(1:blk%e,i,1 )+ shape(:,i)*WdetJ*qdi(:,1 )
           ql(1:blk%e,i,2 ) = ql(1:blk%e,i,2 )+ shape(:,i)*WdetJ*qdi(:,2 )
           ql(1:blk%e,i,3 ) = ql(1:blk%e,i,3 )+ shape(:,i)*WdetJ*qdi(:,3 )

        enddo
c
c.... compute and assemble the element contribution to the lumped
c     mass matrix
c
c
c.... row sum technique
c
        if ( idiff == 1 ) then
           do i=1,blk%s
              rmassl(1:blk%e,i) = rmassl(1:blk%e,i) + shape(:,i)*WdetJ
           enddo
        endif
c
c.... "special lumping technique" (Hughes p. 445)
c
        if ( idiff == 3 ) then
           shpsum = zero
           do i=1,blk%s
              shpsum = shpsum + shape(:,i)*shape(:,i)
              rmassl(1:blk%e,i)=rmassl(1:blk%e,i)+shape(:,i)*shape(:,i)*WdetJ
           enddo
           alph1 = alph1+WdetJ
           alph2 = alph2+shpsum*WdetJ
        endif
c
c.... end of the loop over integration points
c
      enddo
c
c.... normalize the mass matrix for idiff == 3
c
      if ( idiff == 3 ) then
         do i=1,blk%s
            rmassl(1:blk%e,i) = rmassl(1:blk%e,i)*alph1/alph2
         enddo
      endif
c
c.... return
c
       return
       end

