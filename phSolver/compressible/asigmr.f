        subroutine AsIGMR (blk,y,       ac,  xl,dwl,        xmudmi,   
     &                     shp,     shgl,    ien,     
     &                     mater,   res,     rmes,    
     &                     BDiag,   qres,    EGmass,   rerr)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  interior elements.
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        use rlssave     ! Use the resolved Leonard stresses at the nodes.
        use timedataC    ! time series
        use specialBC    ! get ytarget to localize and send down
        use eblock
        include "common.h"
      type (LocalBlkData) blk
c
        dimension y(nshg,ndofl),            ac(nshg,ndofl),
     &            dwl(blk%e,blk%n),         shp(nshl,MAXQPT),  
     &            shgl(nsd,nshl,MAXQPT),
     &            ien(npro,nshl),  
     &            mater(npro),              res(nshg,nflow),
     &            rmes(nshg,nflow),         BDiag(nshg,nflow,nflow),
     &            qres(nshg,idflx)

c
        dimension ycl(npro,nshl,ndofl),     acl(npro,nshl,ndof),
     &            xl(npro,nenl,nsd),        ytargetl(npro,nshl,nflow),
     &            rl(npro,nshl,nflow),      rml(npro,nshl,nflow),
     &            BDiagl(npro,nshl,nflow,nflow),
     &            ql(npro,nshl,idflx)
c        
        dimension  xmudmi(npro,ngauss)
        dimension sgn(npro,nshl),  EGmass(npro,nedof,nedof)

        dimension rlsl(npro,nshl,6) 
        real*8 rerrl(npro,nshl,6), rerr(nshg,numerr)

c
c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
c
        if (ipord .gt. 1) then
           call getsgn(blk, ien,sgn)
        endif
c
c.... gather the variables
c
        call localy(blk,y,      ycl,     ien,    ndofl,  'gather  ')
        call localy(blk,ac,    acl,     ien,    ndofl,  'gather  ')
        call local (blk,qres,   ql,     ien,    idflx,  'gather  ')

        if(matflg(5,1).ge.4 ) 
     &   call localy (blk,ytarget,   ytargetl,  
     &                ien,   nflow,  'gather  ')


        if( (iLES.gt.10).and.(iLES.lt.20)) then  ! bardina 
           call local (blk,rls, rlsl,     
     &     ien,       6, 'gather  ')  
        else
           rlsl = zero
        endif      
c
c.... get the element residuals, LHS matrix, and preconditioner
c
        rl     = zero
        BDiagl = zero

        if(ierrcalc.eq.1) rerrl = zero
        ttim(31) = ttim(31) - secs(0.0)

        call e3  (blk,ycl,     ycl,     acl,     shp,
     &            shgl,    xl,      rl,      rml,   xmudmi,
     &            BDiagl,  ql,      sgn,     rlsl,  EGmass,
     &            rerrl,   ytargetl, dwl)

        ttim(31) = ttim(31) + secs(0.0)
c
c.... assemble the residual and modified residual
c
        call local (blk,res,    rl,     ien,    nflow,  'scatter ')
c
        if ( ierrcalc .eq. 1 ) then
           call local (blk,rerr, rerrl,  ien, 6, 'scatter ')
        endif
c
c.... extract and assemble the Block-Diagonal (see note in elmgmr, line 280)
c
        if (iprec .ne. 0) then 
           do i = 1, nshl
              do j = 1, nflow
                 i0 = (i - 1) * nflow + j
                 do k = 1, nflow
                    j0 = (i - 1) * nflow + k
                    BDiagl(:,i,j,k) = EGmass(:,i0,j0)
                 enddo
              enddo
           enddo
           call local (blk,BDiag,  BDiagl, ien, nflow*nflow, 'scatter ')
        endif
        
c
c... call timeseries
c

        if (exts) then
           if ((iter.eq.1).and.(mod(lstep,freq).eq.0)) then
              call timeseries(blk,ycl,xl,ien,sgn)
           endif
        endif
        
c
c.... end
c
        return
        end
c
c
c
        subroutine AsIGMRSclr (blk, y,       ac,
     &                         xl,       dwl,
     &                         shp,     shgl,      ien,     
     &                         mater,   rest,      rmest,    
     &                         qrest,   EGmasst,   Diag)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  interior elements.
c
c Zdenek Johan, Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        use eblock
        use turbSA
        include "common.h"
      type (LocalBlkData) blk
c
        dimension y(nshg,ndof),             
     &            ac(nshg,ndof),
     &            shp(nshl,MAXQPT),        shgl(nsd,nshl,MAXQPT),
     &            ien(npro,nshl),
     &            mater(npro),            rest(nshg),
     &            rmest(nshg),            Diag(nshg),
     &            qrest(nshg)

c
        dimension ycl(npro,nshl,ndof),       
     &            acl(npro,nshl,ndof),    dwl(npro,nenl),
     &            xl(npro,nenl,nsd),      Diagl(npro,nshl),
     &            rtl(npro,nshl),         rmtl(npro,nshl),
     &            qtl(npro,nshl),         sgn(npro,nshl)
c        
        dimension EGmasst(npro,nshape, nshape)
c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
        call getsgn(blk, ien,sgn)
c
c
c.... gather the variables
c
        call localy (blk, y,       ycl,      ien,    ndof,  'gather  ')
        call localy (blk, ac,      acl,     ien,    ndof,  'gather  ')
c       call local (qrest,   qtl,     ien,    1,     'gather  ')
c
c.... get the element residuals, LHS matrix, and preconditioner
c
        rtl     = zero
        Diagl   = zero
        
        ttim(31) = ttim(31) - tmr()

        call e3Sclr (blk, ycl,     acl,  
     &               dwl,     dwl,   shp,
     &               sgn,    shgl,    xl,
     &               rtl,    rmtl,
     &               qtl,    EGmasst )

        ttim(31) = ttim(31) + tmr()
c
c.... assemble the residual and modified residual
c
        call local (blk, rest,    rtl,     ien,    1,  'scatter ')
c
c.... extract and assemble the Diagonal
c
        if (iprec .ne. 0) then 
           do i=1,nshl
              Diagl(:,i)=EGmassT(:,i,i)
           enddo
           call local(blk,Diag, Diagl, ien, 1, 'scatter ')
        endif
c
c.... end
c
        return
        end


