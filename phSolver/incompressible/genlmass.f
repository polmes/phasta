        subroutine genlmass (x, shp,shgl, iBC, iper, ilwork)
c
        use pointer_data
c
        include "common.h"
      include "eblock.h"
      type (LocalBlkData) blk

        include "mpif.h"
c
        real*8 x(numnp,nsd)
c
        real*8 shp(MAXTOP,maxsh,MAXQPT),   
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT) 
        integer iBC(nshg), iper(nshg), ilwork(nlwork)
c
        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        
c
c gmass came in via pointer_data and will 
c be available wherever it is included  (allocate it now).
c

        allocate (gmass(nshg))  
        gmass=zero
c
c.... loop over the element-blocks
c
        do iblk = 1, nelblk
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk) ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          npro   = lcblk(1,iblk+1) - iel 
          inum   = iel + npro - 1
          ngauss = nint(lcsyst)
          blk%n   = lcblk(5,iblk) ! no. of vertices per element
          blk%s   = lcblk(10,iblk)
          blk%e   = lcblk(1,iblk+1) - iel
          blk%g = nint(lcsyst)
          blk%l = lcblk(3,iblk)
          blk%o = lcblk(4,iblk)


c
c
c.... compute and assemble the residual and tangent matrix
c
          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))
          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)
         

          call AsImass (blk,x,       tmpshp,     
     &                  tmpshgl, mien(iblk)%p,
     &                  gmass)

          deallocate ( tmpshp )
          deallocate ( tmpshgl )
c
c.... end of interior element loop
c
       enddo
c
c.... -------------------->   communications <-------------------------
c
         if (numpe > 1) then
            call commu (gmass, ilwork, 1, 'in ')
         end if
c
c  take care of periodic boundary conditions on mass matrix
c
       do j= 1,nshg
         if ((btest(iBC(j),10))) then
           i = iper(j)
           gmass(i) = gmass(i) + gmass(j)
         endif
       enddo

       do j= 1,nshg
         if ((btest(iBC(j),10))) then
           i = iper(j)
           gmass(j) = gmass(i)
         endif
       enddo
c
      return
      end

        subroutine AsImass (blk,x,      shp,
     &                     shgl,    ien,     
     &                     gmass)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the mass corresponding to the
c each node.
c
c Ken Jansen, Winter 2000.  (Fortran 90)
c----------------------------------------------------------------------
c
      include "common.h"
      include "eblock.h"
      type (LocalBlkData) blk
c
        real*8 x(numnp,nsd),              
     &         shp(nshl,maxsh),       shgl(nsd,nshl,ngauss),
     &         gmass(nshg)

        integer ien(npro,nshl)

c
        real*8    xl(bsz,nenl,nsd),    WdetJ(npro), 
     &            sgn(npro,nshl),       shape(npro,nshl),          
     &            locmass(bsz,nshl),   shg(npro,nshl,nsd),
     &            fmstot(npro),         temp(npro),
     &            dxidx(npro,nsd,nsd),  shdrv(npro,nsd,nshl)

        integer aa
c        
c
c
c.... gather the variables
c
c
c.... get the matrix of mode signs for the hierarchic basis functions. 
c
        if (ipord .gt. 1) then
           write(*,*) 'blk not plumbed this far'

           call getsgn(blk,ien,sgn)
        endif
        
        call localx(blk,x,      xl,     ien,    nsd,    'gather  ')
c
c.... zero the matrices if they are being recalculated
c

        locmass=zero
        fmstot=zero

        do intp = 1, ngauss

           if (Qwt(lcsyst,intp) .eq. zero) cycle ! precaution
c
c.... get the hierarchic shape functions at this int point
c
           call getshp(blk,intp,shp,          shgl,      sgn, 
     &                 shape,        shdrv)

c
c.... --------------------->  Element Metrics  <-----------------------
c
           call e3metric(blk,intp, xl,         shdrv,       dxidx,  
     &                    shg,        WdetJ)

c
c  get this quad points contribution to the integral of the square of  the 
c  shape function
c
           do aa = 1,nshl
              locmass(1:blk%e,aa)= locmass(1:blk%e,aa) 
     &             + shape(:,aa)*shape(:,aa)*WdetJ
           enddo
c
c also accumulate this quad points contribution to the integral of the element
c volume (integral Na^2 d Omega)
c 
           fmstot= fmstot + WdetJ ! intregral  d Omega
c
c.... end of integration loop
c
        enddo
c
c.... lumped mass if needed   Note that the locmass factors accumulated
c     over integration points and weighted with WdetJ already.
c

c.... scale the LHS matrix contribution with special lumping weighting
c
c  The first term we collect is the trace of integral Na^2 d Omega
c
        temp = zero
        do aa = 1, nshl
           temp = temp + locmass(1:blk%e,aa) !reusing temp to save memory
        enddo

c
c scale the diagonal so that the trace will still yield Omega^e (the volume
c of the element)
c
        do aa = 1, nshl
           locmass(1:blk%e,aa) = locmass(1:blk%e,aa) * fmstot / temp
        enddo
c
c.... assemble the residual
c
        call local (blk,gmass,    locmass,     ien,    1,  'scatter ')

c
c.... end
c
        return
        end

      subroutine lmassadd ( ac,       res,
     &                      rowp,     colm,    
     &                      lhs16,     gmass)
c     
      include "common.h"
c     
      real*8 ac(nshg,ndof), res(nshg,4), tmp,tmp1
      real*8 lhs16(16,nnz_tot), gmass(nshg), rho(nshg)
      integer rowp(nnz*nshg),  colm(nshg+1)
      integer	n,	k
c
      integer sparseloc
c
c
      rho=datmat(1,1,1)  ! needs to be generalized for VOF or level set
      tmp1=flmpl*almi
      if((flmpl.ne.0).and.(lhs.eq.1)) then
c
c.... Add lmass to diag of lhs16
c
         do n = 1, nshg
	    k = sparseloc( rowp(colm(n)), colm(n+1)-colm(n), n )
     &       + colm(n)-1
            tmp=gmass(n)*tmp1*rho(n)
	    lhs16(1,k) = lhs16(1,k) + tmp
	    lhs16(6,k) = lhs16(6,k) + tmp
	    lhs16(11,k) = lhs16(11,k) + tmp
         enddo
      endif

      tmp1=flmpr

      if(flmpr.ne.0) then
         rho=rho*gmass*tmp1  ! reuse rho
         res(:,1)=res(:,1)-ac(:,1)*rho(:)
         res(:,2)=res(:,2)-ac(:,2)*rho(:)
         res(:,3)=res(:,3)-ac(:,3)*rho(:)
      endif
     
      return
      end

      subroutine lmassaddSclr ( ac,       res,
     &                          rowp,     colm,    
     &                          lhsS,     gmass)
c     
      include "common.h"
c     
      real*8 ac(nshg),       res(nshg), tmp, tmp1
      real*8 lhsS(nnz_tot), gmass(nshg), rho(nshg)
      integer rowp(nnz*nshg),  colm(nshg+1)
      integer	n,	k
c
      integer sparseloc
c
c
      rho=datmat(1,1,1)  ! needs to be generalized for VOF or level set
      tmp1=flmpl*almi
      if((flmpl.ne.0).and.(lhs.eq.1)) then
c
c.... Add lmass to diag of lhsK
c
         do n = 1, nshg
	    k = sparseloc( rowp(colm(n)), colm(n+1)-colm(n), n )
     &       + colm(n)-1
            tmp=gmass(n)*tmp1*rho(n)
	    lhsS(k) = lhsS(k) + tmp
         enddo
      endif

      tmp1=flmpr
      if(flmpr.ne.0) then
         rho=rho*gmass*tmp1  ! reuse rho
         res(:)=res(:)-ac(:)*rho(:)
      endif
     
      return
      end
