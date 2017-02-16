        subroutine e3dcSclr ( blk, gradS,    giju,     gGradS,
     &                        rLS,      tauS,     srcR,
     &                        dcFct)
c
c
c----------------------------------------------------------------------
c
c This routine calculates the contribution of the Discontinuity-
c Capturing operator to RHS and preconditioner for the scalar solve.
c
c  g1yti   (nflow,blk%e)           : grad-y in direction 1
c  g2yti   (nflow,blk%e)           : grad-y in direction 2
c  g3yti   (nflow,blk%e)           : grad-y in direction 3
c  A0     (nsymdf,blk%e)          : A0 matrix (Symm. storage)
c  raLS   (blk%e)                 : square of LS residual (A0inv norm)
c  rtLS   (blk%e)                 : square of LS residual (Tau norm)
c  giju    (6,blk%e)              : metric matrix
c  DC     (ngauss,blk%e)          : discontinuity-capturing factor
c  intp				 : integration point number
c
c output:
c  ri     (nflow*(nsd+1),blk%e)   : partial residual
c  rmi    (nflow*(nsd+1),blk%e)   : partial modified residual
c  stiff  (nsymdf,6,blk%e)       : diffusivity matrix
c  DC     (blk%e)                : discontinuity-capturing factor
c
c
c Zdenek Johan, Summer 1990. (Modified from e2dc.f)
c Zdenek Johan, Winter 1991. (Recoded)
c Zdenek Johan, Winter 1991. (Fortran 90)
c----------------------------------------------------------------------
c
      use eblock
      include "common.h"
      type (LocalBlkData) blk

c
        dimension gradS(blk%e,nsd),            gGradS(blk%e,nsd),
     &            rLS(blk%e),                  tauS(blk%e),
     &            giju(blk%e,6),               dcFct(blk%e),
     &            srcR(blk%e)
c
c.... Form GijUp gradS and  gradS . GijUp gradS (store in dcFct)
c
	
	    gGradS(:,1) = GijU(:,1) * gradS(:,1)
     1			+ GijU(:,4) * gradS(:,2)
     2			+ GijU(:,6) * gradS(:,3)
	    gGradS(:,2) = GijU(:,4) * gradS(:,1)
     1			+ GijU(:,2) * gradS(:,2)
     2			+ GijU(:,5) * gradS(:,3)
	    gGradS(:,3) = GijU(:,6) * gradS(:,1)
     1			+ GijU(:,5) * gradS(:,2)
     2			+ GijU(:,3) * gradS(:,3)
c
	    dcFct(:)    = gradS(:,1) * gGradS(:,1)
     1		        + gradS(:,2) * gGradS(:,2)
     2		        + gradS(:,3) * gGradS(:,3)
     3		        + epsM
	
	    dcFct(:) = 1.0/ dcFct(:)
c
c.... Form pdeRes 2-norm / gradT 2-norm
c

	    dcFct  = dcFct * (rLS - srcR) ** 2 
c
c.... ------------------------->  DC factor  <------------------------
c
c.... DC-mallet
c
	    if (idcsclr(1) .eq. 1) then
c       
	       fact = one
	       if (ipord .eq. 2)  fact = 0.9
	       if (ipord .eq. 3) fact = 0.75
c       
c$$$  dcFct(:)=dim((fact*sqrt(dcFct(:))),(tauS(:)*dcFct(:))) !not work
                                                          !with all compilers
	       dcFct(:)=max(zero,(fact*sqrt(dcFct(:)))-(tauS(:)*dcFct(:)))
c
	    endif
c       
c       
c....   DC-quadratic
c       
	    if (idcsclr(1) .eq. 2) then
c       
	       dcFct(:) = two * tauS(:) * dcFct(:)
c       
	    endif
c       
c....   DC-min
c       
	    if (idcsclr(1) .eq. 3) then
c       
	       fact = one
	       if (ipord .eq. 2)  fact = 0.9
c       
          dcFct(:) = min( max(zero, (fact * sqrt(dcFct(:)) -
     &	             tauS(:)*dcFct(:)) ), two * tauS(:) * dcFct(:))
c       
	    endif
c
c.... Scale the gGradT for residual formation
c	
	    gGradS(:,1) = dcFct(:) * gGradS(:,1)
	    gGradS(:,2) = dcFct(:) * gGradS(:,2)
	    gGradS(:,3) = dcFct(:) * gGradS(:,3)
	


	return
	end
c
