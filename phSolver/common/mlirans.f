      module mlirans

      real*8, dimension(:),   allocatable :: dns_x
      real*8, dimension(:),   allocatable :: dns_y
      real*8, dimension(:,:), allocatable :: dns_data
      real*8, dimension(:,:), allocatable :: dns_flow

      integer :: dns_x_length
      integer :: dns_y_length
      integer :: dns_total_length
      integer :: dns_num_fields
      integer :: disc_flag
      
      end module

      module neural_network

      real*8, dimension(:,:), allocatable :: w1,w2,w3,w4,w5,w6
      real*8, dimension(:),   allocatable :: b1,b2,b3,b4,b5,b6
      real*8, dimension(:),   allocatable :: y_hat
      real*8, dimension(:,:), allocatable :: a1,a2,a3,a4,a5,a6,feat

      real*8, parameter :: y1min =   -0.3258,     y1max = 1.8953
      real*8, parameter :: y2min =   -1.6015,     y2max = 0.5858
      real*8, parameter :: y3min =   -0.4918,     y3max = 0.4950
            
      integer :: nn = 250
      integer :: input_number = 10
      integer :: output_number = 3

      end module

      module polynomial_expansion

      real*8,  dimension(:),  allocatable :: y_hat_pe,feat_pe
      real*8,  dimension(:),  allocatable :: coeff1,coeff2,coeff3,psi
      integer, dimension(:,:),allocatable :: inds1,inds2,inds3

      integer :: n_basis = 21
      integer :: input_number_pe = 10
      integer :: output_number_pe = 3

      end module
 
      subroutine  load_dns()
      use mlirans
      include "common.h"
      include "mpif.h"

      logical exlog

c.... set discrepancy flag equal to 0 for first time step
      disc_flag = 0

c.... check to see if x_interp.dat file exists        
      inquire(file="input_files/x_interp.dat",exist=exlog)
      if(exlog) then

c.... write to screen if file exist
      if(myrank == 0)write(*,*)"x_interp.dat exists"

c.... read vector of x coordinates
         OPEN(unit=556,file = "input_files/x_interp.dat")
         read(556,*)dns_x_length
         allocate(dns_x(dns_x_length))
         do i = 1,dns_x_length
            read(556,*)dns_x(i)
         enddo
         close(556)

c.... end conditional if x_interp.dat file exists
      endif

c.... check to see if y_interp.dat file exists        
      inquire(file="input_files/y_interp.dat",exist=exlog)
      if(exlog) then

c.... write to screen if file exist
      if(myrank == 0)write(*,*)"y_interp.dat exists"

c.... read matrix of y coordinates 
c.... where y(i*dns_x_length+j)
c.... i = y height and j = x length
         OPEN(unit=556,file = "input_files/y_interp.dat")
         read(556,*)dns_total_length
         dns_y_length = dns_total_length/dns_x_length
         allocate(dns_y(dns_total_length))
         do i = 1,dns_total_length
            read(556,*)dns_y(i)
         enddo
         close(556)

c.... end conditional if y_interp.dat file exists
      endif

c.... check to see if dns_interp.dat file exists        
      inquire(file="input_files/xy_data.dat",exist=exlog)
      if(exlog) then

c.... write to screen if file exist
         if(myrank == 0)write(*,*)"xy_data.dat exists"

c.... typical dns field data includes
c.... 1.)  u
c.... 2.)  v
c.... 3.)  p
c.... 4.)  s11
c.... 5.)  s22
c.... 6.)  s12
c.... 7.)  omega12
c.... 8.)  eps
c.... 9.)  tke
c.... 10.) nut
c.... 11.) p,x
c.... 12.) p,y
c.... 13.) tke,x
c.... 14.) tke,y
c.... 15.) DU/Dt
c.... 16.) delta11
c.... 17.) delta22
c.... 18.) delta12

c.... read matrix of y coordinates 
c.... where y(i*dns_x_length+j)
c.... i = y height and j = x length
         OPEN(unit=556,file = "input_files/xy_data.dat")
         read(556,*)dns_total_length,dns_num_fields
         allocate(dns_data(dns_total_length,dns_num_fields))
         do i = 1,dns_total_length
            read(556,*)(dns_data(i,j),j=1,dns_num_fields)
         enddo
         close(556)

c.... end conditional if y_interp.dat file exists
      endif
 
      return 

      end subroutine

      subroutine set_dns_bc(point2x, BC, ndsurf)

      use mlirans
      include "common.h"
      include "mpif.h"

      real*8  :: point2x(numnp,nsd)
      real*8  :: BC(nshg,ndofBC)
      integer :: ndsurf(nshg)
      real*8  :: xi,eta
      integer :: isrfidmatch = 555,iBCon,i,k,j,l,m

c.... ordering of BC's is p,T,u,v,w for 1-5
      if(myrank == 0)write(*,*)"setting inlet BC's"

c.... loop through all nodal points on processor
      do i=1,nshg                    
              
c.... check to see if surf ID matches inlet ID
         iBcon = 0
         if((ndsurf(i).eq.isrfidmatch)) iBCon=1

c.... if ID matches interpolate BC
         if(iBCon == 1)then

c.... set x coordinate of inflow
            xi = point2x(i,1)

c.... find floor of x coordinate
            l=1
            do m = 1,dns_x_length-1
              if(xi > dns_x(m))then
                 l = m
              endif
            enddo

c.... set y coordinate of inflow
            eta = point2x(i,2)

c.... find floor of y coordinate
            k=1
            do j = 1,dns_y_length
               if(eta > dns_y(j*dns_x_length+0))then
                  k = j
               endif
            enddo
            if(k > dns_y_length-2)k = dns_y_length-2 

c.... set BC's
c.... set x velocity on face with 555 surf ID
            BC(i,3) = dns_data(k*dns_x_length+l,1) +
     &             (eta - dns_y(k*dns_x_length+l))*
     & ((dns_data((k+1)*dns_x_length+l,1)-dns_data(k*dns_x_length+l,1))/
     &  (dns_y((k+1)*dns_x_length+l)-dns_y(k*dns_x_length+l)))

c.... set y velocity on face with 555 surf ID
            BC(i,4) = dns_data(k*dns_x_length+l,2) +
     &             (eta - dns_y(k*dns_x_length+l))*
     & ((dns_data((k+1)*dns_x_length+l,2)-dns_data(k*dns_x_length+l,2))/
     &  (dns_y((k+1)*dns_x_length+l)-dns_y(k*dns_x_length+l)))


      if (SARS_flag == 0 )then    
c.... set eddy viscosity on face with 555 surf ID
            BC(i,7) = dns_data(k*dns_x_length+l,10) +
     &             (eta - dns_y(k*dns_x_length+l))*
     & ((dns_data((k+1)*dns_x_length+l,10)-dns_data(k*dns_x_length+l,10))/
     &  (dns_y((k+1)*dns_x_length+l)-dns_y(k*dns_x_length+l)))
      endif


c.... end condional if surf ID matches
         endif

c.... loop over nodal points on processor
      enddo

      return

      end subroutine

      subroutine set_dns_ic(x,y)

      use mlirans
      include "common.h"
      include "mpif.h"

      real*8,intent(inout) :: y(nshg,ndof),x(numnp,nsd)  
      integer :: i,j,k,l,m
      real*8 :: xi,eta,sol
      real*8 :: a0,a1,a2,a3,denom,denom2
      real*8 :: f11,f12,f21,f22,x1,y1,x2,y2

c.... this particular bilinear interpolation will only work for tensor product grids
c.... interpolation coordinate input files must be in increasing order

      if(myrank==0)write(*,*)"interpolating dns data"

c.... get dns flow features at nodal coordinates at the same time as setting IC's
      allocate(dns_flow(nshg,dns_num_fields))

      do l = 1,nshg
 
         xi  = x(l,1)
         eta = x(l,2)

c.... find floor of x coordinate
         j=1
         do i = 1,dns_x_length-1
            if(xi > dns_x(i))then
               j = i
            endif
         enddo

c.... find floor of y coordinate
         k=1
         do i = 1,dns_y_length-1
            if(eta > dns_y(i*dns_x_length+j))then
               k = i
            endif
         enddo

c.... compute interpolation constants for given point
         x1 = dns_x(j)
         x2 = dns_x(j+1)
         y1 = dns_y(k*dns_x_length+j)
         y2 = dns_y((k+1)*dns_x_length+j)
         denom = (x1-x2)*(y1-y2)
         denom2 = (x1-x2)*(y2-y1)
 
c.... loop over fields
         do m = 1,dns_num_fields

            f11 = dns_data((k-1)*dns_x_length+j,m)
            f12 = dns_data((k-1+1)*dns_x_length+j,m)
            f21 = dns_data((k-1)*dns_x_length+j+1,m)
            f22 = dns_data((k-1+1)*dns_x_length+j+1,m)

c.... compute coefficients for bi-linear interpolation
            a0 = f11*x2*y2/denom+f12*x2*y1/denom2+f21*x1*y2/denom2+f22*x1*y1/denom
            a1 = f11*y2/denom2+f12*y1/denom+f21*y2/denom+f22*y1/denom2
            a2 = f11*x2/denom2+f12*x2/denom+f21*x1/denom+f22*x1/denom2
            a3 = f11/denom+f12/denom2+f21/denom2+f22/denom

c.... compute solution from interpolation
            dns_flow(l,m) = a0 + a1*xi + a2*eta + a3*xi*eta

c.... if outside interpolation range use closet values
            if(j == dns_x_length-1)dns_flow(l,m) = dns_data((k-1)*dns_x_length+j,m)
            if(k == dns_y_length-1)dns_flow(l,m) = dns_data((k-1+1)*dns_x_length+j,m)
            if(j == dns_x_length-1 .and. k == dns_y_length-1) dns_flow(l,m) = dns_data((k-1+1)*dns_x_length+j+1,m)

c.... end loop over fields
         enddo

c.... end loop over points
      enddo

c.... y values = u,v,w,p,T,scalar for 1-6
      if(lstep == 0 )then
         if(myrank==0)write(*,*)"setting IC's in domain"
            do i = 1,nshg
               y(i,1) = dns_flow(i,1)
               y(i,2) = dns_flow(i,2)
               y(i,4) = dns_flow(i,3)
            enddo
      endif

c.... deallocate dns data that was used to get nodal data
      deallocate(dns_x)
      deallocate(dns_y)
      deallocate(dns_data) 

      return 
 
      end subroutine

      subroutine SDC(dns_flow_l, Dij, nut, ql, gradYl, disc_flag, xl)
    
      use neural_network
      use polynomial_expansion 
      include "common.h"

      real*8, intent(inout) :: dns_flow_l(npro,nshl,18)
      real*8, intent(inout) :: Dij(npro,nshl,6)
      real*8, intent(inout) :: nut(npro,nshl)
      real*8, intent(inout) :: ql(npro,nshl,9)
      real*8, intent(inout) :: gradYl(npro,nshl,12)
      real*8, intent(inout) :: xl(npro,nshl,nsd)
      integer, intent(in)   :: disc_flag
      real*8  :: nuu,sI2,slam1
      real*8  :: lam1,lam2,lam3,alpha,beta,gamma,tke,eps
      real*8  :: in1,in2,in3,in4,in5
      real*8,dimension(3,3) :: R1,R2,R3,R4,Dterm,Tone,Ttwo,Tthree,Tfour
      real*8,dimension(3,3) :: Sij,Oij,GPij,Gkij
      real*8 :: lam1s,lam2s,lam3s,alphas,betas,gammas
      real*8 :: lam1d,lam2d,lam3d,alphad,betad,gammad
      real*8 :: err1,err2,err3,sup1,sup2,sup3
      real*8 :: S_denom,O_denom,gradp_denom,gradk_denom
      real*8 :: ex1_denom,ex2_denom,energy
      integer :: i,j,k

c.... set kinematic viscosity in routine for clarity
      nuu = datmat(1,2,1)

c.... zero out error counters
      err1 = 0.00d0
      err2 = 0.00d0
      err3 = 0.00d0
      sup1 = 0.00d0
      sup2 = 0.00d0
      sup3 = 0.00d0
 
c.... get element nodal values of discrepancy terms
c.... loop over elements in block
      do i = 1,npro

c.... loop over number of shape functions on element
         do j = 1,nshl

c.... compute eddy viscosity from dns k and epsilon
            nut(i,j) = dns_flow_l(i,j,10)

c.... if discrepancy flag = 0 then set Dij (Dij will remain constant after first step)
            if(disc_flag == 0)then

c.... compute current mean strain rate tensor from phasta variables
c                Sij(:,:) = 0.00d0
c                Sij(1,1) = dns_flow_l(i,j,4)
c                Sij(2,2) = dns_flow_l(i,j,5)
c                Sij(1,2) = dns_flow_l(i,j,6)

c                Oij(:,:) = 0.00d0
c                Oij(1,2) = dns_flow_l(i,j,7)

                GPij(:,:) = 0.00d0
                GPij(2,3) = dns_flow_l(i,j,11)
                GPij(3,1) = dns_flow_l(i,j,12)

                Sij(:,:) = 0.00d0
                Sij(1,1) = 1.00d0/2.00d0*(gradYl(i,j,1)+gradYl(i,j,1))
                Sij(2,2) = 1.00d0/2.00d0*(gradYl(i,j,5)+gradYl(i,j,5))
                Sij(1,2) = 1.00d0/2.00d0*(gradYl(i,j,2)+gradYl(i,j,4))

                Oij(:,:) = 0.00d0
                Oij(1,2) = 1.00d0/2.00d0*(gradYl(i,j,2)-gradYl(i,j,4))

c                GPij(:,:) = 0.00d0
c                GPij(2,3) = gradYl(i,j,10)
c                GPij(3,1) = gradyl(i,j,11)

                Gkij(:,:) = 0.00d0
                Gkij(2,3) = dns_flow_l(i,j,13)
                Gkij(3,1) = dns_flow_l(i,j,14)

c.... compute principle components of the mean strain rate tensor           
!               call principle(Sij(1,1),         Sij(2,2)         ,0.00d0,
!     &                        Sij(1,2),         0.00d0,           0.00d0,
!     &                        lam1s,            lam2s,            lam3s,
!     &                        alphas,           betas,            gammas)

c.... re-align Eigen-vectors for S such that they are most closely aligned to grad P
               call order_tensor(alphas,      betas,        gammas,
     &                           GPij(2,3),   GPij(3,1))

c.... compute tensors from dns or mean flow properties
c.... compute denominators for non-dimensionalization
               S_denom = dns_flow_l(i,j,8)/dns_flow_l(i,j,9)
               O_denom = sqrt(2.00d0*dns_flow_l(i,j,7)**2)
               gradp_denom = dns_flow_l(i,j,15)
               gradk_denom = dns_flow_l(i,j,8)/sqrt(dns_flow_l(i,j,9))
               ex1_denom = 10000.00*nuu*sqrt(Sij(1,1)**2+2.00d0*Sij(1,2)**2+Sij(2,2)**2)
               ex2_denom = 1.00d0/sqrt(Sij(1,1)**2+2.00d0*Sij(1,2)**2+Sij(2,2)**2)

c.... S tensor
               Tone(1,1) = Sij(1,1)
               Tone(1,2) = Sij(1,2)
               Tone(1,3) = 0.00d0
               Tone(2,1) = Sij(1,2)
               Tone(2,2) = Sij(2,2)
               Tone(2,3) = 0.00d0
               Tone(3,1) = 0.00d0
               Tone(3,2) = 0.00d0
               Tone(3,3) = 0.00d0

c.... non-dimensionlize S
               Tone = Tone/(abs(Tone)+abs(S_denom))

c.... Omega tensor
               Ttwo(1,1) = 0.00d0
               Ttwo(1,2) = Oij(1,2)
               Ttwo(1,3) = 0.00d0
               Ttwo(2,1) = -Oij(1,2)
               Ttwo(2,2) = 0.00d0
               Ttwo(2,3) = 0.00d0
               Ttwo(3,1) = 0.00d0
               Ttwo(3,2) = 0.00d0
               Ttwo(3,3) = 0.00d0

c.... non-dimensionalize Omega
               Ttwo = Ttwo/(abs(Ttwo)+abs(O_denom))

c.... gradient of P tensor eps_ijk p_,k
               Tthree(1,1) = 0.00d0
               Tthree(1,2) = 0.00d0
               Tthree(1,3) = -GPij(3,1)
               Tthree(2,1) = 0.00d0
               Tthree(2,2) = 0.00d0
               Tthree(2,3) = GPij(2,3)
               Tthree(3,1) = GPij(3,1)
               Tthree(3,2) = -GPij(2,3)
               Tthree(3,3) = 0.00d0

c.... non-dimensionalize grad p tensor
               Tthree = Tthree/(abs(Tthree)+abs(gradp_denom))

c.... gradient of k tensor eps_ijk k_,k
               Tfour(1,1) = 0.00d0
               Tfour(1,2) = 0.00d0
               Tfour(1,3) = -Gkij(3,1)
               Tfour(2,1) = 0.00d0
               Tfour(2,2) = 0.00d0
               Tfour(2,3) = Gkij(2,3)
               Tfour(3,1) = Gkij(3,1)
               Tfour(3,2) = -Gkij(2,3)
               Tfour(3,3) = 0.00d0

c.... non-dimensionalize grad k tensor
               Tfour = Tfour/(abs(Tfour) + abs(gradk_denom))

c.... set input features for NN
               call double_contraction(feat(1,1),Tone,Tone)
               !call double_contraction(feat(1,2),Ttwo,Ttwo)
               call double_contraction(feat(1,2),Tthree,Tthree)
               call double_contraction(feat(1,3),Tfour,Tfour)
               call triple_contraction(feat(1,4),Tthree,Tthree,Tone)
               call triple_contraction(feat(1,5),Tfour,Tfour,Tone)
               call double_contraction(feat(1,6),Tthree,Tfour)
               !call quadruple_contraction(feat(1,8),Tthree,Tthree,Ttwo,Tone)
               !call quadruple_contraction(feat(1,9),Tfour,Tfour,Ttwo,Tone)
               call triple_contraction(feat(1,7),Tthree,Tfour,Tone)
               call triple_contraction(feat(1,8),Ttwo,Tthree,Tfour)
               !call quadruple_contraction(feat(1,12),Ttwo,Tthree,Tfour,Tone)
               feat(1,9) = dns_flow_l(i,j,9)/
     &                     (abs(dns_flow_l(i,j,9))+abs(ex1_denom))
c               feat(1,14) = (dns_flow_l(i,j,9)/dns_flow_l(i,j,8))/
c     &                      (abs(dns_flow_l(i,j,9)/dns_flow_l(i,j,8))+abs(ex2_denom))                   
               feat(1,10) = nut(i,j)/
     &                      (abs(nut(i,j))+10.00d0*abs(nuu))

c.... check to see if Nans 
               do k = 1,10
                  if(isNAN(feat(1,k)))then
                     feat(1,k) = 0.00d0
                  endif
               enddo

c.... predict SARS variables with Neural network
               call nn_predict()

c.... re-normalize outputs
               call normalize(y_hat(1),y1min,y1max,0.00d0,1.00d0)
               call normalize(y_hat(2),y2min,y2max,0.00d0,1.00d0)
               call normalize(y_hat(3),y3min,y3max,0.00d0,1.00d0)
      
c.... re-dimensionalize
               y_hat(1) = y_hat(1)*2.00d0*abs(dns_flow_l(i,j,9))
               y_hat(2) = y_hat(2)*2.00d0*abs(dns_flow_l(i,j,9))
               y_hat(3) = y_hat(3)*2.00d0*abs(dns_flow_l(i,j,9))

c.... check for energy stability
c               energy = -2.0000*(nuu+nut(i,j))*lams1
c               if(y_hat(1) >= energy)then
c                  y_hat(1) = energy
c                  y_hat(3) = 0.00d0
c               endif
c               energy = -2.0000*(nuu+nut(i,j))*lams2
c               if(y_hat(2) >= energy)then
c                  y_hat(2) = energy
c                  y_hat(3) = 0.00d0
c               endif
c               y_hat(:) = 0.000d0

c.... compute discrepancy from SARS
               Dterm(:,:) = 0.00d0
               Dterm(1,1) = y_hat(1)!dns_flow_l(i,j,16)
               Dterm(2,2) = y_hat(2)!dns_flow_l(i,j,17)
               Dterm(1,2) = y_hat(3)!dns_flow_l(i,j,18)
               Dterm(2,1) = y_hat(3)!dns_flow_l(i,j,18)

c.... if statement for bad dns data on separation bubble
               if(xl(i,j,1) < -4.00d0)then
                  Dterm(:,:) = 0.00d0
                  Dterm(1,1) = dns_flow_l(i,j,16)
                  Dterm(2,2) = dns_flow_l(i,j,17)
                  Dterm(1,2) = dns_flow_l(i,j,18)
                  Dterm(2,1) = dns_flow_l(i,j,18)
               endif

c.... compute error 
            err1 = err1 + ((dns_flow_l(i,j,16)-y_hat(1)))**2 
            err2 = err2 + ((dns_flow_l(i,j,17)-y_hat(2)))**2
            err3 = err3 + ((dns_flow_l(i,j,18)-y_hat(3)))**2
            sup1 = sup1 + dns_flow_l(i,j,16)**2
            sup2 = sup2 + dns_flow_l(i,j,17)**2
            sup3 = sup3 + dns_flow_l(i,j,18)**2      

c.... compute rotation from Euclidean space to Eigen-vector space of S
               call rotation(R1,0.00d0,0.00d0,gammas)

c.... rotate discrepancy tensor in principle S frame back to Euclidean
c.... D_e = R D_s R^T
               Dterm = matmul(R1,Dterm)
               R1 = transpose(R1)
               Dterm = matmul(Dterm,R1)

c.... load D into global arrays
               dns_flow_l(i,j,1) = Dterm(1,1)
               dns_flow_l(i,j,2) = Dterm(2,2)
               dns_flow_l(i,j,3) = Dterm(1,2)

c.... end conditional for computing discrepancy term on first time step
            endif

c.... compute discrepancy tensor from rotation
             Dij(i,j,1) = dns_flow_l(i,j,1) + 2.00d0/3.00d0*dns_flow_l(i,j,9)
             Dij(i,j,2) = dns_flow_l(i,j,2) + 2.00d0/3.00d0*dns_flow_l(i,j,9)
             Dij(i,j,3) = 0.00d0
             Dij(i,j,4) = dns_flow_l(i,j,3) 
             Dij(i,j,5) = 0.00d0
             Dij(i,j,6) = 0.00d0

c.... end loop over number of shape functions on element
         enddo

c.... end loop over elements on block
      enddo

c.... print error to screen
c      write(*,*)"training error: d11",sqrt(err1)/sqrt(sup1),"d22",sqrt(err2)/sqrt(sup2),"d12",sqrt(err3)/sqrt(sup3)

      return 

      end subroutine

      subroutine rotation(R,alpha,beta,gamma)
      real*8, intent(in) :: alpha, beta, gamma
      real*8, intent(inout), dimension(3,3) :: R
      real*8, dimension(3,3) :: R1,R2,R3

c.... rotation about x-axis
      R1(1,1) = 1.00d0
      R1(1,2) = 0.00d0
      R1(1,3) = 0.00d00
      R1(2,1) = 0.00d0
      R1(2,2) = cos(alpha)
      R1(2,3) = -sin(alpha)
      R1(3,1) = 0.00d0
      R1(3,2) = sin(alpha)
      R1(3,3) = cos(alpha)

c.... rotation about y-axis
      R2(1,1) = cos(beta)
      R2(1,2) = 0.00d0
      R2(1,3) = sin(beta)
      R2(2,1) = 0.00d0
      R2(2,2) = 1.00d0
      R2(2,3) = 0.00d0
      R2(3,1) = -sin(beta)
      R2(3,2) = 0.00d0
      R2(3,3) = cos(beta)

c.... rotation about z-axis
      R3(1,1) = cos(gamma)
      R3(1,2) = -sin(gamma)
      R3(1,3) = 0.00d0
      R3(2,1) = sin(gamma)
      R3(2,2) = cos(gamma)
      R3(2,3) = 0.00d0
      R3(3,1) = 0.00d0
      R3(3,2) = 0.00d0
      R3(3,3) = 1.00d0

c.... total euler angle rotation matrix
      R = matmul(R1,R2)
      R = matmul(R,R3)

      return

      end subroutine

      subroutine euler_angles(alpha,beta,gamma,R)
      real*8, intent(in), dimension(3,3) :: R
      real*8, intent(inout) :: alpha,beta,gamma
      real*8 :: temp,pi = 3.141592653589793238462643383279502884
      real*8 :: theta1,theta2,psi1,psi2,phi1,phi2

c.... given eigen vectors in rotation matrix find euler angles
      temp = abs(R(3,1))-1.0000d0
      if(temp < 1e-12)then
         theta1 = -asin(R(3,1))
         theta2 = pi-theta1
         psi1 = atan2(R(3,2)/cos(theta1),R(3,3)/cos(theta1))
         psi2 = atan2(R(3,2)/cos(theta2),R(3,3)/cos(theta2))
         phi1 = atan2(R(2,1)/cos(theta1),R(1,1)/cos(theta1))
         phi2 = atan2(R(2,1)/cos(theta2),R(1,1)/cos(theta2))
      else
         phi1 = 0.0000
         if(R(3,1) <= -.9999)then
            theta1 = pi/2
            psi1 = phi1+atan2(R(1,2),R(1,3))
         else
            theta1 = -pi/2
            psi1 = -phi1+atan2(-R(1,2),-R(1,3))
         endif
      endif

      alpha = psi1
      beta = theta1
      gamma = phi1

      return

      end subroutine

      subroutine cross_product(a,b,c)
      real*8, intent(inout), dimension(3) :: a,b
      real*8, intent(inout), dimension(3) :: c

c.... compute components of vector resulting for cross product
      c(1) = a(2) * b(3) - a(3) * b(2)
      c(2) = a(3) * b(1) - a(1) * b(3)
      c(3) = a(1) * b(2) - a(2) * b(1)

      return

      end subroutine

      subroutine load_network()
   
      use neural_network

      allocate(feat(1,input_number))
      allocate(y_hat(output_number))

      allocate(a1(1,nn))
      allocate(a2(1,nn))
      allocate(a3(1,nn))
      allocate(a4(1,nn))
      allocate(a5(1,nn))
      allocate(a6(1,output_number))

      allocate(w1(input_number,nn))
      call parse_w(w1,input_number,nn,'input_files/network/w1.dat')

      allocate(w2(nn,nn))
      call parse_w(w2,nn,nn,'input_files/network/w2.dat')

      allocate(w3(nn,nn))
      call parse_w(w3,nn,nn,'input_files/network/w3.dat')

      allocate(w4(nn,nn))
      call parse_w(w4,nn,nn,'input_files/network/w4.dat')

      allocate(w5(nn,nn))
      call parse_w(w5,nn,nn,'input_files/network/w5.dat')

      allocate(w6(nn,output_number))
      call parse_w(w6,nn,output_number,'input_files/network/w6.dat')

      allocate(b1(nn))
      call parse_b(b1,nn,'input_files/network/b1.dat')

      allocate(b2(nn))
      call parse_b(b2,nn,'input_files/network/b2.dat')

      allocate(b3(nn))
      call parse_b(b3,nn,'input_files/network/b3.dat')

      allocate(b4(nn))
      call parse_b(b4,nn,'input_files/network/b4.dat')

      allocate(b5(nn))
      call parse_b(b5,nn,'input_files/network/b5.dat')

      allocate(b6(output_number))
      call parse_b(b6,output_number,'input_files/network/b6.dat')

      return

      end subroutine

      subroutine parse_w(weights,m,n,file_name)
      integer :: m,n
      real*8 :: weights(m,n)
      character(len=26), intent(in) :: file_name

      open(unit=10,file = file_name)
      do i = 1,m
         do j = 1,n
            read(10,*)weights(i,j)
         enddo
      enddo
      close(10)

      return

      end subroutine

      subroutine parse_b(bias,m,file_name)
      integer :: m
      real*8 :: bias(m)
      character(len=26), intent(in) :: file_name

      open(unit=10,file = file_name)
      do i = 1,m
         read(10,*)bias(i)
      enddo
      close(10)

      return

      end subroutine

      subroutine relu(x,m)
      integer, intent(in) :: m
      real *8, intent(inout) :: x(m)

      do i = 1,m
         if(x(i) < 0.0000000000000000)x(i) = 0.0000000000000000
      enddo

      return

      end subroutine

      subroutine tanhh(x,m)
      integer, intent(in) :: m
      real *8, intent(inout) :: x(m)
     
      do i = 1,m
         x(i) = (exp(2.00d0*x(i))-1.00d0)/(exp(2.00d0*x(i))+1.00d0) 
      enddo
  
      return
 
      end subroutine

      subroutine sigmoid(x,m)
      integer, intent(in) :: m
      real *8, intent(inout) :: x(m)

      do i = 1,m
         x(i) = 1.0000000000000000/(1.0000000000000000+exp(-x(i)))
      enddo

      return

      end subroutine

      subroutine normalize(val,min_new,max_new,min_old,max_old)
      real*8, intent(in) :: min_old,max_old,min_new,max_new
      real*8, intent(inout) :: val

      val = (max_new - min_new)/(max_old-min_old)*(val-max_old)+max_new

      return

      endsubroutine

      subroutine nn_predict()
      use neural_network

      a1 = matmul(feat,w1)
      a1(1,:) = a1(1,:) + b1
      call relu(a1(1,:),nn)
      a2 = matmul(a1,w2)
      a2(1,:) = a2(1,:) + b2
      call relu(a2(1,:),nn)
      a3 = matmul(a2,w3)
      a3(1,:) = a3(1,:) + b3
      call relu(a3(1,:),nn)
      a4 = matmul(a3,w4)
      a4(1,:) = a4(1,:) + b4
      call sigmoid(a4(1,:),nn)
      a5 = matmul(a4,w5)
      a5(1,:) = a5(1,:) + b5
      call sigmoid(a5(1,:),nn)
      a6 = matmul(a5,w6)
      a6(1,:) = a6(1,:) + b6
      call sigmoid(a6(1,:),output_number)
      y_hat(:) = a6(1,:)

      return
      end subroutine

      subroutine load_expansion()
      use polynomial_expansion

      allocate(feat_pe(input_number_pe))
      allocate(y_hat_pe(output_number_pe))
      allocate(psi(n_basis))

      allocate(coeff1(n_basis))
      call parse_coeff(coeff1,n_basis,'input_files/expansion/coeff1.dat')

      allocate(inds1(n_basis,input_number_pe))
      call parse_inds(inds1,n_basis,input_number_pe,
     & 'input_files/expansion/inds1.dat')

      allocate(coeff2(n_basis))
      call parse_coeff(coeff2,n_basis,'input_files/expansion/coeff2.dat')

      allocate(inds2(n_basis,input_number_pe))
      call parse_inds(inds2,n_basis,input_number_pe,
     & 'input_files/expansion/inds2.dat')

      allocate(coeff3(n_basis))
      call parse_coeff(coeff3,n_basis,'input_files/expansion/coeff3.dat')

      allocate(inds3(n_basis,input_number_pe))
      call parse_inds(inds3,n_basis,input_number_pe,
     & 'input_files/expansion/inds3.dat')

      return

      endsubroutine

      subroutine parse_coeff(coeff,m,file_name)
      integer :: m
      real*8 :: coeff(m)
      character(len=32), intent(in) :: file_name

      open(unit=10,file = file_name)
      do i = 1,m
         read(10,*)coeff(i)
      enddo
      close(10)

      return

      end subroutine

      subroutine parse_inds(inds,m,n,file_name)
      integer :: m,n
      integer :: inds(m,n)
      character(len=31), intent(in) :: file_name

      open(unit=10,file = file_name)
      do i = 1,m
         read(10,*)(inds(i,j),j=1,n)
      enddo
      close(10)

      return

      end subroutine


      subroutine pe_predict(input_name)
      use polynomial_expansion
      character(len=6), intent(in) :: input_name

      if(input_name == 'angle3')then
         call build_psi(feat_pe,psi,inds1,n_basis,input_number_pe)
         y_hat_pe(1) = dot_product(psi,coeff1)
      elseif(input_name == 'delta1')then
         call build_psi(feat_pe,psi,inds2,n_basis,input_number_pe)
         y_hat_pe(2) = dot_product(psi,coeff2)
      elseif(input_name == 'delta2')then
         call build_psi(feat_pe,psi,inds3,n_basis,input_number_pe)
         y_hat_pe(3) = dot_product(psi,coeff3)
      else
         write(*,*)"invalid input into prediction"
      endif

      return

      end subroutine

      subroutine build_psi(x,psi,inds,n_basis,input_number)
      integer, intent(in) :: n_basis,input_number
      integer, intent(in), dimension(n_basis,input_number) :: inds
      real*8, intent(in), dimension(input_number) :: x
      real*8, intent(inout),dimension(n_basis) :: psi
      real*8, allocatable, dimension(:) :: Lvec,LL
      integer :: i,j

      allocate(Lvec(input_number))

      do i = 1, n_basis
         do j = 1,input_number
            allocate(LL(inds(i,j)+1))
            call legendre_poly(x(j),inds(i,j),LL)
            Lvec(j) = LL(inds(i,j)+1)
            deallocate(LL)
         enddo
         psi(i) = product(Lvec)
      enddo

      deallocate(Lvec)

      return

      end subroutine

      subroutine legendre_poly(val,p,LL)
      real*8, intent(in) :: val
      integer, intent(in) :: p
      real*8, intent(inout),dimension(p+1) :: LL
      integer :: n,k,tmp
      real*8 :: magL

      LL(:) = 0.00d0

      do n = 0,p
         do k = 0,n
            call nchoosek(n,k,tmp)
            LL(n+1) = LL(n+1) + tmp**2.00*(val-1.00d0)**(n-k)*(val+1.00d0)**k
         enddo
         LL(n+1) = LL(n+1)/(2.00d0**n)
         magL = sqrt(2.00d0/(2.00d0*n+1.00d0))
         LL(n+1) = LL(n+1)/magL
      enddo

      return

      end subroutine

      subroutine nchoosek(n,k,val)
      integer,intent(in) :: n,k
      integer,intent(out) :: val
      integer :: a,c,f,i

      a = product((/(i, i=1,n)/))
      c = product((/(i, i=1,k)/))
      f = product((/(i, i=1,n-k)/))
      val = a/(c*f)

      return

      end subroutine

      subroutine order_tensor(alpha,beta,gamma,a,b)
      real(kind = 8), intent(inout) :: alpha,beta,gamma,a,b
      real(kind = 8) :: mag,tmp,tmp1,tmp2
      real(kind = 8), dimension(3,3) :: R
      real(kind = 8), dimension(3) :: v1,v2,v3,nv,gp,cp,v1n,v2n

      !compute rotation matrix from Euler angles
      call rotation(R,alpha,beta,gamma)

      !compute normalized vector to align with Eigen-vector
      mag = a**2 + b**2
      a = a/sqrt(mag)
      b = b/sqrt(mag)

      !load components in to normalized vector
      gp(1) = a
      gp(2) = b
      gp(3) = 0.00d0

      !choose Eigen vector which is closest to normalized vector
      tmp = dot_product(gp,R(:,1))
      v1n = -R(:,1)
      if(dot_product(gp,v1n)<tmp)R(:,1) = v1n

      tmp = dot_product(gp,R(:,2))
      v2n = -R(:,2)

      if(dot_product(gp,v2n)<tmp)R(:,2) = v2n

      call euler_angles(alpha,beta,gamma,R)

      !need to change back because type inout
      a = a*sqrt(mag)
      b = b*sqrt(mag)

      return

      end subroutine

      subroutine double_contraction(invariant,one,two)
      real(kind=8),intent(in),dimension(3,3) :: one,two
      real(kind=8),intent(inout) :: invariant
      integer :: i,j,dim
      dim = 3

      invariant = 0.00d0
      do i = 1,dim
         do j = 1,dim
            invariant = invariant + one(i,j)*two(j,i)
         enddo
      enddo

      return

      endsubroutine

      subroutine triple_contraction(invariant,one,two,three)
      real(kind=8),intent(in),dimension(3,3) :: one,two,three
      real(kind=8),intent(inout) :: invariant
      integer :: i,j,k,dim
      dim = 3

      invariant = 0.00d0
      do i = 1,dim
         do j = 1,dim
            do k = 1,dim
               invariant = invariant + one(i,j)*two(j,k)*three(k,i)
            enddo
         enddo
      enddo

      return

      endsubroutine

      subroutine quadruple_contraction(invariant,one,two,three,four)
      real(kind=8),intent(in),dimension(3,3) :: one,two,three,four
      real(kind=8),intent(inout) :: invariant
      integer :: i,j,k,l,dim
      dim = 3

      invariant = 0.00d0
      do i = 1,dim
         do j = 1,dim
            do k = 1,dim
               do l = 1,dim
                  invariant = invariant + one(i,j)*two(j,k)*three(k,l)*four(l,i)
               enddo
            enddo
         enddo
      enddo

      return

      endsubroutine

      subroutine quintuple_contraction(invariant,one,two,three,four,five)
      real(kind=8),intent(in),dimension(3,3) :: one,two,three,four,five
      real(kind=8),intent(inout) :: invariant
      integer :: i,j,k,l,m,dim
      dim = 3

      invariant = 0.00d0
      do i = 1,dim
         do j = 1,dim
            do k = 1,dim
               do l = 1,dim
                  do m = 1,dim
                     invariant = invariant + one(i,j)*two(j,k)*three(k,l)*four(l,m)*five(m,i)
                  enddo
               enddo 
            enddo
         enddo
      enddo

      return

      endsubroutine

      subroutine sextuple_contraction(invariant,one,two,three,four,five,six)
      real(kind=8),intent(in),dimension(3,3) :: one,two,three,four,five,six
      real(kind=8),intent(inout) :: invariant
      integer :: i,j,k,l,m,n,dim
      dim = 3

      invariant = 0.00d0
      do i = 1,dim
         do j = 1,dim
            do k = 1,dim
               do l = 1,dim
                  do m = 1,dim
                     do n = 1,dim
                        invariant = invariant + one(i,j)*two(j,k)*three(k,l)*four(l,m)*five(m,n)*six(m,i)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      return

      endsubroutine

      subroutine set_wmles_ic(x,y)

      include "common.h"
      include "mpif.h"

      real*8,intent(inout) :: y(nshg,ndof),x(numnp,nsd)
      integer ic_length, num_fields, ind, l, ii
      real*8 xi, yi, r, rtmp
      real*8,  dimension(:,:),  allocatable :: ic_data
      logical exlog
 
      inquire(file="IC_interp.dat",exist=exlog)
      if(exlog) then
        if(myrank == 0)write(*,*)"IC_interp.dat exists"
        OPEN(unit=556,file = "IC_interp.dat")
        read(556,*)ic_length,num_fields
        allocate(ic_data(ic_length,num_fields))
        do i = 1,ic_length
            read(556,*)(ic_data(i,j),j=1,num_fields)
        enddo
        close(556)
      else
        if(myrank == 0) write(*,*) "IC_interp.dat not found"
      endif
 
      if (exlog) then
      if(myrank.eq.0)write(*,*)"Setting IC to 2D field in IC_interp.dat"
      do l = 1,nshg

c.... do a closest point search
         ! coordinates of mesh node
         xi  = x(l,1)
         yi = x(l,2)
         if (xi.lt.zero) then
         ! loop through points in imported data to find closest
         r = 100.0
         do ii=1,ic_length
            rtmp = sqrt((xi-ic_data(ii,1))**2+(yi-ic_data(ii,2))**2)
            if (rtmp<r) then
               r = rtmp
               ind = ii
            endif
            if (r.lt.1.0d-10) exit
         enddo
         ! assign the flow vars at the closest point to y
         y(l,1) = ic_data(ind,4) ! u
         y(l,2) = ic_data(ind,5) ! v
         y(l,4) = ic_data(ind,3) ! p
         endif

      enddo
      deallocate(ic_data)
      endif
     
      return

      endsubroutine



