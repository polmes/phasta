
c-----------------------------------------------------------------------
!   Menter SST(revised 2003) Turbulence model constants
c
c-----------------------------------------------------------------------
      module turbKW
      real*8 CmuKW, kappa, a1, CDES1, CDES2,  Cd1, Cd2, tenpowerminusten
      real*8 alp1, alp2, beta1, beta2, sigk1, sigk2, sigw1, sigw2
      real*8 tenpowerminustwenty
      parameter (
     &  CmuKW                   = 0.0900d0, !!This variable is \beta^* in the NASA document
     &  kappa                   = 0.4100d0,
     &  a1                      = 0.3100d0,
     &  CDES1                   = 0.7800d0,
     &  CDES2                   = 0.6100d0,
     &  Cd1                     = 20.000d0,
     &  Cd2                     = 3.0000d0,

     &  alp1                    = 0.5555555555d0,
     &  beta1                   = 0.075000d0,
     &  sigk1                   = 0.850000d0,
     &  sigw1                   = 0.500000d0,
     &  alp2                    = 0.440000d0,
     &  beta2                   = 0.082800d0,
     &  sigk2                   = 1.00000d0,
     &  sigw2                   = 0.85600d0,
     &  gam1                    = 0.555555555d0,
     &  gam2                    = 0.44000000d0,
     &  tenpowerminusten     = 1.0d-10,
     &  tenpowerminustwenty     = 1.0d-20
     &     )

      end module
