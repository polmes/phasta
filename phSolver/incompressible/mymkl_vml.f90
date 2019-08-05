! file: mkl_vml.fi
!===============================================================================
! Copyright 2006-2018 Intel Corporation.
!
! This software and the related documents are Intel copyrighted  materials,  and
! your use of  them is  governed by the  express license  under which  they were
! provided to you (License).  Unless the License provides otherwise, you may not
! use, modify, copy, publish, distribute,  disclose or transmit this software or
! the related documents without Intel's prior written permission.
!
! This software and the related documents  are provided as  is,  with no express
! or implied  warranties,  other  than those  that are  expressly stated  in the
! License.
!===============================================================================

!++
!  Fortran 90 VML interface.
!--

!++
!  PARAMETER DEFINITIONS
!  Parameter definitions for VML mode and VML error status.
!
!  VML mode controls VML function accuracy, floating-point settings (rounding
!  mode and precision) and VML error handling options. Default VML mode is
!  VML_HA | VML_ERRMODE_DEFAULT, i.e. VML high accuracy functions are
!  called, and current floating-point precision and the rounding mode is used.
!
!  Error status macros are used for error classification.
!
!  NOTE: A | B means bitwise OR operation with A and B
!--

!  VML FUNCTION ACCURACY CONTROL
!  VML_HA - when VML_HA is set, high accuracy VML functions are called
!  VML_LA - when VML_LA is set, low accuracy VML functions are called
!  VML_EP - when VML_EP is set, enhanced performance VML functions are called
!
!  NOTE: VML_HA, VML_LA and VML_EP must not be used in combination
!      INTEGER(KIND=4) VML_LA
!      INTEGER(KIND=4) VML_HA
!      INTEGER(KIND=4) VML_EP
!      PARAMETER (VML_LA = Z"00000001")
!      PARAMETER (VML_HA = Z"00000002")
!      PARAMETER (VML_EP = Z"00000003")
!
!  SETTING OPTIMAL FLOATING-POINT PRECISION AND ROUNDING MODE
!  Definitions below are to set optimal floating-point control word
!  (precision and rounding mode).
!
!  For their correct work, VML functions change floating-point precision and
!  rounding mode (if necessary). Since control word changing is typically
!  expensive operation, it is recommended to set precision and rounding mode
!  to optimal values before VML function calls.
!
!  VML_FLOAT_CONSISTENT  - use this value if the calls are typically to single
!                          precision VML functions
!  VML_DOUBLE_CONSISTENT - use this value if the calls are typically to double
!                          precision VML functions
!  VML_RESTORE           - restore original floating-point precision and
!                          rounding mode
!  VML_DEFAULT_PRECISION - use default (current) floating-point precision and
!                          rounding mode
!  NOTE: VML_FLOAT_CONSISTENT, VML_DOUBLE_CONSISTENT, VML_RESTORE and
!        VML_DEFAULT_PRECISION must not be used in combination
!      INTEGER(KIND=4) VML_DEFAULT_PRECISION
!      INTEGER(KIND=4) VML_FLOAT_CONSISTENT
!      INTEGER(KIND=4) VML_DOUBLE_CONSISTENT
!      INTEGER(KIND=4) VML_RESTORE
!      PARAMETER (VML_DEFAULT_PRECISION = Z"00000000")
!      PARAMETER (VML_FLOAT_CONSISTENT  = Z"00000010")
!      PARAMETER (VML_DOUBLE_CONSISTENT = Z"00000020")
!      PARAMETER (VML_RESTORE           = Z"00000030")
!
!  VML ERROR HANDLING CONTROL
!  Macros below are used to control VML error handler.
!
!  VML_ERRMODE_IGNORE   - ignore errors
!  VML_ERRMODE_ERRNO    - errno variable is set on error
!  VML_ERRMODE_STDERR   - error description text is written to stderr on error
!  VML_ERRMODE_EXCEPT   - exception is raised on error
!  VML_ERRMODE_CALLBACK - user's error handler function is called on error
!  VML_ERRMODE_NOERR    - ignore errors and do not update status
!  VML_ERRMODE_DEFAULT  - errno variable is set, exceptions are raised and
!                         user's error handler is called on error
!  NOTE: VML_ERRMODE_IGNORE must not be used in combination with
!        VML_ERRMODE_ERRNO, VML_ERRMODE_STDERR, VML_ERRMODE_EXCEPT,
!        VML_ERRMODE_CALLBACK and VML_ERRMODE_DEFAULT.
!  NOTE: VML_ERRMODE_NOERR must not be used in combination with any
!        other VML_ERRMODE setting.
!      INTEGER(KIND=4) VML_ERRMODE_IGNORE
!      INTEGER(KIND=4) VML_ERRMODE_ERRNO
!      INTEGER(KIND=4) VML_ERRMODE_STDERR
!      INTEGER(KIND=4) VML_ERRMODE_EXCEPT
!      INTEGER(KIND=4) VML_ERRMODE_CALLBACK
!      INTEGER(KIND=4) VML_ERRMODE_NOERR
!      INTEGER(KIND=4) VML_ERRMODE_DEFAULT
!      PARAMETER (VML_ERRMODE_IGNORE   = Z"00000100")
!      PARAMETER (VML_ERRMODE_ERRNO    = Z"00000200")
!      PARAMETER (VML_ERRMODE_STDERR   = Z"00000400")
!      PARAMETER (VML_ERRMODE_EXCEPT   = Z"00000800")
!      PARAMETER (VML_ERRMODE_CALLBACK = Z"00001000")
!      PARAMETER (VML_ERRMODE_NOERR    = Z"00002000")
!      PARAMETER (VML_ERRMODE_DEFAULT  = IOR(VML_ERRMODE_ERRNO, IOR(VML_ERRMODE_CALLBACK,VML_ERRMODE_EXCEPT)))

!  ACCURACY, FLOATING-POINT CONTROL AND ERROR HANDLING MASKS
!  Accuracy, floating-point and error handling control are packed in
!  the VML mode variable. Macros below are useful to extract accuracy and/or
!  floating-point control and/or error handling control settings.
!
!  VML_ACCURACY_MASK           - extract accuracy bits
!  VML_FPUMODE_MASK            - extract floating-point control bits
!  VML_ERRMODE_MASK            - extract error handling control bits
!                                (including error callback bits)
!  VML_ERRMODE_STDHANDLER_MASK - extract error handling control bits
!                                (not including error callback bits)
!  VML_ERRMODE_CALLBACK_MASK   - extract error callback bits
!      INTEGER(KIND=4) VML_ACCURACY_MASK
!      INTEGER(KIND=4) VML_FPUMODE_MASK
!      INTEGER(KIND=4) VML_ERRMODE_MASK
!      INTEGER(KIND=4) VML_ERRMODE_STDHANDLER_MASK
!      INTEGER(KIND=4) VML_ERRMODE_CALLBACK_MASK
!      PARAMETER (VML_ACCURACY_MASK = Z"0000000f")
!      PARAMETER (VML_FPUMODE_MASK  = Z"000000f0")
!      PARAMETER (VML_ERRMODE_MASK  = Z"0000ff00")
!      PARAMETER (VML_ERRMODE_STDHANDLER_MASK = Z"00000f00")
!      PARAMETER (VML_ERRMODE_CALLBACK_MASK = Z"0000f000")
!
!  ERROR STATUS PARAMETER DEFINITIONS
!  VML_STATUS_OK        - no errors
!  VML_STATUS_BADSIZE   - array dimension is not positive
!  VML_STATUS_BADMEM    - invalid pointer passed
!  VML_STATUS_ERRDOM    - at least one of arguments is out of function domain
!  VML_STATUS_SING      - at least one of arguments caused singularity
!  VML_STATUS_OVERFLOW  - at least one of arguments caused overflow
!  VML_STATUS_UNDERFLOW - at least one of arguments caused underflow
!  VML_STATUS_ACCURACYWARNING - function doesn't support set accuracy mode,
!                               lower accuracy mode was used instead
!      INTEGER(KIND=4) VML_STATUS_OK
!      INTEGER(KIND=4) VML_STATUS_BADSIZE
!      INTEGER(KIND=4) VML_STATUS_BADMEM
!      INTEGER(KIND=4) VML_STATUS_ERRDOM
!      INTEGER(KIND=4) VML_STATUS_SING
!      INTEGER(KIND=4) VML_STATUS_OVERFLOW
!      INTEGER(KIND=4) VML_STATUS_UNDERFLOW
!      INTEGER(KIND=4) VML_STATUS_ACCURACYWARNING
!      PARAMETER (VML_STATUS_OK        = 0)
!      PARAMETER (VML_STATUS_BADSIZE   = -1)
!      PARAMETER (VML_STATUS_BADMEM    = -2)
!      PARAMETER (VML_STATUS_ERRDOM    = 1)
!      PARAMETER (VML_STATUS_SING      = 2)
!      PARAMETER (VML_STATUS_OVERFLOW  = 3)
!      PARAMETER (VML_STATUS_UNDERFLOW = 4)
!      PARAMETER (VML_STATUS_ACCURACYWARNING = 1000)
!
!++
!  TYPE DEFINITIONS
!--

!  ERROR CALLBACK CONTEXT.
!  Error callback context structure is used in a user's error callback
!  function with the following interface:
!
!  Error callback context fields:
!  ICODE        - error status
!  IINDEX       - index of bad argument
!  DBA1         - 1-st argument value, at which error occured
!  DBA2         - 2-nd argument value, at which error occured
!                 (2-argument functions only)
!  DBR1         - 1-st resulting value
!  DBR2         - 2-nd resulting value (2-result functions only)
!  CFUNCNAME    - function name, for which error occured
!  IFUNCNAMELEN - length of function name
!      !dec$ options /warn=noalignment
!      TYPE ERROR_STRUCTURE
!            SEQUENCE
!            INTEGER(KIND=4) ICODE
!            INTEGER(KIND=4) IINDEX
!            REAL(KIND=8)    DBA1
!            REAL(KIND=8)    DBA2
!            REAL(KIND=8)    DBR1
!            REAL(KIND=8)    DBR2
!            CHARACTER(64)   CFUNCNAME
!            INTEGER(KIND=4) IFUNCNAMELEN
!            REAL(KIND=8)    DBA1IM
!            REAL(KIND=8)    DBA2IM
!            REAL(KIND=8)    DBR1IM
!            REAL(KIND=8)    DBR2IM
!      END TYPE ERROR_STRUCTURE
!      !dec$ end options
!
!     INTERFACE
!       SUBROUTINE vdmul(n,a,b,r)
!         INTEGER,INTENT(IN) :: n
!         REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
!         REAL(KIND=8),INTENT(OUT)   :: r(n)
!       END SUBROUTINE
!     END INTERFACE
