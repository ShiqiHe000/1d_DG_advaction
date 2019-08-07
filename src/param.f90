
MODULE PARAM

    USE MPI

    IMPLICIT NONE
    
    ! INITIAL MESH FILE-------------------------------------------------
    CHARACTER(LEN=*), PARAMETER :: MESHFILE = "two_boundary.msh"    
    !-------------------------------------------------------------------
    
    ! SET POLYNOMIAL ORDER----------------------------------------------
    INTEGER :: N = 6    ! POLYNOMIAL DEGREE IN X DIRECTION
    INTEGER :: M = 4    ! POLYNOMIAL DEGREE IN Y DIRECTION
    INTEGER :: MNMAX = 12   ! MAXIMUM POLYNOMIAL DEGREE
    !-------------------------------------------------------------------
    
    ! TIME STEP---------------------------------------------------------
!    DOUBLE PRECISION :: DELTA_T = 1.0D0
    !-------------------------------------------------------------------
    
    ! WAVE SPEED -------------------------------------------------------
    DOUBLE PRECISION, PARAMETER :: C=1.0D0
    !-------------------------------------------------------------------
    
    ! TIME--------------------------------------------------------------
    DOUBLE PRECISION :: T_TOTAL = 0.5D0     ! TOTAL TIME INTEGRAL
    INTEGER :: NT = 0.3333E4                    ! TIME STEP NUMBER
    !-------------------------------------------------------------------
    
    ! SET ADAPATION ----------------------------------------------------
    INTEGER, PARAMETER :: SPLIT_MAX_NUM = 3     ! MAXIMUM SPLIT NUMBER
    !-------------------------------------------------------------------
    
    ! SET MPI ----------------------------------------------------------
    INTEGER :: RANK     ! PROCESSOR RANK
    INTEGER :: NUM_PROC ! NUMBER OF PROCESSOR
    INTEGER :: IERROR 
    !-------------------------------------------------------------------
    
    ! OUTPUT FILE ------------------------------------------------------
    INTEGER :: FRAME=1  
    !-------------------------------------------------------------------
    
    ! OTHER PARAMETER --------------------------------------------------
    DOUBLE PRECISION, PARAMETER :: PI=4.0D0*DATAN(1.0D0)
    !-------------------------------------------------------------------

END MODULE PARAM
