
MODULE DG_WAVE

USE MPI
USE BASIS
USE PARAM, ONLY: N, C, T_TOTAL, NT
USE USER_DEFINED

IMPLICIT NONE

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: DER_TRANSPOSE   ! THE TRANSPOSE OF THE FIRST DERIVATIVE MATRIX
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SOLUTION   ! SOLUTION AT GL POINTSJ
DOUBLE PRECISION :: DELTA_T

CONTAINS

SUBROUTINE DG_WAVE_PROCEDURE
! ----------------------------------------------------------------------
! SOLVE A SCALER ADVACTION PROBLEM
! INTEGRATE THE SYSTEM FROM INITIAL TO FINAL TIME 
! USE GAUSS-LENGENDRE POINTS
! MODIFIED ALORITHM 51
!-----------------------------------------------------------------------

    IMPLICIT NONE
    
    INTEGER :: K
    
    DOUBLE PRECISION :: TN = 0.0D0 ! CURRENT TIME
    DOUBLE PRECISION :: ERROR(0:N)
    DOUBLE PRECISION :: EXACT(0:N)

    
    !-------------------------------------------------------------------
    ! ALLOCATE
    ALLOCATE(SOLUTION(0:N))
    SOLUTION=0.0D0
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    ! COMPUTE GAUSS-LEGENDRE POINTS AND WEIGHTS
    ! COMPUTE THE FIRST ORDER PLOYMOIAL DERIVATIVE MATRIX
    ! COMPUTE THE LAGANDRE POLYNOMIAL AT THE TWO ENDS OF THE BOUNDARIES
    CALL CONSTRUCT_DG
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    ! USER DEFINED INITIAL CONDITION
    CALL INITIAL_CONDITION(N, SOLUTION)
    !-------------------------------------------------------------------
    
!    PRINT *, SOLUTION
    
    !-------------------------------------------------------------------
    ! CONFIRM THE TIME STEP
    CALL DECIDE_TIME_STEP(N)
!    DELTA_T = 1.5E-4
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
!    DO K=0, NT-1
    DO K=0,0

        CALL DG_STEP_BY_RK3(TN, N)
        TN = (K+1)*DELTA_T
    ENDDO
    !-------------------------------------------------------------------
    
    ! ------------------------------------------------------------------
    ! GET ERROR
    CALL GET_ERROR(N, ERROR, EXACT, TN)
    !-------------------------------------------------------------------
    
!    PRINT *, MAXVAL(DABS(ERROR))
    
    !-------------------------------------------------------------------
    ! WRITE SOLUTION TO FILE
!    OPEN(UNIT=6, FILE="advaction_solution.txt")
!    DO K=0, N
!        WRITE(6, 30) K, GL_POINT(K), SOLUTION(K), EXACT(K), ERROR(K)
!30 FORMAT(I2, 2X, F10.5, 2X, E20.10, 2X, E20.10, 2X, E20.10)        
    
!    ENDDO
!    CLOSE(UNIT=6)
    !-------------------------------------------------------------------
    
    ! WRITE DATA TO FILE------------------------------------------------
    CALL WRITE_DATA(N, EXACT, TN)   ! FILE = "advaction_solution.dat"
    !-------------------------------------------------------------------
    
    ! GET ERROR MAX-----------------------------------------------------
!    OPEN(UNIT = 6, FILE="max_error.dat", ACCESS="APPEND")
!    WRITE(6, 70) N, MAXVAL(DABS(ERROR))
!70 FORMAT(I3, 2X, E20.10)
!    CLOSE(UNIT = 6)
    !-------------------------------------------------------------------

END SUBROUTINE DG_WAVE_PROCEDURE

SUBROUTINE CONSTRUCT_DG
!-----------------------------------------------------------------------
! CONSTRUCTOR FOR THE DISCONTINUOUS GALERKIN CLASS
! ALGORITHM 59
!-----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER :: I, J
    
    !-------------------------------------------------------------------
    ALLOCATE(LAGRANGE_LEFT(0: N), LAGRANGE_RIGHT(0: N))
    ALLOCATE(DER_TRANSPOSE(0: N, 0: N))
    
    LAGRANGE_LEFT = 0.0D0
    LAGRANGE_RIGHT = 0.0D0
    DER_TRANSPOSE = 0.0D0
    !-------------------------------------------------------------------
    
    ! GL POINTS AND WEIGHTS---------------------------------------------
    CALL GL(N)
    !-------------------------------------------------------------------
    
    ! FIRST DERIVATIVE MATRIX-------------------------------------------
    CALL mth_Order_Polynomial_Derivative_Matrix(N,1, GL_POINT)    ! USE GL POINTS, 
    !--------------------------------------------------------------------
    
    ! LAGRANGE POLY AT POINT -1.0 AND +1.0 -----------------------------
    CALL LAGRANGE_INTERPOLATING_POLYNOMIAL(N, -1.0D0, GL_POINT, BARY, LAGRANGE_LEFT)
    CALL LAGRANGE_INTERPOLATING_POLYNOMIAL(N, 1.0D0, GL_POINT, BARY, LAGRANGE_RIGHT)
    !-------------------------------------------------------------------
    
!    PRINT *, LAGRANGE_LEFT
    
    ! COMPUTE DERIVATIVE TRANSPOSE--------------------------------------
    DO J=0, N
        DO I=0, N
            DER_TRANSPOSE(I, J) = - DER(J, I)*GL_W(J)/GL_W(I)
        ENDDO
    
    ENDDO
    !-------------------------------------------------------------------
    

    
    
END SUBROUTINE CONSTRUCT_DG

SUBROUTINE DG_DERIVATIVE(N, SOLUTION_LEFT, SOLUTION_RIGHT, SOLUTION_DER)
!-----------------------------------------------------------------------
! FIRST SPATIAL DERIVATIVE VIA THE GALERKIN APPROXIMATION
! ALGORITHM 60
!----------------------------------------------------------------------- 

    IMPLICIT NONE
    
    INTEGER :: N    ! POLY ORDER
    INTEGER :: J
    
    DOUBLE PRECISION :: SOLUTION_LEFT   ! SOLUTION ON THE LEFT BOUNDARY
    DOUBLE PRECISION :: SOLUTION_RIGHT   ! SOLUTION ON THE RIGHT BOUNDARY
    DOUBLE PRECISION :: SOLUTION_DER(0:N)
    
    CALL MATRIX_VECTOR_DERIVATIVE(N, DER_TRANSPOSE, SOLUTION, SOLUTION_DER)
    

    DO J=0, N
        SOLUTION_DER(J) = SOLUTION_DER(J)+&
                        (SOLUTION_RIGHT*LAGRANGE_RIGHT(J) - &
                        SOLUTION_LEFT*LAGRANGE_LEFT(J))/GL_W(J)
    ENDDO
    

END SUBROUTINE DG_DERIVATIVE

SUBROUTINE EVALUATE_BOUNDARY(N, T, C, SOLUTION_LEFT, SOLUTION_RIGHT)
!!----------------------------------------------------------------------
!! TIME DERIVATIVE FOR DG
!! ALGORITHM 61 (PART 1)
!! INPUT : T: TIME; C: WAVE SPEED
!!----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER:: N ! POLY ORDER
    
    DOUBLE PRECISION :: T   ! TIME
    DOUBLE PRECISION :: C   ! WAVE SPEED
    DOUBLE PRECISION :: SOLUTION_LEFT, SOLUTION_RIGHT   ! APPROXIMATE SOLUTION ON THE LEFT AND RIGHT END OF THE BOUNDARY


!    PRINT *, SOLUTION
    
    ! GET SOLUTION ON THE LEFT AND RIGHT BOUNDARY-----------------------
    IF(C > 0) THEN
        SOLUTION_LEFT = BOUNDARY_SOLUTION(T)
        CALL INTERPOLATE_TO_BOUNDARY(N, SOLUTION, LAGRANGE_RIGHT, SOLUTION_RIGHT)
    ELSE 
        SOLUTION_RIGHT = BOUNDARY_SOLUTION(T)
        CALL INTERPOLATE_TO_BOUNDARY(N, SOLUTION, LAGRANGE_LEFT, SOLUTION_LEFT)
    ENDIF
    !-------------------------------------------------------------------
    
!    PRINT *, LAGRANGE_RIGHT
    
!    PRINT *, SOLUTION_RIGHT

!    PRINT *, SOLUTION

END SUBROUTINE EVALUATE_BOUNDARY

SUBROUTINE DG_TIME_DERIVATIVE(N, SOLUTION_LEFT, SOLUTION_RIGHT, C, &
                                SOLUTION_TIME_DER)
!-----------------------------------------------------------------------
! TIME DERIVATIVE
! ALORITHM 61 (PART 2)
!-----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER :: N
    
    DOUBLE PRECISION :: SOLUTION_LEFT, SOLUTION_RIGHT
    DOUBLE PRECISION :: SOLUTION_TIME_DER(0:N)
    DOUBLE PRECISION :: C   ! WAVE SPEED
    DOUBLE PRECISION :: SOLUTION_DER(0:N)
    
    CALL DG_DERIVATIVE(N, SOLUTION_LEFT, SOLUTION_RIGHT, SOLUTION_DER)
    
    SOLUTION_TIME_DER = -C * SOLUTION_DER

END SUBROUTINE DG_TIME_DERIVATIVE


SUBROUTINE INTERPOLATE_TO_BOUNDARY(N, SOLUTION1, LAGRANGE1, INTERPOLATE_VALUE)
!!----------------------------------------------------------------------
!! TO COMPUTE THE SOLUTION AT THE OTHER BOUNDARY
!! ALORITHM 61 (PART 3)
!!----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER :: N    ! POLY ORDER
    INTEGER :: J
    
    DOUBLE PRECISION :: INTERPOLATE_VALUE
    DOUBLE PRECISION, DIMENSION(0:N) :: SOLUTION1 
    DOUBLE PRECISION, DIMENSION(0:N) :: LAGRANGE1
    
    INTERPOLATE_VALUE = 0.0D0
    
    DO J=0, N
        INTERPOLATE_VALUE = INTERPOLATE_VALUE + LAGRANGE1(J)*SOLUTION1(J)
    
    ENDDO
    

END SUBROUTINE INTERPOLATE_TO_BOUNDARY

SUBROUTINE DG_STEP_BY_RK3(TN, N)
! ----------------------------------------------------------------------
! THRID ORDER RUNGE-KUTTA INTEGRATION 
! INPUT: TN: CURRENT TIME; TIME_DER: DG TIME DERIVATIVE
! ALGORITHM 62
!-----------------------------------------------------------------------

    IMPLICIT NONE
    
    INTEGER :: N
    INTEGER :: I, J
    
    DOUBLE PRECISION :: TN  ! CURRENT TIME
    DOUBLE PRECISION :: T
    DOUBLE PRECISION :: SOLUTION_TIME_DER(0:N)   ! DG TIME DERIVATIVE AT CURRENT TIME
    DOUBLE PRECISION :: G(0:N)  ! INTERMEDIATE ARRAY
    DOUBLE PRECISION :: AM(3), BM(3), GM(3) ! COEFFICIENT
    DOUBLE PRECISION :: SOLUTION_LEFT, SOLUTION_RIGHT

    
    AM = (/0.0D0, (-5.0D0/9.0D0), (-153D0/128.0D0)/)
    BM = (/0.0D0, (1.0D0/3.0D0), (3.0D0/4.0D0)/)
    GM = (/(1.0D0/3.0D0), (15.0D0/16.0D0), (8.0D0/15.0D0)/)
                                    
    SOLUTION_TIME_DER = 0.0D0
    
    DO I=1,3
        T=TN+BM(I)*DELTA_T
!        PRINT *, T
        
        CALL EVALUATE_BOUNDARY(N, T, C, SOLUTION_LEFT, SOLUTION_RIGHT) 
        CALL DG_TIME_DERIVATIVE(N, SOLUTION_LEFT, SOLUTION_RIGHT, C, SOLUTION_TIME_DER)
        
!        PRINT *, SOLUTION_LEFT
        
        DO J=0,N
            G(J) = AM(I)*G(J) + SOLUTION_TIME_DER(J)
            SOLUTION(J) = SOLUTION(J) + GM(I)*DELTA_T*G(J)
        ENDDO
    ENDDO
    


END SUBROUTINE DG_STEP_BY_RK3


SUBROUTINE DECIDE_TIME_STEP(N)
    
    USE PARAM, ONLY: NT, T_TOTAL

    IMPLICIT NONE
    
    INTEGER :: N    ! POLY ORDER
    
    DOUBLE PRECISION :: DEL_T_MAX   ! MAXIMUM DELTA_T FOR STABILITY
    
    !-------------------------------------------------------------------
    IF(N < 32) THEN
        DEL_T_MAX=DBLE(2.5D0/N)
    
    ELSE
        DEL_T_MAX=DBLE(38.0D0*2.51D0/(N)**2)
    ENDIF
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    DELTA_T = DBLE(T_TOTAL/NT)
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    IF(DEL_T_MAX < DELTA_T) THEN

        PRINT *, "TIME STEP NUMBER TN IS TOO SMALL."
        PRINT *, "TN >= ", INT(T_TOTAL/DEL_T_MAX)
        
        STOP
                
    ENDIF
    !-------------------------------------------------------------------

END SUBROUTINE DECIDE_TIME_STEP

SUBROUTINE GET_ERROR(N, ERROR, EXACT, T)

    IMPLICIT NONE
    
    INTEGER :: N
    INTEGER :: I
    
    DOUBLE PRECISION :: EXACT(0: N)
    DOUBLE PRECISION :: ERROR(0: N)
    DOUBLE PRECISION :: T   ! TIME T
    
    ERROR = 0.0D0
    
    CALL EXACT_SOLUTION(T, GL_POINT, N, EXACT)
    
    DO I=0, N
        ERROR(I) = EXACT(I) - SOLUTION(I)
    
    ENDDO


END SUBROUTINE GET_ERROR

SUBROUTINE WRITE_DATA(N, EXACT, TN)

    INTEGER :: N
    INTEGER :: K
    
    DOUBLE PRECISION :: EXACT(0:N)  ! EXACT SOLUTION
    DOUBLE PRECISION :: TN

    OPEN(UNIT=7, FILE="advaction_solution.dat")
    WRITE(7, FMT='(''VARIABLES = "X", "SOLUTION", "EXACT"'')')
    WRITE(7, 60) 'ZONE DATAPACKING=BLOCK, T= ', '"', TN, '"', 'I=',N+1
60 FORMAT(A27, A1, F5.2, A1, 2X, A2, I3)

    ! WRITE SOLUTION----------------------------------------------------
    DO K=0, N
        WRITE(7, 40, ADVANCE='NO') GL_POINT(K)
    ENDDO
40 FORMAT(F10.5)
    
    DO K=0, N
        WRITE(7, 50, ADVANCE='NO') SOLUTION(K)
    ENDDO
50 FORMAT(E20.10)

    DO K=0, N
        WRITE(7, 50, ADVANCE='NO') EXACT(K)
    ENDDO
    !-------------------------------------------------------------------


    CLOSE(UNIT=7)



END SUBROUTINE WRITE_DATA

END MODULE DG_WAVE
