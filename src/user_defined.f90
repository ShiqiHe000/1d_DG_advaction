
MODULE USER_DEFINED

USE MPI
USE BASIS, ONLY: GL_POINT

IMPLICIT NONE

CONTAINS

SUBROUTINE INITIAL_CONDITION(N, SOLUTION1)
!-----------------------------------------------------------------------
! USER DEFINE INITIAL CONDITION FOR THE WAVE FUNCTION
!-----------------------------------------------------------------------

    IMPLICIT NONE
    
    INTEGER :: N    ! PLOYNOMIAL ORDER
    INTEGER :: I
    
    DOUBLE PRECISION :: X
    DOUBLE PRECISION :: SOLUTION1(0:N)
    DOUBLE PRECISION :: SIGMA = 0.2D0
                
    DO I=0, N
        
        X = GL_POINT(I)
        
        SOLUTION1(I) = DEXP((-DLOG(2.0D0)*(X+1.0D0)**2)/(SIGMA**2))
!       SOLUTION1(I) = 1.0D0
    
    ENDDO

END SUBROUTINE INITIAL_CONDITION

FUNCTION BOUNDARY_SOLUTION(T)
!!----------------------------------------------------------------------
!! SOLUTION ON THE BOUNDARY
!! B.C.
!!----------------------------------------------------------------------
    USE PARAM, ONLY: C

    IMPLICIT NONE
    
    DOUBLE PRECISION :: BOUNDARY_SOLUTION
    DOUBLE PRECISION :: T   ! TIME
    DOUBLE PRECISION :: SIGMA = 0.2D0
    
    BOUNDARY_SOLUTION = DEXP(-DLOG(2.0D0)*((C*T)**2)/(SIGMA**2))
!    BOUNDARY_SOLUTION = 1.0D0
    
    RETURN
    
END FUNCTION BOUNDARY_SOLUTION

SUBROUTINE EXACT_SOLUTION(T, X, N, EXACT_SOLU)

    USE PARAM, ONLY: C

    IMPLICIT NONE
    
    INTEGER :: N    ! POLY ORDER
    INTEGER :: I
    
    DOUBLE PRECISION :: X(0:N)  ! NODES
    DOUBLE PRECISION :: T        ! END TIME
    DOUBLE PRECISION :: EXACT_SOLU(0:N) ! RESULTS
    DOUBLE PRECISION :: SIGMA = 0.2D0
    
    DO I=0, N
        EXACT_SOLU(I) = DEXP(-DLOG(2.0D0)*&
                        (X(I)-C*T+1.0D0)**2/(SIGMA**2))

    ENDDO
    
    
    

END SUBROUTINE EXACT_SOLUTION

END MODULE USER_DEFINED
