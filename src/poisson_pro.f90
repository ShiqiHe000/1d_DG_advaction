
MODULE POISSON_PRO

USE MPI
USE BASIS
USE PARAM, ONLY: M, N

IMPLICIT NONE

CONTAINS

SUBROUTINE SOLVE_POISSON

    ! CREATE A NODAL GALERKIN CONSTRUCTOR-------------------------------
    CALL NODAL_POTENTIAL_CLASS
    !-------------------------------------------------------------------


END SUBROUTINE SOLVE_POISSON

SUBROUTINE NODAL_POTENTIAL_CLASS

    IMPLICIT NONE
    
    ! GET GAUSS_LOBBATO_LEGENDRE POINTS AND WEIGHTS---------------------
    CALL GLL(M) ! X
    CALL GLL(N) ! Y
    !-------------------------------------------------------------------
    
    ! GET FIRST DERIVATIVE MATRIX---------------------------------------
    CALL mth_Order_Polynomial_Derivative_Matrix(M,1)    ! X DIRECTION
    CALL mth_Order_Polynomial_Derivative_Matrix(N,1)    ! Y DIRECTION
    !-------------------------------------------------------------------
    
    ! GET G MATRIX (INV(MASS) * (STIFFNESS))----------------------------
    CALL CG_DERIVATIVE_MATRIX(M)    ! X
    CALL CG_DERIVATIVE_MATRIX(N)    ! Y
    !-------------------------------------------------------------------


END SUBROUTINE NODAL_POTENTIAL_CLASS

SUBROUTINE LAPLACE_ON_THE_SQUARE

 IMPLICIT NONE
 
 INTEGER :: I, J
 
 
 


END SUBROUTINE LAPLACE_ON_THE_SQUARE



END MODULE POISSON_PRO
