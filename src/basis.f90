MODULE BASIS

USE MPI
USE PARAM, ONLY: N, M

IMPLICIT NONE

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GLL_POINT, GLL_W  ! GAUSS LEGENDRE LOBATTO POINTS AND WEIGHTS
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GL_POINT, GL_W  ! GAUSS LEGENDRE POINTS AND WEIGHTS
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BARY     ! BARYCENTRIC WEIGHTS
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: DER              ! DERIVATIVE MATRIX
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: G     ! INV(MASS MATRIX) * (STIFFNESS MATRIX)
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: LAGRANGE_LEFT   ! LAGRANGE POLY AT LEFT END POINT -1
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: LAGRANGE_RIGHT  ! LAGRANGE POLY AT RIGHT END POINT +1


CONTAINS

SUBROUTINE MATRIX_VECTOR_DERIVATIVE(N, D, F, DER1)
!!----------------------------------------------------------------------
!! COMPUTE THE DERIVATIVE BY MATRIX MULTIPLICATION
!! ALGORITHM 19
!!----------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: N    ! POLY ORDER
    INTEGER :: I, J
    
    DOUBLE PRECISION :: D(0:N, 0:N) ! DERVATIVE MATRIX
    DOUBLE PRECISION :: F(0:N)
    DOUBLE PRECISION :: T
    DOUBLE PRECISION :: DER1(0:N)
    
    !-------------------------------------------------------------------
    DO I=0, N
        
        T=0.0D0
        
        DO J=0, N
            T =T + D(I,J)*F(J)
        ENDDO
        DER1(I) = T
        
    ENDDO
    !-------------------------------------------------------------------
    
    
END SUBROUTINE MATRIX_VECTOR_DERIVATIVE

SUBROUTINE LEGENDRE_POLYNOMIAL_AND_DERIVATIVE(N,X,Q,DQ)
!!--------------------------------------------------------------------
!!       Algorithm 22 Kopriva - compute L_n (Q) and L_n' (DQ)
!!
!! Legendre polynomial of degree k and its derivative using the three
!!                        term recursive.
!!--------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: N
    INTEGER :: K 
    
    DOUBLE PRECISION, INTENT(IN) :: X   ! GLL POINT 
    DOUBLE PRECISION :: Q   ! LEGENDRE POLY OF DEGREE K
    DOUBLE PRECISION :: DQ  ! DERIVETIVE OF LEGENDRE POLY 
    DOUBLE PRECISION :: Q_M1, Q_M2
    DOUBLE PRECISION :: DQ_M1, DQ_M2

    IF (N==0) THEN
        Q = 1.0D0
        DQ = 0.0D0
    ELSEIF (N==1) THEN
        Q = X
        DQ = 1.0D0
    ELSE
        Q_M2 = 1.0D0
        Q_M1 = X
        DQ_M2 = 0.0D0
        DQ_M1 = 1.0D0

        DO K = 2,N
            Q = DBLE(2*K-1)*X*Q_M1/DBLE(K) - DBLE(K-1)*Q_M2/DBLE(K)
            DQ = DQ_M2 + DBLE(2*K-1)*Q_M1
            Q_M2 = Q_M1
            Q_M1 = Q
            DQ_M2 = DQ_M1
            DQ_M1 = DQ
        END DO
    END IF

END SUBROUTINE LEGENDRE_POLYNOMIAL_AND_DERIVATIVE

SUBROUTINE GL(N)
!!----------------------------------------------------------------------
!! COMPUTE THE GAUSS LEGENDRE NODES AND WEIGHTS
!! ALGORITHM 23
!!----------------------------------------------------------------------
    USE PARAM, ONLY : PI
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: N    ! POLYNOMIAL ORDER
    INTEGER :: J
    
    DOUBLE PRECISION :: DELTA
    DOUBLE PRECISION :: Q, DQ, TOL
    
    TOL = 4.0D0*EPSILON(1.0D0)
    Q = 0.0D0
    DQ = 0.0D0
    
    ALLOCATE(GL_POINT(0:N), GL_W(0:N))
    
    IF(N==0)THEN
        GL_POINT(0)=0.0D0
        GL_W(0)=2.0D0
    ELSEIF (N==1) THEN
        GL_POINT(0)=-DSQRT(1.0D0/3.0D0)
        GL_W(0)=1.0D0
        GL_POINT(1)= DSQRT(1.0d0/3.0d0)
        GL_W(1)=1.0D0  
    ELSE
        DO J=0,((N+1)/2) -1
            !Initial guess
            GL_POINT(J)=-DCOS(PI*DBLE(2*J+1)/DBLE(2*N+2))
              
            !Iterative method
            DELTA=1.0D30
            DO WHILE(DABS(DELTA).GE.TOL*DABS(GL_POINT(J)))
                CALL LEGENDRE_POLYNOMIAL_AND_DERIVATIVE(N+1,GL_POINT(J),Q,DQ)   ! ALGORITHM 22
                DELTA=-Q/DQ
                GL_POINT(J)=GL_POINT(J)+DELTA
            END DO
              
            CALL LEGENDRE_POLYNOMIAL_AND_DERIVATIVE(N+1,GL_POINT(J),Q,DQ)
            GL_POINT(N-J)=-GL_POINT(J)
            GL_W(J)=2.0D0/((1.0d0 - GL_POINT(J)**2)*(DQ**2))
            GL_W(N-J) = GL_W(J)
        END DO
    END IF

    IF(MOD(N,2)==0)THEN
        CALL LEGENDRE_POLYNOMIAL_AND_DERIVATIVE(N+1,0.0d0,Q,DQ)
        GL_POINT(N/2)=0.0D0
        GL_W(N/2)=2.0D0/(DQ**2)
    END IF


END SUBROUTINE GL


SUBROUTINE q_And_L_Evaluation(N, X,POL,DPOL,Q,DQ)
!!--------------------------------------------------------------------
!!       Algorithm 24 Kopriva - combined algorithm to compute L_n (POL) 
!!       L'_n (DPOL), Q = L_n+1 - L_n-1, and Q'
!!
!!--------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: N
    INTEGER :: K
    
    DOUBLE PRECISION, INTENT(IN) :: X
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: VAL,DVAL   ! L_n and L'_n
    DOUBLE PRECISION :: POL,DPOL,Q,DQ
    
    ! ALLOCATE----------------------------------------------------------
    ALLOCATE(VAL(0:N+1), DVAL(0:N+1))
    !-------------------------------------------------------------------

    K=2
    VAL(K-2)=1.0D0
    VAL(K-1)=X
    DVAL(K-2)=0.0D0
    DVAL(K-1)=1.0D0

    DO K=2,N+1
        VAL(K) = DBLE(2*K-1)/DBLE(K)*X*VAL(K-1)-DBLE(K-1)/DBLE(K)*VAL(K-2)
        DVAL(K)= DVAL(K-2)+DBLE(2*K-1)*VAL(K-1)
    ENDDO
    
    POL=VAL(N)
    DPOL=DVAL(N)

    Q=VAL(N+1)-VAL(N-1)
    DQ=DVAL(N+1)-DVAL(N-1)
    
    ! DEALLOCATE--------------------------------------------------------
    DEALLOCATE(VAL, DVAL)
    !-------------------------------------------------------------------

END SUBROUTINE q_And_L_Evaluation


SUBROUTINE GLL(N)
!!----------------------------------------------------------------------
!!
!!       Compute Gauss Legendre Lobatto nodes and weights. 
!!       Algorithm taken from algorithm 25 in Kopriva 
!!----------------------------------------------------------------------
    USE PARAM, ONLY: PI
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: N ! POLY ORDER
    INTEGER :: J, ITER

    DOUBLE PRECISION :: DELTA,VAL,DVAL,Q,DQ
    DOUBLE PRECISION :: FACTOR
    DOUBLE PRECISION,PARAMETER:: TOL=1.0D-12


    ! ALLOCATE----------------------------------------------------------
    ALLOCATE(GLL_POINT(0:N))
    ALLOCATE(GLL_W(0:N))
    !-------------------------------------------------------------------

    ! GLL POINTS AND WEIGHTS--------------------------------------------
    IF(N==1)THEN
        GLL_POINT(0)=-1.0D0; GLL_W(0)=1.0D0
        GLL_POINT(1)= 1.0D0; GLL_W(1)=GLL_W(0)
    ELSE
        GLL_POINT(0)=-1.0D0; GLL_W(0)=2.0D0/(DBLE(N*(N+1)))
        GLL_POINT(N)= 1.0D0; GLL_W(N)=GLL_W(0)
        
        DO J=1,(N+1)/2-1
        
            ! TO FIND THE QUADRATURE POINTS, USE NEWTON METHOD----------
            ! Initial guess
            FACTOR=PI*(DBLE(J)+0.25D0)
            GLL_POINT(J)=-DCOS(FACTOR/DBLE(N)-3.0D0/8.0D0/DBLE(N)/FACTOR)
            !-----------------------------------------------------------
      
            !Iterative method
            DELTA=1.0D30
            ITER=0
            DO WHILE(DABS(DELTA).GE.TOL*DABS(GLL_POINT(J)))
                ITER=ITER+1
                CALL q_And_L_Evaluation(N,GLL_POINT(J),VAL,DVAL,Q,DQ)
                DELTA=-Q/DQ
                GLL_POINT(J)=GLL_POINT(J)+DELTA
            ENDDO

            CALL q_And_L_Evaluation(N,GLL_POINT(J),VAL,DVAL,Q,DQ)
            GLL_POINT(N-J)=-GLL_POINT(J)
            GLL_W(J)=2.0D0/(DBLE(N*(N+1))*VAL**2)
            GLL_W(N-J)=GLL_W(J)
        ENDDO
    ENDIF

    IF(MOD(N,2)==0)THEN
        CALL q_And_L_Evaluation(N,0.0D0,VAL,DVAL,Q,DQ)
        GLL_POINT(N/2)=0.0D0
        GLL_W(N/2)=2.0D0/(DBLE(N*(N+1))*VAL**2)
    ENDIF
    !-------------------------------------------------------------------
    
!    PRINT *, GLL_POINT

END SUBROUTINE GLL


SUBROUTINE BARW(N, X)
!!--------------------------------------------------------------------
!!       Barycentric weights
!!       Algorithm 30 in Kopriva
!!-------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: N
    INTEGER :: I,J,K
    
    DOUBLE PRECISION :: X(0:N)    ! SPECTRAL POINTS
    
    ! ALLOCATE----------------------------------------------------------
    ALLOCATE(BARY(0:N))    ! BARYCENTRIC WEIGHTS
    !-------------------------------------------------------------------

    DO I=0,N
        BARY(I)=1.0D0
    ENDDO

    DO J=1,N
        DO K=0,J-1
            BARY(K)=BARY(K)*(X(K)-X(J))
            BARY(J)=BARY(J)*(X(J)-X(K))
        ENDDO
    ENDDO
    
    
    DO J=0,N
        BARY(J)=1.0D0/BARY(J)
    ENDDO
    
!    PRINT *, BARY
    
END SUBROUTINE BARW

SUBROUTINE LAGRANGE_INTERPOLATING_POLYNOMIAL(N, POINT_X, X, BARY, LAGRANGE)
!!---------------------------------------------------------------------
!! LAGRANGE INTERPOLATING POLYNOMIAL VALUES AT POINT X
!! INPUT : GL_POINTS AND BARYCENTRIC WEIGHTS
!! ALGORITHM 34
!!----------------------------------------------------------------------
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: N    ! PLOY ORDER
    INTEGER :: J
    
    DOUBLE PRECISION :: POINT_X ! THE AIMED POINT
    DOUBLE PRECISION :: X(0:N)  ! GL POINTS
    DOUBLE PRECISION :: BARY(0:N) ! BARYCENTRIC WEIGHTS
    DOUBLE PRECISION :: S, T
    DOUBLE PRECISION, DIMENSION(0:N) :: LAGRANGE ! LAGRANGE INTERPOLATING POLYNOMIAL VALUE AT POINT X
    
    LOGICAL :: X_MATCHES_NODE = .FALSE.
    LOGICAL :: FLAG = .FALSE.
    
    LAGRANGE=0.0D0
    S=0.0D0
    
    DO J=0, N
        LAGRANGE(J) = 0.0D0
        
        CALL ALMOSTEQUAL(FLAG, POINT_X, X(J))
        
        IF(FLAG) THEN
            LAGRANGE(J) = 1.0D0
            X_MATCHES_NODE = .TRUE.
        ENDIF
    ENDDO
    
    IF (X_MATCHES_NODE) RETURN
    
    DO J=0, N
        T = BARY(J) / (POINT_X - X(J))
        LAGRANGE(J) = T
        S = S+T
    ENDDO
    
    DO J=0, N
        LAGRANGE(J) = LAGRANGE(J)/S
    
    ENDDO
    

END SUBROUTINE LAGRANGE_INTERPOLATING_POLYNOMIAL

SUBROUTINE mth_Order_Polynomial_Derivative_Matrix(N,MTH_DER, X)
!!----------------------------------------------------------------------
!!       M-th order derivative matrix
!!       Algorithm 37-38 in Kopriva
!!--------------------------------------------------------------------        
    IMPLICIT NONE

    INTEGER :: I,J,K
    INTEGER, INTENT(IN) :: MTH_DER ! M-TH ORDER POLY DERIVATIVE
    INTEGER, INTENT(IN) ::N
    
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: AUX
    DOUBLE PRECISION :: X(0:N)  ! SPECTRAL POINTS
    
    ! ALLOCATE----------------------------------------------------------
    ALLOCATE(DER(0:N, 0:N))
    !-------------------------------------------------------------------
    
    ! FIRST DERIVATIVE APPROXIMATION MATRIX-----------------------------
    CALL BARW(N, X)    ! GET BARYCENTRIC WEIGHTS
    
    DER=0.0D0
    
    DO I=0,N
        DER(I,I)=0.0D0
        DO J=0,N
            IF(J /= I)THEN
                DER(I,J)=BARY(J)/BARY(I)/(X(I)-X(J))
                DER(I,I)=DER(I,I)-DER(I,J)
            ENDIF
        ENDDO
    ENDDO
    !-------------------------------------------------------------------
    
    
    ! IF MTH_DER == 1 STOP----------------------------------------------
    IF(MTH_DER==1) RETURN
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    ALLOCATE(AUX(0:N, 0:N))
    !-------------------------------------------------------------------
    
    ! IF MTH_DER > 1 ---------------------------------------------------
    AUX=DER
    DO K=2,MTH_DER
        DO I=0,N
            DER(I,I)=0.0D0
            DO J=0,N
                IF(J /= I)THEN
                    DER(I,J)=DBLE(K)/(X(I)-X(J))*(BARY(J)/BARY(I)*AUX(I,I)-AUX(I,J))
                    DER(I,I)=DER(I,I)-DER(I,J)
                ENDIF
            ENDDO
        ENDDO
    ENDDO
    !-------------------------------------------------------------------
    
    ! DEALLOCATE--------------------------------------------------------
    DEALLOCATE(AUX)
    !-------------------------------------------------------------------

END SUBROUTINE mth_Order_Polynomial_Derivative_Matrix

SUBROUTINE CG_DERIVATIVE_MATRIX(N)
!!----------------------------------------------------------------------
!!  MATRIX FOR LEGENDRE GALERKIN APPROXIMATION
!!  G = INV(MASS MATRIX) * (STIFFNESS MATRIX)
!!  ALGORITHM 57
!!  NEED FIRST DERIVATIVE MATRIX
!!----------------------------------------------------------------------

    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: N
    INTEGER :: I, J, K
    
    DOUBLE PRECISION :: S
    
    ! COMPUTE G---------------------------------------------------------
    ALLOCATE(G(0:N, 0:N))
    
    DO J=0, N
        DO I=0, N
            S=0.0D0
            
            DO K=0,N
                S = S + DER(K, I) * DER(K, J) * GLL_W(K)
            ENDDO
            G(J,I) = S
        ENDDO
    
    ENDDO
    !-------------------------------------------------------------------
    

END SUBROUTINE CG_DERIVATIVE_MATRIX

SUBROUTINE ALMOSTEQUAL(FLAG,A,B)
!!----------------------------------------------------------------------
!! TESTING EQUALITY OF TWO FLOATING POINT NUMBER
!! ALGORITHM 139
!!----------------------------------------------------------------------
    IMPLICIT NONE

    LOGICAL FLAG
    DOUBLE PRECISION A,B

    IF ((A.EQ.0.0d0).OR.(B.EQ.0.0d0)) THEN
        IF (DABS(A-B).LE.2.0d0*EPSILON(A)) THEN
            FLAG = .TRUE.
        ELSE
            FLAG = .FALSE.
        ENDIF
    ELSE
        IF ((DABS(A-B).LE.EPSILON(A)*DABS(A)).AND.(DABS(A-B).LE.EPSILON(A)*DABS(B))) THEN
            FLAG = .TRUE.
        ELSE
            FLAG = .FALSE.
        ENDIF
    ENDIF

END SUBROUTINE ALMOSTEQUAL



END MODULE BASIS
