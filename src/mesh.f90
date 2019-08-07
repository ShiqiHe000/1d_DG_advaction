

MODULE MESH

USE MPI
USE PARAM, ONLY: MESHFILE, RANK, NUM_PROC

IMPLICIT NONE

INTEGER, ALLOCATABLE, DIMENSION(:) :: ELE_PROC_NUM   ! ELEMENT PROCCESSOR NUMBER
INTEGER :: NELMAX       ! MAX NUMBER OF ELEMENT EACH PROCESSOR
INTEGER :: NEL_LOCAL ! NUMBER OF ELEMENT IN EACH PROC
INTEGER :: NEL_TOTAL                ! TOTAL NUMBER OF ELEMENTS
INTEGER :: NDIM                     ! DIMENSION
INTEGER :: TOTAL_NODE      ! TOTAL NODE NUMBER
INTEGER, ALLOCATABLE, DIMENSION(:, :) :: GRAPH_LAPLACIAN    
INTEGER, ALLOCATABLE, DIMENSION(:, :) :: NEIGHBOR_ELEM
INTEGER, ALLOCATABLE, DIMENSION(:, :) :: NEIGHBOR_ELEM_SIDE
INTEGER, ALLOCATABLE, DIMENSION(:) :: BOUNDARY_LINE_TAG ! PRO-DEFINED BOUNDARY LINE TAG
INTEGER, ALLOCATABLE, DIMENSION(:, :) :: QUAD_NODE  ! QUADRANGLE NODES
   

REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: X_LOCAL, Y_LOCAL ! NODE COORDINATES 
REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: X_GLOBAL, Y_GLOBAL    ! ELEMENT GLOBAL COORD

CHARACTER(1), ALLOCATABLE, DIMENSION(:, :) :: BOUNDARY_TYPE
CHARACTER(1), DIMENSION(:), ALLOCATABLE :: PHYSICAL_NAME

CONTAINS

SUBROUTINE READ_MESH(MESHFILE)
    IMPLICIT NONE
    
    ! VARIABLES---------------------------------------------------------
    INTEGER :: NUM_OF_PHY_NAME  ! NUMBER OF PHYSICAL NAMES
    INTEGER :: NUM_OF_ELEMENT   ! NUMBER OF ELEMENT IN .MSH FILE(LINE + QUAD)
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: BOUNDARY_LINE
    INTEGER :: TOTAL_B_LINE     ! TOTAL BOUNDARY LINES  
    INTEGER :: TOTAL_QUAD   ! TOTAL QUADRANGLE 
    INTEGER :: ELEM_TYPE    ! ELEMENT TYPE
    INTEGER :: I, A, B, C, D, E, J
    INTEGER, ALLOCATABLE, DIMENSION(:) :: PHYSICAL_TAG
    INTEGER :: TAG
!    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: ELEM_EDGE   ! THE NODES NUMBER OF EACH EDGE OF EACH ELEMENT
!    INTEGER, ALLOCATABLE, DIMENSION(:) :: BOUNDARY_LINE_TAG
    
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: NODE_XY ! NODE COORDINATES 
    
    CHARACTER(LEN=*), INTENT(IN) :: MESHFILE
    CHARACTER(LEN=80) :: CHARLINE    ! DUMMY LINE
    CHARACTER(LEN=50) :: LINE   ! STORE ELEMENT LINES
    CHARACTER(LEN=10) :: PHYSICAL_NAME1  ! STORE THE PHYSICAL NAME (BOUNDARY TYPE) OF THE ELEMENT SIDE
    !-------------------------------------------------------------------
    
    ! ONLY RANK0 READ MESHFILE------------------------------------------
    IF(RANK == 0) THEN
    
        PRINT *, "------------------------------------------------------"
        PRINT *, "START TO READ ", MESHFILE
        PRINT *, "------------------------------------------------------"
        
        OPEN (UNIT=1, FILE=MESHFILE, STATUS='OLD')
        
        ! GET DIMENSION-------------------------------------------------
        DO WHILE (.TRUE.)
            READ(1, *) CHARLINE
            CHARLINE = TRIM(CHARLINE)
            IF(CHARLINE == "$PhysicalNames") EXIT
        ENDDO
        
        READ(1, *) NUM_OF_PHY_NAME
        
!        PRINT *, NUM_OF_PHY_NAME
        
        ALLOCATE(PHYSICAL_TAG(NUM_OF_PHY_NAME))   ! NUM_OF_PHY_NAME -1
        ALLOCATE(PHYSICAL_NAME(NUM_OF_PHY_NAME))
        
        ! INITIALIZE----------------------------------------------------
        PHYSICAL_TAG = 0
        PHYSICAL_NAME = "N"
        !---------------------------------------------------------------
        
        ! RECORD THE PHYSICAL_NAME--------------------------------------
        DO I=1, NUM_OF_PHY_NAME
        
            READ(1, *) NDIM, TAG, PHYSICAL_NAME1   !READ PHYSICAL DIMENSION
            
            PHYSICAL_NAME1 = TRIM(PHYSICAL_NAME1)
            
            IF (NDIM == 1) THEN ! IF DIMENSION == 1
                PHYSICAL_TAG(I) = TAG
                
                IF(LEN(TRIM(PHYSICAL_NAME1)) /= 1) THEN
                    PRINT *, "$PhysicalNames in .msh file are in wrong format. "
                    STOP
                ENDIF
                
                PHYSICAL_NAME(TAG) = PHYSICAL_NAME1
            ENDIF
        
        ENDDO
        !---------------------------------------------------------------

        !---------------------------------------------------------------
        READ(1, *)  ! READ DUMMY LINE
        READ(1, *)
        !---------------------------------------------------------------
        
        ! READ TOTAL NODES----------------------------------------------
        READ(1, *) TOTAL_NODE
        !---------------------------------------------------------------
        
        ! READ NODE COORDINATES-----------------------------------------
        IF(NDIM == 2) THEN
        
            ! ALLOCATE--------------------------------------------------
            ALLOCATE(NODE_XY(2, TOTAL_NODE))
            !-----------------------------------------------------------
            
            ! READ------------------------------------------------------
            DO I=1, TOTAL_NODE
                READ(1, *) A, NODE_XY(1, I), NODE_XY(2, I)
                
            ENDDO
            !-----------------------------------------------------------
            
            !-----------------------------------------------------------
            READ(1, *)  ! READ DUMMY LINE, $EndNodes
            READ(1, *)  ! $Elements
            !-----------------------------------------------------------
            
            ! NUM OF ELEMENT--------------------------------------------
            READ(1, *) NUM_OF_ELEMENT
            !-----------------------------------------------------------
            
!            PRINT *, "NUM_OF_ELE", NUM_OF_ELEMENT
            
            !-----------------------------------------------------------
            ALLOCATE(BOUNDARY_LINE(2, NUM_OF_ELEMENT))
            ALLOCATE(BOUNDARY_LINE_TAG(NUM_OF_ELEMENT))
            ALLOCATE(QUAD_NODE(4, NUM_OF_ELEMENT))
            !-----------------------------------------------------------
            
            ! INITIALIZE------------------------------------------------
            BOUNDARY_LINE_TAG = 0
            BOUNDARY_LINE = 0
            TOTAL_B_LINE = 0
            TOTAL_QUAD = 0
            !-----------------------------------------------------------
            
            !-----------------------------------------------------------
            DO I=1, NUM_OF_ELEMENT
                
                READ(1, '(A)') LINE     ! STORE THE WHOLE LINE AS A STRING
                
                LINE = TRIM(LINE)
                
                READ(LINE, *) A, ELEM_TYPE  ! GET THE ELEMENT TYPE
                
                !-------------------------------------------------------
                IF (ELEM_TYPE == 1) THEN    ! 2-NODE LINE
                
                    TOTAL_B_LINE = TOTAL_B_LINE+1
                
                    READ(LINE, *) A, B, C, &
                                BOUNDARY_LINE_TAG(TOTAL_B_LINE), D, &
                                BOUNDARY_LINE(1, TOTAL_B_LINE), &
                                BOUNDARY_LINE(2, TOTAL_B_LINE)
                                

                ELSEIF (ELEM_TYPE == 3) THEN    ! 4-NODE QUAD
                    TOTAL_QUAD = TOTAL_QUAD+1
                    
                    READ(LINE, *) A, B, C, D, E, &
                                QUAD_NODE(1, TOTAL_QUAD), &
                                QUAD_NODE(2, TOTAL_QUAD), &
                                QUAD_NODE(3, TOTAL_QUAD), &
                                QUAD_NODE(4, TOTAL_QUAD)
                    
                            
                ELSE
                    PRINT *, "ELEMENT TYPE TAG INCORRECT"
                    STOP
                ENDIF
                !-------------------------------------------------------
            
            ENDDO
            !-----------------------------------------------------------

            ! STORE GLOBAL ELEMENT NODE COORDINATES---------------------
            NEL_TOTAL = TOTAL_QUAD
            NEL_LOCAL = TOTAL_QUAD
            
            ALLOCATE(X_GLOBAL(4, TOTAL_QUAD), Y_GLOBAL(4, TOTAL_QUAD))
            
            DO I=1, TOTAL_QUAD
            
                DO J=1, 4
                    X_GLOBAL(J, I) = NODE_XY(1, QUAD_NODE(J, I))
                    Y_GLOBAL(J, I) = NODE_XY(2, QUAD_NODE(J, I))
                
                ENDDO
                
            ENDDO
            !-----------------------------------------------------------
            
            ! GENERATE GRAPH LAPLACIAN----------------------------------
            CALL GEN_GRAPH_LAPLACIAN_QUAD(TOTAL_NODE)
            !-----------------------------------------------------------
            
            ! ----------------------------------------------------------
            CALL FIND_NEIGHBOR(TOTAL_B_LINE, BOUNDARY_LINE)
            !-----------------------------------------------------------
            
            ! WRITE .REA FILE-------------------------------------------
            CALL WRITE_REA
            !-----------------------------------------------------------

        ELSE
        
            PRINT *, "THE DIMENSION OF THE MESH FILE DOES NOT FIT."
            STOP
            
        ENDIF
        !---------------------------------------------------------------
        
        CLOSE(UNIT=1)
        
        ! DEALLOCATE----------------------------------------------------
        DEALLOCATE(BOUNDARY_LINE)
        DEALLOCATE(QUAD_NODE)
        DEALLOCATE(PHYSICAL_TAG)
        DEALLOCATE(NODE_XY)
        DEALLOCATE(BOUNDARY_LINE_TAG)
        DEALLOCATE(BOUNDARY_TYPE)
        DEALLOCATE(PHYSICAL_NAME)
        !---------------------------------------------------------------
    
    ENDIF
    !-------------------------------------------------------------------
    
END SUBROUTINE READ_MESH

SUBROUTINE WRITE_REA

    INTEGER :: GROUP, IEL, ISIDE
    
    ! WRITE DATA TO .REA FILE-------------------------------------------
    OPEN(UNIT=3, FILE='mesh.rea')
    
    WRITE(3, *) "**MESH DATA** 1st line is X of corner 1,2,3,4. 2nd line is Y."
    
    WRITE(3 ,FMT='(5X,I4,9X,I4,9X,I4,1X,''NEL,'',''NDIM,'',''NGROUPS'')') &
                NEL_TOTAL,NDIM,NUM_PROC
    
    ! WRITE NODE COORDS-------------------------------------------------
    DO GROUP=1, NUM_PROC
        DO IEL=1, NEL_LOCAL
            WRITE(3,FMT='(5X,''ELEMENT'',5X,I5,'' [    1A]'',3X,''GROUP'',1X,I4)') &
                            IEL, GROUP
            
            WRITE(3,FMT='(1X,F10.7,6X,F10.7,6X,F10.7,6X,F10.7)') &
                    X_GLOBAL(1,IEL),X_GLOBAL(2,IEL),&
                    X_GLOBAL(3,IEL),X_GLOBAL(4,IEL)
            
            WRITE(3,FMT='(1X,F10.7,6X,F10.7,6X,F10.7,6X,F10.7)') &
                    Y_GLOBAL(1,IEL),Y_GLOBAL(2,IEL),&
                    Y_GLOBAL(3,IEL),Y_GLOBAL(4,IEL)

        ENDDO
    
    ENDDO
    !-------------------------------------------------------------------
    
    WRITE(3,FMT='('' ***** CURVED SIDE DATA *****'')')
    
    WRITE(3,FMT='(''       0  Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE '')')
    
    WRITE(3,FMT='('' ***** GROUP CONNECTIVITY ***** '')')
    
    ! WRITE GROUP CONNECTIVITY------------------------------------------
    IF (NEL_TOTAL < 10000) THEN
        DO IEL=1, NEL_TOTAL
            DO ISIDE=1, 4
                
!                WRITE(3,FMT='(1X,A3,I3,I4,I3,5G14.6)')
                WRITE(3,FMT='(A1, 2X, I4, 2X, I1, 2X, I3, 2X, I1)') &
                        BOUNDARY_TYPE(ISIDE, IEL), &
                        IEL, ISIDE, NEIGHBOR_ELEM(ISIDE, IEL), &
                        NEIGHBOR_ELEM_SIDE(ISIDE, IEL)
            ENDDO
        
        ENDDO
    ELSE
        PRINT *, "TOTAL ELEMNT > 10000, TOO LARGE."
        STOP
    ENDIF
    !-------------------------------------------------------------------

    
    
    CLOSE(UNIT=3)
    !-------------------------------------------------------------------

END SUBROUTINE WRITE_REA

SUBROUTINE GEN_GRAPH_LAPLACIAN_QUAD(TOTAL_NODE)
!-----------------------------------------------------------------------
! GENERATE THE GRAPH LAPLACIAN FOR ELEMENT TYPE = QUAD
! GRAPH LAPLACIAN : TOTAL_NODE * TOTAL_NODE
!-----------------------------------------------------------------------

    INTEGER, INTENT(IN) :: TOTAL_NODE
    INTEGER, DIMENSION(8) :: NODE_SEQUENCE
    INTEGER :: IEL, INODE, I
    INTEGER :: NODE1, NODE2
    
    !-------------------------------------------------------------------
    ALLOCATE(GRAPH_LAPLACIAN(TOTAL_NODE, TOTAL_NODE))
    
    GRAPH_LAPLACIAN = 0
    !-------------------------------------------------------------------
    
    ! ADJACENCY MATRIX * -1 --------------------------------------------
    DO IEL=1, NEL_TOTAL
    
        NODE_SEQUENCE(1:4) = QUAD_NODE(:, IEL)
        NODE_SEQUENCE(5:8) = QUAD_NODE(:, IEL)

        DO INODE=1, 4
            
            NODE1 = NODE_SEQUENCE(INODE)
            NODE2 = NODE_SEQUENCE(INODE+1)
            
            IF (GRAPH_LAPLACIAN(NODE1, NODE2) == 0 ) THEN
                GRAPH_LAPLACIAN(NODE1, NODE2) = -1
            ENDIF
            
            NODE2 = NODE_SEQUENCE(INODE+3)
            
            IF (GRAPH_LAPLACIAN(NODE1, NODE2) == 0 ) THEN
                GRAPH_LAPLACIAN(NODE1, NODE2) = -1
            ENDIF
            
        ENDDO
    
    ENDDO
    !-------------------------------------------------------------------
    
    ! ADD DEGREE MATRIX-------------------------------------------------
    DO I=1, TOTAL_NODE
        GRAPH_LAPLACIAN(I,I) = ABS(SUM(GRAPH_LAPLACIAN(I, :)))
    
    ENDDO
    !-------------------------------------------------------------------
    
    ! WRITE GRAPH LAPLACIAN TO FILE-------------------------------------
    OPEN(UNIT=4, FILE='graph_laplacian.txt')
    DO I=1, TOTAL_NODE
        WRITE(4, *) GRAPH_LAPLACIAN(I, :)
    ENDDO
    CLOSE(UNIT=4)
    !-------------------------------------------------------------------
    

END SUBROUTINE GEN_GRAPH_LAPLACIAN_QUAD

SUBROUTINE FIND_NEIGHBOR(TOTAL_B_LINE, BOUNDARY_LINE)
!-----------------------------------------------------------------------
! FINE THE NEIGHBOR OF EACH ELEMENT
! NODE NUMBERLING CAN BE CLOCKWISE OR COUNTER-CLOCKWISE OR MIXED
!-----------------------------------------------------------------------

    INTEGER :: IEL, ISIDE, I, J, IBLINE
    INTEGER, ALLOCATABLE, DIMENSION(:) :: NODE_SEQUENCE
    INTEGER :: INODE, NODE1, NODE2, NODE3, NODE4
    INTEGER, INTENT(IN) :: TOTAL_B_LINE ! TOTAL PRE-DEFINED BOUNDARY LINE
    INTEGER, DIMENSION(2, TOTAL_B_LINE) :: BOUNDARY_LINE
    !-------------------------------------------------------------------
    ALLOCATE(NEIGHBOR_ELEM(4,NEL_TOTAL))  ! (SIDE, NUM_OF_ELEMENT)
    ALLOCATE(NEIGHBOR_ELEM_SIDE(4,NEL_TOTAL))  ! (SIDE, NUM_OF_ELEMENT)
    ALLOCATE(NODE_SEQUENCE(5))
!    ALLOCATE(BOUNDARY_TYPE(4, NEL_TOTAL))
    !-------------------------------------------------------------------
    
    ! INITIALIZE -------------------------------------------------------
    NEIGHBOR_ELEM = 0
    NEIGHBOR_ELEM_SIDE = 0
    !-------------------------------------------------------------------
    
    ! FIND NEIGHBOR AND NEIGHBOR SIDE-----------------------------------
    DO IEL=1, NEL_TOTAL
        !---------------------------------------------------------------
        DO ISIDE=1,4    ! SIDE
        
            INODE = QUAD_NODE(ISIDE, IEL)

            IF(GRAPH_LAPLACIAN(INODE, INODE) < 2 ) THEN
                PRINT *, "BUGS IN GRAPH_LAPLACIAN"
                STOP
            ELSEIF(GRAPH_LAPLACIAN(INODE, INODE) > 2) THEN  ! GRAPH_LAPLACIAN == 2 NO NEIGHBOR
                CALL GET_NODE_SEQUENCE(QUAD_NODE(:, IEL), NODE_SEQUENCE)
                
                NODE1 = NODE_SEQUENCE(ISIDE)
                NODE2 = NODE_SEQUENCE(ISIDE+1)
                ! FIND THE NEIGHBOR ELEMENT NUM AND SIDE NUM------------
                DO I=IEL+1, NEL_TOTAL
                    
                    CALL GET_NODE_SEQUENCE(QUAD_NODE(:, I), NODE_SEQUENCE)
                    
                    DO J=1,4    ! SIDE
                    
                        NODE3 = NODE_SEQUENCE(J)
                        NODE4 = NODE_SEQUENCE(J+1)
                        
                        IF((NODE1 == NODE4 .AND. NODE2 == NODE3) .OR. &
                            (NODE1 ==NODE3 .AND. NODE2 == NODE4)) THEN
                            NEIGHBOR_ELEM(ISIDE, IEL) = I   ! NEIGHBOR ELEMENT NUM
                            NEIGHBOR_ELEM_SIDE(ISIDE, IEL) = J  ! NEIGHBOR ELEM SIDE
                            
                            NEIGHBOR_ELEM(J, I) = IEL
                            NEIGHBOR_ELEM_SIDE(J, I) = ISIDE
                            EXIT
                        ENDIF
                    ENDDO
                ENDDO
                !-------------------------------------------------------
                
            ENDIF
            
        ENDDO
        !---------------------------------------------------------------
    
    ENDDO
    !-------------------------------------------------------------------
    
    ! FIND THE FOUR BOUNDARY TYPE --------------------------------------
    CALL FIND_BOUDARY_TYPE(TOTAL_B_LINE, BOUNDARY_LINE)
    !-------------------------------------------------------------------

    DEALLOCATE(NODE_SEQUENCE)

END SUBROUTINE FIND_NEIGHBOR

SUBROUTINE GET_NODE_SEQUENCE(QUAD_NODE, NODE_SEQUENCE)

    INTEGER, DIMENSION(5) :: NODE_SEQUENCE
    INTEGER, DIMENSION(4) :: QUAD_NODE
    
    NODE_SEQUENCE(1:4) = QUAD_NODE(:)
    
    NODE_SEQUENCE(5) = QUAD_NODE(1)
    
END SUBROUTINE GET_NODE_SEQUENCE

SUBROUTINE FIND_BOUDARY_TYPE(TOTAL_B_LINE, BOUNDARY_LINE)

    INTEGER, INTENT(IN) :: TOTAL_B_LINE
    INTEGER :: IEL, IBLINE, ISIDE
    INTEGER :: NODE1, NODE2, NODE3, NODE4
    INTEGER, DIMENSION(5) :: NODE_SEQUENCE
    INTEGER, DIMENSION(2, TOTAL_B_LINE) :: BOUNDARY_LINE
    
    !-------------------------------------------------------------------
    ALLOCATE(BOUNDARY_TYPE(4, NEL_TOTAL))
    
    BOUNDARY_TYPE = "N" ! NO SPECIFIED BOUNDARY CONDITION
    !-------------------------------------------------------------------
    
    DO IEL=1, NEL_TOTAL
        !---------------------------------------------------------------
        DO ISIDE=1, 4
        
            CALL GET_NODE_SEQUENCE(QUAD_NODE(:, IEL), NODE_SEQUENCE)
            NODE1 = NODE_SEQUENCE(ISIDE)
            NODE2 = NODE_SEQUENCE(ISIDE+1)
            
            DO IBLINE=1, TOTAL_B_LINE
                NODE3 = BOUNDARY_LINE(1, IBLINE)
                NODE4 = BOUNDARY_LINE(2, IBLINE)
                
                IF((NODE1 == NODE3 .AND. NODE2 == NODE4) .OR. &
                    (NODE1 == NODE4 .AND. NODE2 == NODE3)) THEN
                    BOUNDARY_TYPE(ISIDE, IEL) = &
                    PHYSICAL_NAME(BOUNDARY_LINE_TAG(IBLINE))
                    EXIT
                ENDIF
            
            ENDDO
            
            ! NO PHYSICAL TAG == INTERNAL(ELEMENT CONNECTIVITY)
            IF(BOUNDARY_TYPE(ISIDE, IEL) == "N") THEN
                BOUNDARY_TYPE(ISIDE, IEL) = "E"
            ENDIF
            
        ENDDO
        !---------------------------------------------------------------
    ENDDO

END SUBROUTINE FIND_BOUDARY_TYPE




END MODULE MESH
