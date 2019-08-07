
MODULE GRAPH_PARTITION

USE MPI
USE MESH
USE PARAM, ONLY: RANK, NUM_PROC

IMPLICIT NONE

INTEGER, ALLOCATABLE, DIMENSION(:) :: GROUP ! GROUP(ELEMENT_I) = NUM_OF_PROC

INTEGER, ALLOCATABLE, DIMENSION(:, :) :: DUAL_LAPLACIAN

CONTAINS

SUBROUTINE BFS
!-----------------------------------------------------------------------
! BRESDTH-FIRST-SEARCH
! (1) BFS ASSIGNS EVERY VERTEX TO A LEVEL
! (2) STOPS WHEN ABOUT 1/2VERTICES ARE VISITED
! NUM OF PARITION SHOULD BE 2^N
!-----------------------------------------------------------------------
    INTEGER :: PARTITION_NUM    ! NUM OF PARTIITON
    INTEGER :: SEARCH_ELEM_NUM   ! SEARCH FOR 1/2VERTICES
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ELEM_INDEX  ! THE INDEX OF THE ELEM IN THE SAME GROUP
    INTEGER :: ILAYER, I, INEIGHBOR, J, IROW
    INTEGER :: ELEM_FIND, ELEM_SUBTRACT
    INTEGER :: IGROUP
    
    ! DO THE DOMAIN PARTITION ON PROC 1---------------------------------
    IF (RANK == 0) THEN
    
        ! GENERATE THE DUAL GRAPH LAPLACIAN-----------------------------
        CALL GEN_DUAL_LAPLACIAN
        !---------------------------------------------------------------
        
        ! --------------------------------------------------------------
        ALLOCATE(GROUP(NEL_TOTAL))
        !---------------------------------------------------------------
    
        ! NUM OF PARITIOTN  = PROC NUM----------------------------------
!        PARITIOTN_NUM = NUM_PROC
        PARTITION_NUM = 4
        !---------------------------------------------------------------
        
        IF (PARTITION_NUM == 1) THEN
            GROUP = 0
        ELSE    ! PARTITION NUMBER >= 2
            !-----------------------------------------------------------
            SEARCH_ELEM_NUM = NEL_TOTAL/2     ! SUBTRUCT TEH FISRT ELEM
            !-----------------------------------------------------------
            
            ! ALLOCATE--------------------------------------------------
            ALLOCATE(ELEM_INDEX(NEL_TOTAL))
            !-----------------------------------------------------------
            
            !-----------------------------------------------------------
            GROUP = 0
            ELEM_INDEX = 0
            ELEM_FIND = 1
            IGROUP = 1
            IROW = 1
            !-----------------------------------------------------------
            
            !-----------------------------------------------------------
            GROUP(1) = 1
            ELEM_INDEX(1) = 1
            !-----------------------------------------------------------
            
            !-----------------------------------------------------------
            DO ILAYER=2, NEL_TOTAL
                INEIGHBOR = 0
                !-------------------------------------------------------
                DO I=1, NEL_TOTAL
                    IF(DUAL_LAPLACIAN(IROW, I) == -1) THEN
                    
                        INEIGHBOR = INEIGHBOR+1
                        
                        IF (GROUP(I) == 0 ) THEN
!                            LAYER(INEIGHBOR, IEL+1) = I
                            ELEM_FIND = ELEM_FIND+1
                            GROUP(I) = IGROUP 
                            ELEM_INDEX(ELEM_FIND) = I
                        ENDIF
                        
                        IF (INEIGHBOR == DUAL_LAPLACIAN(IROW, IROW)) EXIT
                    ENDIF
                ENDDO
                !-------------------------------------------------------
                
                IROW = ELEM_INDEX(ILAYER)
                
                !-------------------------------------------------------
                IF(ELEM_FIND > SEARCH_ELEM_NUM) THEN
                    ELEM_SUBTRACT = ELEM_FIND-SEARCH_ELEM_NUM
                    DO J=1, ELEM_SUBTRACT
                        GROUP(ELEM_INDEX(ELEM_FIND)) = IGROUP-1
                        ELEM_FIND = ELEM_FIND-1
                    ENDDO
                    
                ENDIF
                !-------------------------------------------------------
                
                IF(ELEM_FIND == SEARCH_ELEM_NUM) THEN
                    EXIT
                ELSEIF(ELEM_FIND > SEARCH_ELEM_NUM) THEN
                    PRINT *, "BUGS IN 'ELEM_FIND'"
                    STOP 
                ENDIF
                
            ENDDO
            !-----------------------------------------------------------
            
            ! THE OTHER HALF--------------------------------------------
            DO I=1, NEL_TOTAL
                IF(GROUP(I) /= IGROUP) THEN
                    GROUP(I) = IGROUP+1
                ENDIF 
            ENDDO
            !-----------------------------------------------------------
        ENDIF
    ENDIF
    !-------------------------------------------------------------------


END SUBROUTINE BFS

SUBROUTINE BFS_NEW

    INTEGER :: PARTITION_NUM    ! NUM OF PARTIITON
    INTEGER :: SEARCH_ELEM_NUM   ! SEARCH FOR 1/2VERTICES
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: ELEM_INDEX  ! THE INDEX OF THE ELEM IN THE SAME GROUP
    INTEGER :: ILAYER, I, INEIGHBOR, J, IROW
    INTEGER :: ELEM_FIND, ELEM_SUBTRACT
    INTEGER :: IGROUP, IPARTITION, IGROUP_OLD
    INTEGER, DIMENSION(NEL_TOTAL, NEL_TOTAL) :: DUAL_LAPLACIAN_NEW
    INTEGER, ALLOCATABLE, DIMENSION(:) :: GROUP_ELEM_NUM
!    INTEGER :: K
    
    LOGICAL :: START = .TRUE.
    
    ! DO THE DOMAIN PARTITION ON PROC 1---------------------------------
    IF (RANK == 0) THEN
    
        ! GENERATE THE DUAL GRAPH LAPLACIAN-----------------------------
        CALL GEN_DUAL_LAPLACIAN
        !---------------------------------------------------------------
        
        ! --------------------------------------------------------------
        ALLOCATE(GROUP(NEL_TOTAL))
        !---------------------------------------------------------------
        
        !---------------------------------------------------------------
        GROUP = 0
        !---------------------------------------------------------------
        
        ! NUM OF PARITIOTN  = PROC NUM----------------------------------
!        PARITIOTN_NUM = NUM_PROC
        PARTITION_NUM = 4
        !---------------------------------------------------------------
        
        IF (PARTITION_NUM == 1) THEN
            GROUP = 1
        ELSE    ! PARTITION NUMBER >= 2
            
            ! ALLOCATE--------------------------------------------------
            ALLOCATE(ELEM_INDEX(NEL_TOTAL, PARTITION_NUM))
            ALLOCATE(GROUP_ELEM_NUM(PARTITION_NUM))
            !-----------------------------------------------------------
        
            !-----------------------------------------------------------
            ELEM_INDEX = 0
            GROUP_ELEM_NUM = 0
            GROUP_ELEM_NUM(1) = NEL_TOTAL
            GROUP(:) = 1
            !-----------------------------------------------------------
            
            !-----------------------------------------------------------
            DO I=1, NEL_TOTAL
                ELEM_INDEX(I, 1) = I
            ENDDO
            !-----------------------------------------------------------
            
            !-----------------------------------------------------------
            DO WHILE(START)
!            DO K=1, 1
                !-------------------------------------------------------
                ELEM_FIND = 1
                !-------------------------------------------------------
                
                IGROUP_OLD = MAXLOC(GROUP_ELEM_NUM, 1)
            
                CALL GEN_SUB_DUAL_GRAPH(IGROUP_OLD, &
                                        DUAL_LAPLACIAN_NEW)
                                        
                
                
                IROW = ELEM_INDEX(1, IGROUP_OLD)
                
                IGROUP = MAXVAL(GROUP) + 1
                
                ELEM_INDEX(1, IGROUP) = IROW
                
                GROUP(IROW) = IGROUP
                
                SEARCH_ELEM_NUM = GROUP_ELEM_NUM(IGROUP_OLD)/2
                
!                PRINT *, "IGROUP_OLD", IGROUP_OLD, "IGROUP", IGROUP, "IROW", IROW
!                PRINT *, "SEARCH_ELEM", SEARCH_ELEM_NUM
            
                !-------------------------------------------------------
                DO ILAYER=2, NEL_TOTAL
                    INEIGHBOR = 0
                    !---------------------------------------------------
                    DO I=1, NEL_TOTAL
                        IF(DUAL_LAPLACIAN_NEW(IROW, I) == -1) THEN
                        
                            INEIGHBOR = INEIGHBOR+1
                            
                            IF (GROUP(I) == IGROUP_OLD ) THEN
                                ELEM_FIND = ELEM_FIND+1
                                GROUP(I) = IGROUP 
                                ELEM_INDEX(ELEM_FIND, IGROUP) = I
                            ENDIF
                            
!                            PRINT *, "ELEM_FIND", ELEM_FIND
                            
                            IF (INEIGHBOR == DUAL_LAPLACIAN_NEW(IROW, IROW)) EXIT
                            IF (ELEM_FIND == SEARCH_ELEM_NUM) EXIT
                        ENDIF
                    ENDDO
                    !---------------------------------------------------
                    
                    IROW = ELEM_INDEX(ILAYER, IGROUP)
                    
                    !---------------------------------------------------
                    IF(ELEM_FIND > SEARCH_ELEM_NUM) THEN
!                        PRINT *, "ELEM_FIND > SEARCH_ELEM_NUM"
                        ELEM_SUBTRACT = ELEM_FIND-SEARCH_ELEM_NUM
                        DO J=1, ELEM_SUBTRACT
                            GROUP(ELEM_INDEX(ELEM_FIND, IGROUP)) = IGROUP_OLD
                            ELEM_INDEX(ELEM_FIND, IGROUP) = 0
                            ELEM_FIND = ELEM_FIND-1
                        ENDDO
                        
                    ENDIF
                    !---------------------------------------------------
                    
                    !---------------------------------------------------
                    IF(ELEM_FIND == SEARCH_ELEM_NUM) THEN
                        GROUP_ELEM_NUM(IGROUP) = ELEM_FIND
                        
                        CALL UPDATE_ELEM_INDEX(IGROUP, ELEM_INDEX, &
                                PARTITION_NUM, GROUP_ELEM_NUM, &
                                IGROUP_OLD, ELEM_FIND)
                        
                        
                        GROUP_ELEM_NUM(IGROUP_OLD) = GROUP_ELEM_NUM(IGROUP_OLD) &
                                                    - ELEM_FIND

                        EXIT
                    ELSEIF(ELEM_FIND > SEARCH_ELEM_NUM) THEN
                        PRINT *, "BUGS IN 'ELEM_FIND'"
                        STOP 
                    ENDIF
                    !---------------------------------------------------
                    
                ENDDO
                !-------------------------------------------------------
                
!                IF(MAXVAL(GROUP) == 4) THEN
!                    PRINT *, GROUP_ELEM_NUM
!                ENDIF
                
                !-------------------------------------------------------
                IF(MAXVAL(GROUP) == PARTITION_NUM) EXIT
                !-------------------------------------------------------
                
            ENDDO
            !-----------------------------------------------------------
        ENDIF
        
        !---------------------------------------------------------------
        DEALLOCATE(ELEM_INDEX)
        DEALLOCATE(GROUP_ELEM_NUM)
        !---------------------------------------------------------------
        
    ENDIF
    !-------------------------------------------------------------------
    



END SUBROUTINE BFS_NEW

SUBROUTINE GEN_DUAL_LAPLACIAN
!-----------------------------------------------------------------------
! GENERATE THE DUAL GRAPH LAPLACIAN
!-----------------------------------------------------------------------

    INTEGER :: IEL, J
    
    !-------------------------------------------------------------------
    ALLOCATE(DUAL_LAPLACIAN(NEL_TOTAL, NEL_TOTAL))
    
    DUAL_LAPLACIAN = 0
    !-------------------------------------------------------------------

    ! ADJECENCY---------------------------------------------------------
    DO IEL=1, NEL_TOTAL
        DO J=1, 4
            IF(NEIGHBOR_ELEM(J, IEL) /= 0) THEN
                DUAL_LAPLACIAN(IEL, NEIGHBOR_ELEM(J, IEL)) = -1
            ENDIF
        ENDDO
        
    ENDDO
    !-------------------------------------------------------------------
    
    ! DEGREE------------------------------------------------------------
    DO IEL=1, NEL_TOTAL
        DUAL_LAPLACIAN(IEL, IEL) = ABS(SUM(DUAL_LAPLACIAN(IEL, :)))
    ENDDO
    !-------------------------------------------------------------------
    

END SUBROUTINE GEN_DUAL_LAPLACIAN

SUBROUTINE GEN_SUB_DUAL_GRAPH(IGROUP_OLD, DUAL_LAPLACIAN_NEW)

    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: IGROUP_OLD
    INTEGER, DIMENSION(NEL_TOTAL, NEL_TOTAL) :: DUAL_LAPLACIAN_NEW
    INTEGER :: I, J
    
    DUAL_LAPLACIAN_NEW = 0

    
    !-------------------------------------------------------------------
    DO I=1, NEL_TOTAL
        IF(GROUP(I) == IGROUP_OLD) THEN

            DO J=1, NEL_TOTAL
                
                IF((DUAL_LAPLACIAN(I, J) == -1) .AND. &
                    (GROUP(J) == IGROUP_OLD)) THEN
                    
                    DUAL_LAPLACIAN_NEW(I, J) = -1

                ENDIF
                
            ENDDO
            
            DUAL_LAPLACIAN_NEW(I, I) = ABS(SUM(DUAL_LAPLACIAN_NEW(I, :)))
            
        ENDIF
        
    ENDDO
    !-------------------------------------------------------------------

END SUBROUTINE GEN_SUB_DUAL_GRAPH

SUBROUTINE UPDATE_ELEM_INDEX(IGROUP, ELEM_INDEX, &
                            PARTITION_NUM, GROUP_ELEM_NUM, IGROUP_OLD, &
                            ELEM_FIND)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: IGROUP, PARTITION_NUM, IGROUP_OLD, ELEM_FIND
    INTEGER, DIMENSION(NEL_TOTAL, PARTITION_NUM) :: ELEM_INDEX
    INTEGER, DIMENSION(PARTITION_NUM) :: GROUP_ELEM_NUM
    INTEGER, ALLOCATABLE, DIMENSION(:) :: REMAIN
    INTEGER :: I, J, IEL, K
    
    LOGICAL :: STAY
    
    ALLOCATE(REMAIN(GROUP_ELEM_NUM(IGROUP_OLD)-ELEM_FIND))
    
    REMAIN=0
    K=1
    
    !-------------------------------------------------------------------
    DO J=1, GROUP_ELEM_NUM(IGROUP_OLD)
    
        IEL=ELEM_INDEX(J, IGROUP_OLD)
        STAY = .FALSE.
        
        DO I=1, GROUP_ELEM_NUM(IGROUP)
            IF(IEL == ELEM_INDEX(I, IGROUP)) THEN
                STAY = .FALSE.
                EXIT
            ELSE 
                STAY = .TRUE.
            ENDIF
        ENDDO
        
        IF(STAY) THEN
            REMAIN(K) = IEL
            K=K+1
        ENDIF
        
    ENDDO
    !-------------------------------------------------------------------
    
    ELEM_INDEX(:, IGROUP_OLD) = 0
    
    DO I=1, K-1
        ELEM_INDEX(I, IGROUP_OLD) = REMAIN(I)
    ENDDO
    
    DEALLOCATE(REMAIN)
    

END SUBROUTINE UPDATE_ELEM_INDEX
    

END MODULE GRAPH_PARTITION
