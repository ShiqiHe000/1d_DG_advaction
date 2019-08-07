
MODULE FIELDS

USE MESH, ONLY: NEL_TOTAL, X_GLOBAL, Y_GLOBAL
USE GRAPH_PARTITION, ONLY: GROUP
USE PARAM, ONLY: RANK
USE MPI

IMPLICIT NONE

INTEGER :: FRAME = 1

CONTAINS

SUBROUTINE WRITE_TO_FILES

    INTEGER :: IEL
    INTEGER :: ELEM
    
    CHARACTER(LEN=16) :: FILENAME

    ! ONLY PROC 1 WRITE-------------------------------------------------
    IF (RANK ==0 ) THEN
    
        WRITE(FILENAME,FMT='(''aoutput'',I5.5,''.dat'')') FRAME
        
        OPEN(5,FILE=FILENAME)
        
        WRITE(5, FMT='(''TITLE = "MESH AND SOLUTIONS"'')')
        WRITE(5, FMT='(''VARIABLES = "X", "Y", "GROUP"'')')
        
        ELEM=1
        
        DO IEL=1, NEL_TOTAL

            WRITE(5, 20) "ZONE T= ", '"', "IEL", ELEM, '"', &
                        "I=2, J=2","DATAPACKING = POINT" 
20 FORMAT(A8, A1, A3, I6, A1, 2X, A8, 2X, A19)

            ELEM = ELEM+1
            
            
            WRITE(5, 10) X_GLOBAL(1, IEL), Y_GLOBAL(1, IEL), GROUP(IEL)
            WRITE(5, 10) X_GLOBAL(2, IEL), Y_GLOBAL(2, IEL), GROUP(IEL)
            WRITE(5, 10) X_GLOBAL(4, IEL), Y_GLOBAL(4, IEL), GROUP(IEL)
            WRITE(5, 10) X_GLOBAL(3, IEL), Y_GLOBAL(3, IEL), GROUP(IEL)
            
            
        ENDDO
10 FORMAT(F10.5, 2X, F10.5, 2X, I2)

    CLOSE(UNIT=5)

    ENDIF
    !-------------------------------------------------------------------


END SUBROUTINE WRITE_TO_FILES

END MODULE FIELDS
