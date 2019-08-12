!> mainpage 1d_advection
!> @author Shiqi He
!> @version 1.0


PROGRAM MAIN_LOOP
!-----------------------------------------------------------------------
!   POISSON PROBLEM
!-----------------------------------------------------------------------
    USE MPI
    USE BASIS
    USE PARAM, ONLY: IERROR, MESHFILE
    USE SET_MPI
    USE MESH
    USE GRAPH_PARTITION
    USE FIELDS
    USE DG_WAVE

    IMPLICIT NONE
    
!    INTEGER :: I

    
    ! START MPI---------------------------------------------------------
    CALL START_MPI
    !-------------------------------------------------------------------
    
!    CALL READ_MESH(MESHFILE)
    
!    CALL BFS_NEW
    
!    CALL WRITE_TO_FILES
    !-------------------------------------------------------------------
    
    CALL DG_WAVE_PROCEDURE
    

    
    !-------------------------------------------------------------------
    
    ! END MPI-----------------------------------------------------------
    CALL MPI_FINALIZE(IERROR)
    !-------------------------------------------------------------------


END PROGRAM MAIN_LOOP
