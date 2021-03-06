module petsc_data_structure_mod
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petsclog.h>


    use petscvec
    use parameter_mod
    use data_structure_mod

    implicit none

    PetscMPIInt      :: rank,proc
    Vec          :: p_x, p_y
    Vec          :: p_dq
    Vec          :: p_qm
    Vec          :: p_prim
    PetscLogEvent    :: dq_comm, prim_comm, qm_comm, x_comm, y_comm

    contains

    subroutine init_petsc()
        implicit none
        PetscErrorCode       :: ierr
        if(rank == 0) then
            write(*,*) '%%%%%%%-Setting up parallel vectors-%%%%%%%'
            write(*,*)
        end if
        pghost = pghost - 1

        call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,2*4,2*4*local_points,&
            &PETSC_DECIDE,ghost_points,pghost,point%dq(1,1,1),p_dq,ierr)
        
        call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,2*4,2*4*local_points,&
            &PETSC_DECIDE,ghost_points,pghost,point%qm(1,1,1),p_qm,ierr)

        call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,4,4*local_points,&
            &PETSC_DECIDE,ghost_points,pghost,point%prim(1,1),p_prim,ierr)

        call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,1,local_points,&
            &PETSC_DECIDE,ghost_points,pghost,point%x(1),p_x,ierr)

        call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,1,local_points,&
            &PETSC_DECIDE,ghost_points,pghost,point%y(1),p_y,ierr)

        call VecGetSize(p_prim,plen,ierr)
        plen = plen/4

        call PetscLogEventRegister('dq_comm',  0,dq_comm,ierr);
        call PetscLogEventRegister('qm_comm',  0,qm_comm,ierr);
        call PetscLogEventRegister('prim_comm',  0,prim_comm,ierr);
        call PetscLogEventRegister('x_comm',  0, x_comm,ierr);
        call PetscLogEventRegister('y_comm',  0, y_comm,ierr);

    end subroutine 

    subroutine dest_petsc()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return
            
        call VecDestroy(p_dq,ierr)
        call VecDestroy(p_qm,ierr)
        call VecDestroy(p_prim,ierr)
        call VecDestroy(p_x,ierr)
        call VecDestroy(p_y,ierr)
    end subroutine 


    subroutine update_begin_dq_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return

        call VecGhostUpdateBegin(p_dq,INSERT_VALUES,SCATTER_FORWARD,ierr)

    end subroutine 

    subroutine update_begin_prim_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return
 
        call VecGhostUpdateBegin(p_prim,INSERT_VALUES,SCATTER_FORWARD,ierr)

    end subroutine 

    subroutine update_begin_qm_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return

        call VecGhostUpdateBegin(p_qm,INSERT_VALUES,SCATTER_FORWARD,ierr)

    end subroutine 

    subroutine update_end_dq_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return

        call PetscLogEventBegin(dq_comm, ierr)
        call VecGhostUpdateEnd(p_dq,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call PetscLogEventEnd(dq_comm, ierr)

    end subroutine 

    subroutine update_end_qm_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return

        call PetscLogEventBegin(qm_comm, ierr)
        call VecGhostUpdateEnd(p_qm,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call PetscLogEventEnd(qm_comm, ierr)

    end subroutine 
    
    subroutine update_end_prim_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return

        call PetscLogEventBegin(prim_comm, ierr)
        call VecGhostUpdateEnd(p_prim,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call PetscLogEventEnd(prim_comm, ierr)

    end subroutine

    subroutine update_begin_x_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return

        call PetscLogEventBegin(x_comm, ierr)
        call VecGhostUpdateBegin(p_x,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call PetscLogEventBegin(x_comm, ierr)

    end subroutine 

    subroutine update_end_x_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return

        call PetscLogEventBegin(x_comm, ierr)
        call VecGhostUpdateEnd(p_x,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call PetscLogEventEnd(x_comm, ierr)

    end subroutine

    subroutine update_begin_y_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return

        call PetscLogEventBegin(y_comm, ierr)
        call VecGhostUpdateBegin(p_y,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call PetscLogEventBegin(y_comm, ierr)

    end subroutine 

    subroutine update_end_y_ghost()
        implicit none
        PetscErrorCode      :: ierr
        if (proc==1) return

        call PetscLogEventBegin(y_comm, ierr)
        call VecGhostUpdateEnd(p_y,INSERT_VALUES,SCATTER_FORWARD,ierr)
        call PetscLogEventEnd(y_comm, ierr)

    end subroutine

end module petsc_data_structure_mod
