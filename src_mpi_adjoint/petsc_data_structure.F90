module petsc_data_structure_mod
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petsclog.h>


        use petscvec
        use parameter_mod
        use data_structure_mod_diff

        implicit none

        PetscMPIInt          :: rank,proc
        Vec                  :: p_dq
        Vec                  :: p_ddq
        Vec                  :: p_qm
        Vec                  :: p_prim

                ! Adjoint PETSc variables
        Vec                  :: pb_dq, pb_q
        Vec                  :: pb_prim, pb_x, pb_y, pb_qm
        PetscLogEvent        :: dq_comm, prim_comm, qm_comm, q_comm, x_comm, y_comm

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
                        

                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,2*4,2*4*local_points,&
                        &PETSC_DECIDE,ghost_points,pghost,pointb%dq(1,1,1),pb_dq,ierr)
                
                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,2*4,2*4*local_points,&
                        &PETSC_DECIDE,ghost_points,pghost,pointb%qm(1,1,1),pb_qm,ierr)

                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,4,4*local_points,&
                        &PETSC_DECIDE,ghost_points,pghost,pointb%prim(1,1),pb_prim,ierr)
                

                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,4,4*local_points,&
                        &PETSC_DECIDE,ghost_points,pghost,pointb%q(1,1),pb_q,ierr)

                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,1,local_points,&
                        &PETSC_DECIDE,ghost_points,pghost,pointb%x(1),pb_x,ierr)

                call VecCreateGhostBlockWithArray(PETSC_COMM_WORLD,1,local_points,&
                        &PETSC_DECIDE,ghost_points,pghost,pointb%y(1),pb_y,ierr)


                call VecGetSize(p_prim,plen,ierr)
                plen = plen/4

                call PetscLogEventRegister('dq_comm',  0,dq_comm,ierr);
                call PetscLogEventRegister('qm_comm',  0,qm_comm,ierr);
                call PetscLogEventRegister('prim_comm',  0,prim_comm,ierr);
                call PetscLogEventRegister('q_comm',  0,q_comm,ierr);
                call PetscLogEventRegister('x_comm',  0,x_comm,ierr);
                call PetscLogEventRegister('y_comm',  0,y_comm,ierr);
        end subroutine 

        subroutine dest_petsc()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return
                    

                call VecDestroy(p_dq,ierr)
                call VecDestroy(p_qm,ierr)
                call VecDestroy(p_prim,ierr)
                call VecDestroy(pb_dq,ierr)
                call VecDestroy(pb_qm,ierr)
                call VecDestroy(pb_prim,ierr)
                call VecDestroy(pb_x,ierr)
                call VecDestroy(pb_y,ierr)
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


        !Adjoint Begin Add

        subroutine update_begin_dqb_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateBegin(pb_dq,ADD_VALUES,SCATTER_REVERSE,ierr)
        end subroutine 

        subroutine update_begin_qmb_ghost()
            implicit none
            PetscErrorCode      :: ierr
            if (proc==1) return

            call VecGhostUpdateBegin(pb_qm,ADD_VALUES,SCATTER_REVERSE,ierr)
        end subroutine         
        

        subroutine update_begin_primb_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return
 
                call VecGhostUpdateBegin(pb_prim,ADD_VALUES,SCATTER_REVERSE,ierr)
        end subroutine 

        subroutine update_begin_qb_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return
 
                call VecGhostUpdateBegin(pb_q,ADD_VALUES,SCATTER_REVERSE,ierr)
        end subroutine 

        subroutine update_begin_xb_ghost()
            implicit none
            PetscErrorCode      :: ierr
            if (proc==1) return

            call VecGhostUpdateBegin(pb_x,ADD_VALUES,SCATTER_REVERSE,ierr)
        end subroutine 

        subroutine update_begin_yb_ghost()
            implicit none
            PetscErrorCode      :: ierr
            if (proc==1) return

            call VecGhostUpdateBegin(pb_y,ADD_VALUES,SCATTER_REVERSE,ierr)
        end subroutine 


        !Adjoint Begin Insert

        subroutine update_begin_dqb_ghost_()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call VecGhostUpdateBegin(pb_dq,INSERT_VALUES,SCATTER_FORWARD,ierr)
        end subroutine 
        

        subroutine update_begin_primb_ghost_()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return
 
                call VecGhostUpdateBegin(pb_prim,INSERT_VALUES,SCATTER_FORWARD,ierr)
        end subroutine 

        subroutine update_begin_qb_ghost_()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return
 
                call VecGhostUpdateBegin(pb_q,INSERT_VALUES,SCATTER_FORWARD,ierr)
        end subroutine 

        ! Forward End Insert

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

        ! Adjoint End Add
        subroutine update_end_dqb_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call PetscLogEventBegin(dq_comm, ierr)
                call VecGhostUpdateEnd(pb_dq,ADD_VALUES,SCATTER_REVERSE,ierr)
                call PetscLogEventEnd(dq_comm, ierr)

        end subroutine 

        subroutine update_end_qmb_ghost()
            implicit none
            PetscErrorCode      :: ierr
            if (proc==1) return

            call PetscLogEventBegin(qm_comm, ierr)
            call VecGhostUpdateEnd(pb_qm,ADD_VALUES,SCATTER_REVERSE,ierr)
            call PetscLogEventEnd(qm_comm, ierr)

        end subroutine 


        subroutine update_end_primb_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call PetscLogEventBegin(prim_comm, ierr)
                call VecGhostUpdateEnd(pb_prim,ADD_VALUES,SCATTER_REVERSE,ierr)
                call PetscLogEventEnd(prim_comm, ierr)
        end subroutine
        
        subroutine update_end_qb_ghost()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call PetscLogEventBegin(q_comm, ierr)
                call VecGhostUpdateEnd(pb_q,ADD_VALUES,SCATTER_REVERSE,ierr)
                call PetscLogEventBegin(q_comm, ierr)

        end subroutine

        subroutine update_end_xb_ghost()
            implicit none
            PetscErrorCode      :: ierr
            if (proc==1) return

            call PetscLogEventBegin(x_comm, ierr)
            call VecGhostUpdateEnd(pb_x,ADD_VALUES,SCATTER_REVERSE,ierr)
            call PetscLogEventBegin(x_comm, ierr)

        end subroutine

        subroutine update_end_yb_ghost()
            implicit none
            PetscErrorCode      :: ierr
            if (proc==1) return

            call PetscLogEventBegin(y_comm, ierr)
            call VecGhostUpdateEnd(pb_y,ADD_VALUES,SCATTER_REVERSE,ierr)
            call PetscLogEventBegin(y_comm, ierr)

        end subroutine


        ! Adjoint End Insert

        subroutine update_end_dqb_ghost_()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call PetscLogEventBegin(dq_comm, ierr)
                call VecGhostUpdateEnd(pb_dq,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call PetscLogEventEnd(prim_comm, ierr)

        end subroutine 


        subroutine update_end_primb_ghost_()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call PetscLogEventBegin(prim_comm, ierr)
                call VecGhostUpdateEnd(pb_prim,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call PetscLogEventEnd(prim_comm, ierr)
        end subroutine
        
        subroutine update_end_qb_ghost_()
                implicit none
                PetscErrorCode      :: ierr
                if (proc==1) return

                call PetscLogEventBegin(q_comm, ierr)
                call VecGhostUpdateEnd(pb_q,INSERT_VALUES,SCATTER_FORWARD,ierr)
                call PetscLogEventBegin(q_comm, ierr)

        end subroutine


        subroutine UPDATE_DQB()
                CALL UPDATE_BEGIN_DQB_GHOST()
                CALL UPDATE_END_DQB_GHOST()
                CALL UPDATE_BEGIN_DQB_GHOST_()
                CALL UPDATE_END_DQB_GHOST_()
        end subroutine

        subroutine UPDATE_QB()
                CALL UPDATE_BEGIN_QB_GHOST()
                CALL UPDATE_END_QB_GHOST()
                CALL UPDATE_BEGIN_QB_GHOST_()
                CALL UPDATE_END_QB_GHOST_()
        end subroutine

        subroutine UPDATE_PRIMB()
                CALL UPDATE_BEGIN_PRIMB_GHOST()
                CALL UPDATE_END_PRIMB_GHOST()
                CALL UPDATE_BEGIN_PRIMB_GHOST_()
                CALL UPDATE_END_PRIMB_GHOST_()
        end subroutine

end module petsc_data_structure_mod
