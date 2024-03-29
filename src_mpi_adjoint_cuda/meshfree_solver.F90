program meshfree_solver
#include <petsc/finclude/petscsys.h> 
    
    use petscsys
    use parameter_mod
    use data_structure_mod_diff
    use petsc_data_structure_mod
    use point_preprocessor_mod
    use q_lskum_mod_diff
    use Q_LSKUM_MOD_CHKPTS_DIFF
    use compute_force_coeffs_mod_diff
    use Q_LSKUM_MOD_CHKPTS_CUDA_DIFF
    
    implicit none
    real*8  :: totaltime,runtime
    ! integer :: iterate
    PetscErrorCode  :: ierr
    
    call PetscInitialize('case.in', ierr)
    if(ierr /= 0) stop "Unable to initialize PETSc"
    call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(PETSC_COMM_WORLD, proc, ierr)
    if(rank==0) then
        call execute_command_line('mkdir -p solution')
        call execute_command_line('mkdir -p phistore')
        call execute_command_line('mkdir -p cp')
        call execute_command_line('mkdir -p sensitivity')
    end if
    
    totaltime = MPI_Wtime()
    
    if(rank == 0) then
        write(*,*)
        write(*,*)'%%%%%%%%%%%%%-MPI Meshfree Code-%%%%%%%%%%%'
        write(*,*)
    end if
    
    !       Read the case file
    
    if(rank == 0) then
        write(*,*)'%%%%%%%%%%%-Reading the case file-%%%%%%%%%'
        write(*,*)
    end if
    
    call readcase()
    
    !	Reading the input data ..
    
    if(rank == 0) then
        write(*,*)
        write(*,*)'%%%%%%%%%%%%-Reading HDF5 point data-%%%%%%%%%%%'
        write(*,*)
    end if
    
    call read_input_point_data()

    ! if(rank == 0) then
    !     write(*,*)
    !     write(*,*)'%%%%%%%%%%%%-Reading point data-%%%%%%%%%%%'
    !     write(*,*)
    ! end if

    ! call read_input_point_data()
    
    if(rank == 0) then
        write(*,*)
        write(*,*)'%%%%%%%%%%%%-Reading phi data-%%%%%%%%%%%'
        write(*,*)
    end if
    
    !       Allocate solution variables
    
    call allocate_soln()
    call allocate_soln_b()
    !       Initialize Petsc vectors
    
    if(proc .ne. 1)call init_petsc()
    if(proc == 1) plen = max_points
    if(rank == 0) then
        write(*,*) 'Number of points:         ', plen
        write(*,*)
    end if
    !	Assign the initial conditions for the primitive variables ..	
    call initial_conditions()
    if(rank == 0) then
        write(*,*)'%%%%%%%%%%%-Solution initialised-%%%%%%%%%%'
        write(*,*)
    end if
    
    !	Primal fixed point iterative solver ..
    runtime = MPI_Wtime()
    if(ad_mode == 0) then ! Black Box approach
        if(rank == 0) then
            write(*,*)'%%%%%%%%%-Using Black box approach-%%%%%%%%'
            write(*,*)
        end if
        call q_lskum_b()
    else
        if(rank == 0) then
            write(*,*)'%%%%%%%-Using Checkpointing approach-%%%%%%'
            write(*,*)
        end if
        call Q_LSKUM_CHKPTS_B_CUDA()
    end if
    runtime = MPI_Wtime() - runtime
    
    !       Save solution one last time
    call print_primal_output()
    
    !       destroy petsc vectors and deallocate point/solution vectors
    call dest_petsc()
    call deallocate_soln()
    call deallocate_soln_b()
    call dealloc_points()
    
    totaltime = MPI_Wtime() - totaltime
    
    if(rank == 0) then
        write(*,*)
        write(*,*) '%%%%%%%%%%%-Simulation finished-%%%%%%%%%%%'
        write(*,*) 'Run time:  ',runtime,'seconds'
        write(*,*) 'Total time:',totaltime,'seconds'
    end if
    
    !       stop petsc
    call PetscFinalize(ierr)
    
end program meshfree_solver
