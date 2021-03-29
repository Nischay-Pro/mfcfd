program meshfree_solver
#include <petsc/finclude/petscsys.h>
    
    use petscsys
    use cudafor
    use parameter_mod
    use data_structure_mod
    use petsc_data_structure_mod
    use point_preprocessor_mod
    use initial_conditions_mod
    use q_lskum_mod
    
    implicit none
    integer :: istat, i, j, nDevices=0, local_gpu
    integer :: accessPeer1, accessPeer2
    real*8  :: start,finish, runtime
    type(cudaDeviceProp) :: prop
    real*8  :: totaltime
    MPI_Comm local_comm
    PetscErrorCode  :: ierr
    
    call PetscInitialize('case.in', ierr)
    if(ierr /= 0) stop "Unable to initialize PETSc"
    call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(PETSC_COMM_WORLD, proc, ierr)
    if(rank==0) then
        call execute_command_line('mkdir -p solution')
        call execute_command_line('mkdir -p cp')
    end if
    
    totaltime = MPI_Wtime()
    
    
    if(rank == 0) then
        write(*,*)
        write(*,*)'%%%%%%%%-CUDA Fortran Meshfree Code-%%%%%%%'
        write(*,*)
        write(*,*)'%%%%%%%%%%%%%%%-Device info-%%%%%%%%%%%%%%%'
    end if
    
    istat = cudaGetDeviceCount ( nDevices )
    if(rank == 0) then
        do i = 0, nDevices - 1
            istat = cudaGetDeviceProperties(prop, i)
            write(*,*)'Device Name:               ', trim(prop%name)
            write(*,*)'Compute Capability:        ',prop%major, prop%minor
            write(*,*)'Device number:             ',i
            write(*,*)'MemoryClockRate(KHz):      ',prop%memoryBusWidth
            write(*,*)'PeakMemoryBandwidth(GB/s): ',2.0 *prop%memoryClockRate * &
            & (prop%memoryBusWidth/8) * 1.e-6
            write(*,*)
        end do
        
        if (nDevices .lt. 1) then
            write(*,*) 'ERROR: There are no devices available on this host.  ABORTING.'
        endif
    end if

    call MPI_Comm_split_type(PETSC_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, local_comm, ierr)
    call MPI_Comm_size(local_comm, local_size, ierr)
    call MPI_Comm_rank(local_comm, local_rank, ierr)

    local_gpu = mod(local_rank,local_size)
    WRITE(*,*) local_gpu
    istat = cudaSetDevice(local_gpu)
    if (istat /= 0) then
        print *, 'main: error setting CUDADevice',cudaGetErrorString(istat)
        stop
    endif
    istat = cudaDeviceSetCacheConfig(cudaFuncCachePreferNone)
    if (istat /= 0) then
        print *, 'main: error setting CUDAFuncAttributePreferredSharedMemoryCarveout',cudaGetErrorString(istat)
        stop
    endif
    
    
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
    
    call read_hdf5input_point_data()
    
    ! call read_input_point_data()
    
    !       Allocate solution variables
    
    call allocate_soln()
    call allocate_device_soln()
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
    if(runop == 1)then
        if(rank == 0) then
            write(*,*)'%%%%%%%%%-Using inbuilt solvers-%%%%%%%%%%%'
            write(*,*)
        end if
        call q_lskum()
    end if
    runtime = MPI_Wtime() - runtime
    
    !       Save solution one last time
    
    ! call print_primal_output()
    
    !       destroy petsc vectors and deallocate point/solution vectors
    call dest_petsc()
    call deallocate_soln()
    call deallocate_device_soln()
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
