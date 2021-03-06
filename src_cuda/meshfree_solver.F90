program meshfree_solver
    
    use cudafor
    use parameter_mod
    use data_structure_mod
    use point_preprocessor_mod
    use initial_conditions_mod
    use q_lskum_mod
    use post_processing_mod
    use objective_function_mod
    use adaptation_sensors_mod
    
    
    implicit none
    integer :: istat, i, j, nDevices=0
    integer :: accessPeer1, accessPeer2
    real*8  :: start,finish, runtime
    type(cudaDeviceProp) :: prop
    
    call cpu_time(start)
    
    write(*,*)
    write(*,*)'%%%%%%%%-CUDA Fortran Meshfree Code-%%%%%%%'
    write(*,*)
    write(*,*)'%%%%%%%%%%%%%%%-Device info-%%%%%%%%%%%%%%%'
    istat = cudaGetDeviceCount ( nDevices )
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
    
    istat = cudaDeviceSetCacheConfig(cudaFuncCachePreferShared)
    if (istat /= 0) then
       print *, 'main: error setting cudaFuncAttributePreferredSharedMemoryCarveout',cudaGetErrorString(istat)
       stop
    endif
    
    ! establish the Peer-to-Peer connections
    ! do i=1, istat-1
    !     do j=i+1,istat
    !         istat = cudaDeviceCanAccessPeer(accessPeer1, devices(i), devices(j))
    !         istat = cudaDeviceCanAccessPeer(accessPeer2, devices(j), devices(i))
    !         if (accessPeer1 .eq. 1 .and. accessPeer2 .eq. 1 ) then
    !             write(*,"('Peer-to-Peer capable between ',i0, ' and ', i0)") &
    !             devices(i), devices(j)
    !             ! enable peer access required for copies
    !             istat = cudaSetDevice(devices(i))  
    !             istat = cudaDeviceEnablePeerAccess(devices(j), 0) 
    !             istat = cudaSetDevice(devices(j))
    !             istat = cudaDeviceEnablePeerAccess(devices(i), 0)  
    !         else
    !             write(*,"('Peer-to-Peer not capable between ',i0, ' and ', i0)") &
    !             devices(i), devices(j)
    !             stop
    !         endif
    !     enddo
    ! enddo
    
    !       Set up case input
    call readnml()

    	! Reading the input data ..
    write(*,*)'%%%%%%%%%%%%-Reading HDF5 point data-%%%%%%%%%%%'
    call read_hdf5input_point_data()

    ! write(*,*) '%%%%%%%%%%%%-Reading point file-%%%%%%%%%%%'
    ! call read_input_point_data()
    
    write(*,*) 'Number of points:         ', max_points
    write(*,*) 'Number of wall points:    ', wall_points
    write(*,*) 'Number of shape points:   ', shape_points
    write(*,*) 'Number of interior points:', interior_points
    write(*,*) 'Number of outer points:   ', outer_points
    write(*,*)
    
    !       Allocate solution variables
    call allocate_soln()
    
    !       Allocate device solution variables
    call allocate_device_soln()
    
    !	Assign the initial conditions for the primitive variables ..	
    call initial_conditions()
    write(*,*) '%%%%%%%%%%%-Solution initialized-%%%%%%%%%%'
    write(*,*)
    
    !	Primal fixed point iterative solver ..
    call q_lskum(runtime)
    
    !       Compute sensor values
    call compute_adapt_sensor()
    
    !       Objective function computation
    call objective_function()
    
    !       Save solution one last time
    call print_primal_output()
    
    !       Deallocate point/solution vectors
    call deallocate_soln()
    call dealloc_points()
    call deallocate_device_soln()
    
    call cpu_time(finish) 
    write(*,*) '%%%%%%%%%%%-Simulation finished-%%%%%%%%%%%'
    write(*,*) 'Run time:  ',runtime,'seconds'
    write(*,*) 'Total time:',finish-start,'seconds'
    
end program meshfree_solver