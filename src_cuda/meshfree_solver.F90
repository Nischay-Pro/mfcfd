program meshfree_solver

        use cudafor
        use parameter_mod
        use data_structure_mod
        use point_preprocessor_mod
        use initial_conditions_mod
        use q_lskum_mod
        use post_processing_mod


        implicit none
        integer :: istat
        real*8  :: start,finish
        real*8  :: startr,finishr
        type(cudaDeviceProp) :: prop
        
        call cpu_time(start)

        write(*,*)
        write(*,*)'%%%%%%%%-CUDA Fortran Meshfree Code-%%%%%%%'
        write(*,*)
        write(*,*)'%%%%%%%%%%%%%%%-Device info-%%%%%%%%%%%%%%%'
        istat = cudaGetDeviceProperties(prop, 0)
        write(*,*)'Device Name:', trim(prop%name)
        write(*,*)'Compute Capability: ',prop%major, prop%minor
        write(*,*)

!       Set up case input

        call readnml()

!	Reading the input data ..

        write(*,*) '%%%%%%%%%%%%-Reading point file-%%%%%%%%%%%%'
        call read_input_point_data()
        write(*,*) 'Number of points:         ', max_points
        write(*,*) 'Number of wall points:    ', wall_points
        write(*,*) 'Number of interior points:', interior_points
        write(*,*) 'Number of outer points:   ', outer_points
        write(*,*)

!       Allocate solution variables

        call allocate_soln()

!       Allocate device solution variables
        
        call allocate_device_soln()

!	Assign the initial conditions for the primitive variables ..	

        call initial_conditions()
        write(*,*) '%%%%%%%%%%%-Solution initialized-%%%%%%%%%%%'
        write(*,*)
       
!	Primal fixed point iterative solver ..
        
        call cpu_time(startr) 
        call q_lskum()
        call cpu_time(finishr) 

!       Save solution one last time
        call print_primal_output()


!       Deallocate point/solution vectors
        call deallocate_soln()
        call dealloc_points()
        call deallocate_device_soln()

        call cpu_time(finish) 
        print*,'Simulation completed!!!'
        print*,'runtime:',finishr-startr
        print*,'totaltime:',finish-start
        
end program meshfree_solver