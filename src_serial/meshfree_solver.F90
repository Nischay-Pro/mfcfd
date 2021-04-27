program meshfree_solver

        use parameter_mod
        use data_structure_mod
        use point_preprocessor_mod
        use q_lskum_mod
        use compute_force_coeffs_mod

        implicit none
        real*8  :: totaltime,runtime


        call execute_command_line('mkdir -p solution')
        call execute_command_line('mkdir -p cp')
        

!       Read the case file
        
        call readnml()

!	Reading the input data ..

        
        call read_hdf5input_point_data()
        
!       Allocate solution variables

        call allocate_soln()


        plen = max_points
        
!	Assign the initial conditions for the primitive variables ..	
        call initial_conditions()

!	Primal fixed point iterative solver ..
  
        call q_lskum()

!       Save solution one last time

        ! call print_primal_output()

        call deallocate_soln()
        call dealloc_points()

end program meshfree_solver
