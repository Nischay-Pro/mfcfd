program meshfree_solver

        use parameter_mod
        use data_structure_mod
        use point_preprocessor_mod
        use q_lskum_mod
        use compute_force_coeffs_mod

        implicit none
        
        write(*,*)
        write(*,*)'%%%%%%%%%%%%%-Serial Meshfree Code-%%%%%%%%%%%'
        write(*,*)

!       Read the case file

        write(*,*)'%%%%%%%%%%%-Reading the case file-%%%%%%%%%'
        write(*,*)
        
        call readnml()

!	Reading the input data ..

        write(*,*)
        write(*,*)'%%%%%%%%%%%%-Reading point data-%%%%%%%%%%%'
        write(*,*)

        call read_input_point_data()
        
        write(*,*)
        write(*,*)'%%%%%%%%%%%%-Reading phi data-%%%%%%%%%%%'
        write(*,*)

        if(read_phi_file == 1) then
            call read_phi_data()
        endif
        
        if(read_phi_file == 0) then
            allocate(point%phi1(4, max_points))
            allocate(point%phi2(4, max_points))
        endif

!       Allocate solution variables

        call allocate_soln()

!       Initialize

        plen = max_points
        write(*,*) 'Number of points:         ', plen
        write(*,*)
        
!	Assign the initial conditions for the primitive variables ..	
        call initial_conditions()
        write(*,*)'%%%%%%%%%%%-Solution initialised-%%%%%%%%%%'
        write(*,*)

!	Primal fixed point iterative solver ..
        call q_lskum()

!       Save solution one last time
        call print_primal_output()

!       deallocate point/solution vectors
        call deallocate_soln()
        call dealloc_points()

        write(*,*)
        write(*,*) '%%%%%%%%%%%-Simulation finished-%%%%%%%%%%%'
        write(*,*) 'Run time:  ',runtime,'seconds'
        write(*,*) 'Total time:',totaltime,'seconds'

end program meshfree_solver
