module generate_connectivity_mod
    
    use data_structure_mod
    use petsc_data_structure_mod
    
    contains
    
    subroutine generate_connectivity()
        
        implicit none
        
        integer :: i, k
        real*8 :: nx, ny
        
        do k = 1, interior_points
            i = interior_points_index(k)
            nx = point%nxy(1,i)
            ny = point%nxy(2,i)
            call get_interior_neighbours(i, nx, ny)
            call check_condition_number(i, nx, ny)
        enddo
        
        do k = 1, wall_points
            i = wall_points_index(k)
            nx = point%nxy(1,i)
            ny = point%nxy(2,i)
            call get_wall_boundary_neighbours(i, nx, ny)
        enddo
        
        do k = 1, outer_points
            i = outer_points_index(k)
            nx = point%nxy(1,i)
            ny = point%nxy(2,i)
            call get_outer_boundary_neighbours(i, nx, ny)
        enddo
        
    end subroutine 
    
    subroutine get_interior_neighbours(i, nx, ny)
        
        implicit none
        
        real*8 :: xi, yi, xk, yk
        real*8 :: delx, dely, dels, deln
        real*8 :: nx, ny, tx, ty
        integer :: i, r, count, nbh
        
        xi = point%xy(1,i)
        yi = point%xy(2,i)
        
        tx = ny
        ty = -nx
        
        point%xpos_nbhs(i) = 0
        point%xneg_nbhs(i) = 0
        point%ypos_nbhs(i) = 0
        point%yneg_nbhs(i) = 0
        
        do r=1, point%nbhs(i)
            nbh = point%conn(i,r)
            xk = point%xy(1,nbh)
            yk = point%xy(2,nbh)
            
            nbh = find_loc_f90(point%conn, 20, i, nbh)
            
            delx = xk - xi
            dely = yk - yi
            
            dels = delx*tx + dely*ty
            deln = delx*nx + dely*ny
            
            if(dels .le. 0.0d0) then
                
                point%xpos_nbhs(i) = point%xpos_nbhs(i) + 1;
                
                count = point%xpos_nbhs(i);
                point%xpos_conn(i,count) = nbh;
                
            endif
            
            if(dels .ge. 0.0d0) then
                
                point%xneg_nbhs(i) = point%xneg_nbhs(i) + 1;
                
                count = point%xneg_nbhs(i);
                point%xneg_conn(i,count) = nbh;
                
            endif
            
            if(deln .le. 0.0d0) then
                
                point%ypos_nbhs(i) = point%ypos_nbhs(i) + 1;
                
                count = point%ypos_nbhs(i);
                point%ypos_conn(i,count) = nbh;
                
            endif
            
            if(deln .ge. 0.0d0) then
                
                point%yneg_nbhs(i) = point%yneg_nbhs(i) + 1;
                
                count = point%yneg_nbhs(i);
                point%yneg_conn(i,count) = nbh;
                
            endif
            
        enddo
        
        if(point%xpos_nbhs(i) == 0) then
            print*,"WARNING!!! xpos zero for interior point number:", i,", rank:", rank
            elseif(point%xneg_nbhs(i) == 0) then
                print*,"WARNING!!! xneg zero for interior point number:", i,", rank:", rank
                elseif(point%ypos_nbhs(i) == 0) then
                    print*,"WARNING!!! ypos zero for interior point number:", i,", rank:", rank
                    elseif(point%yneg_nbhs(i) == 0) then
                        print*,"WARNING!!! yneg zero for interior point number:", i,", rank:", rank
                    end if
                    
                end subroutine
                
                subroutine get_wall_boundary_neighbours(i, nx, ny)
                    
                    implicit none
                    
                    real*8 :: xi, yi, xk, yk
                    real*8 :: delx, dely, dels, deln
                    real*8 :: nx, ny, tx, ty
                    integer :: i, r, count, nbh
                    
                    
                    xi = point%xy(1,i)
                    yi = point%xy(2,i)
                    
                    tx = ny
                    ty = -nx
                    
                    point%xpos_nbhs(i) = 0
                    point%xneg_nbhs(i) = 0
                    point%yneg_nbhs(i) = 0
                    
                    do r=1, point%nbhs(i)
                        
                        nbh = point%conn(i,r)
                        
                        xk = point%xy(1,nbh)
                        yk = point%xy(2,nbh)
                        
                        nbh = find_loc_f90(point%conn, 20, i, nbh)
                        
                        delx = xk - xi
                        dely = yk - yi
                        
                        dels = delx*tx + dely*ty
                        deln = delx*nx + dely*ny
                        
                        if(dels .le. 0.0d0) then
                            
                            point%xpos_nbhs(i) = point%xpos_nbhs(i) + 1;
                            
                            count = point%xpos_nbhs(i);
                            point%xpos_conn(i,count) = nbh;
                            
                        endif
                        
                        if(dels .ge. 0.0d0) then
                            
                            point%xneg_nbhs(i) = point%xneg_nbhs(i) + 1;
                            
                            count = point%xneg_nbhs(i);
                            point%xneg_conn(i,count) = nbh;
                            
                        endif
                        
                        point%yneg_nbhs(i) = point%yneg_nbhs(i) + 1;
                        
                        count = point%yneg_nbhs(i);
                        point%yneg_conn(i,count) = nbh;
                        
                    enddo
                    
                    if(point%xpos_nbhs(i) == 0) then
                        print*,"WARNING!!! xpos zero for wall point number:", i,", rank:", rank
                        elseif(point%xneg_nbhs(i) == 0) then
                            print*,"WARNING!!! xneg zero for wall point number:", i,", rank:", rank
                            elseif(point%yneg_nbhs(i) == 0) then
                                print*,"WARNING!!! yneg zero for wall point number:", i,", rank:", rank
                            end if
                            
                        end subroutine
                        
                        
                        subroutine get_outer_boundary_neighbours(i, nx, ny)
                            
                            implicit none
                            
                            real*8 :: xi, yi, xk, yk
                            real*8 :: delx, dely, dels, deln
                            real*8 :: nx, ny, tx, ty
                            integer :: i, r, count, nbh
                            
                            
                            xi = point%xy(1,i)
                            yi = point%xy(2,i)
                            
                            tx = ny
                            ty = -nx
                            
                            point%xpos_nbhs(i) = 0
                            point%xneg_nbhs(i) = 0
                            point%ypos_nbhs(i) = 0
                            
                            do r=1, point%nbhs(i)
                                
                                nbh = point%conn(i,r)
                                
                                xk = point%xy(1,nbh)
                                yk = point%xy(2,nbh)
                                
                                nbh = find_loc_f90(point%conn, 20, i, nbh)
                                
                                delx = xk - xi
                                dely = yk - yi
                                
                                dels = delx*tx + dely*ty
                                deln = delx*nx + dely*ny
                                
                                if(dels .le. 0.0d0) then
                                    
                                    point%xpos_nbhs(i) = point%xpos_nbhs(i) + 1;
                                    
                                    count = point%xpos_nbhs(i);
                                    point%xpos_conn(i,count) = nbh;
                                    
                                endif
                                
                                if(dels .ge. 0.0d0) then
                                    
                                    point%xneg_nbhs(i) = point%xneg_nbhs(i) + 1;
                                    
                                    count = point%xneg_nbhs(i);
                                    point%xneg_conn(i,count) = nbh;
                                    
                                endif
                                
                                
                                point%ypos_nbhs(i) = point%ypos_nbhs(i) + 1;
                                
                                count = point%ypos_nbhs(i);
                                point%ypos_conn(i,count) = nbh;
                                
                            enddo
                            
                            if(point%xpos_nbhs(i) == 0) then
                                print*,"WARNING!!! xpos zero for outer point number:", i,", rank:", rank
                                elseif(point%xneg_nbhs(i) == 0) then
                                    print*,"WARNING!!! xneg zero for outer point number:", i,", rank:", rank
                                    elseif(point%ypos_nbhs(i) == 0) then
                                        print*,"WARNING!!! ypos zero for outer point number:", i,", rank:", rank
                                    end if
                                    
                                end subroutine
                                
                                subroutine check_condition_number(i, nx, ny)
                                    
                                    ! Use lapack_example_aux, Only: nagf_file_print_matrix_real_gen
                                    ! Use lapack_interfaces, Only: dbdsqr, dgebrd, dlacpy, dorgbr
                                    ! Use lapack_precision, Only: dp
                                    
                                    implicit none
                                    
                                    integer :: i
                                    real*8 :: nx, ny
                                    
                                end subroutine
                                
                                integer function find_loc_f90(array, pointcount, pidx, nbhvalue)
                                integer, dimension(:,:) :: array
                                integer :: pointcount, pidx, nbhvalue, i
                                do i = 1, pointcount
                                    if (array(pidx, i) == nbhvalue) then
                                        find_loc_f90 = i
                                        return
                                    endif
                                end do
                                print*, "warning could not find point in conn", pidx, nbhvalue
                                stop
                            end function
                            
                        end module
                        