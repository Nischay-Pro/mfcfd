module point_normals_mod
    
    use data_structure_mod
    
    contains
    
    
    subroutine compute_normals()
        
        implicit none
        
        double precision :: lx, ly, mx, my, rx, ry
        double precision :: nx1, nx2, ny1, ny2, nx, ny
        double precision :: det
        
        integer:: i, j, k, l, m, r
        
        
        !Finding the normals for the points on the shapes ..   
        
        do i = 1, wall_points
            
            m = wall_points_index(i)
            l = point%left(m)
            r = point%right(m)
            
            lx = point%xy(1,l)
            ly = point%xy(2,l)
            
            mx = point%xy(1,m)
            my = point%xy(2,m)
            
            rx = point%xy(1,r)
            ry = point%xy(2,r)
            
            nx1 = my - ly
            nx2 = ry - my
            
            ny1 = mx - lx
            ny2 = rx - mx
            
            nx = 0.5d0*(nx1 + nx2)
            ny = 0.5d0*(ny1 + ny2)
            
            det = dsqrt(nx*nx + ny*ny)
            
            nx = -nx/det
            ny = ny/det
            
            point%nxy(1,m) = nx
            point%nxy(2,m) = ny
            
        enddo
        
        
        !	Finding the normals for the outer boundary points ..
        
        do i = 1, outer_points
            m = outer_points_index(i)
            l = point%left(m)
            r = point%right(m)
            
            lx = point%xy(1,l)
            ly = point%xy(2,l)
            
            mx = point%xy(1,m)
            my = point%xy(2,m)
            
            rx = point%xy(1,r)
            ry = point%xy(2,r)
            
            nx1 = my - ly
            nx2 = ry - my
            
            ny1 = mx - lx
            ny2 = rx - mx
            
            nx = 0.5d0*(nx1 + nx2)
            ny = 0.5d0*(ny1 + ny2)
            
            det = dsqrt(nx*nx + ny*ny)
            
            nx = -nx/det
            ny = ny/det
            
            point%nxy(1,m) = nx
            point%nxy(2,m) = ny
            
        enddo
        
        if(interior_points_normal_flag == 0 .and. format .ne. 2) then
            do i = 1, interior_points
                k = interior_points_index(i)
                point%nxy(1,k) = 0.d0
                point%nxy(2,k) = 1.d0
            enddo
            elseif(interior_points_normal_flag == 1 .and. format .ne. 2) then
                do i = 1, interior_points
                    k = interior_points_index(i)
                    point%nxy(1,k) = 1.d0
                    point%nxy(2,k) = 0.d0
                enddo
            endif
            
        end subroutine 
        
    end module 
