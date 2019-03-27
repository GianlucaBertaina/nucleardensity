   MODULE manage_cubes


   implicit none

   save

   REAL*8  :: cubelmax, dx, dy, dz
 
   contains


   SUBROUTINE set_cube


        use input_def,    only: x_eq_cart, nat
        use input_def,    only: nxpoints, nypoints, nzpoints   

        implicit none
       

        REAL*8 :: barx, bary, barz, maxdist
        REAL*8, DIMENSION(nat) :: dist


        INTEGER :: i, j, k
        REAL*8, PARAMETER :: maxdist_plus=2.d0


!       Get the geometric center of mass
        barx=0
        bary=0
        barz=0
        DO i=1, nat         
           barx = barx + x_eq_cart(3*i - 2) !sum over x coord
           bary = bary + x_eq_cart(3*i - 1) !sum over y coord
           barz = barz + x_eq_cart(3*i    ) !sum over z coord
        END DO
        barx = barx / nat
        bary = bary / nat
        barz = barz / nat 

!       Centre eq. geom at the geometric ceter of mass 
        DO i=1, nat         
           x_eq_cart(3*i - 2) = x_eq_cart(3*i - 2) - barx !x coord
           x_eq_cart(3*i - 1) = x_eq_cart(3*i - 1) - bary !y coord
           x_eq_cart(3*i    ) = x_eq_cart(3*i    ) - barz !z coord
        END DO

!       Compute distances from centre
        DO i=1, nat
           dist(i) = SQRT(      x_eq_cart(3*i - 2)**2 +&
                                x_eq_cart(3*i - 1)**2 +& 
                                x_eq_cart(3*i    )**2        )
        END DO

        maxdist = MAXVAL(dist)
        cubelmax = maxdist + maxdist_plus

!       Set 3d grid (each point is the midpoint of the voxel edge)
!       Boundaries of the box are given by cubelmax:
!       The grid starts at -cubelmax and ends at cubelmax along x,y,z        
!       the central voxel is centered in the geometric center of mass
        !dx = cubelmax*2 / REAL(nxpoints-1, 8)
        !dy = cubelmax*2 / REAL(nypoints-1, 8)
        !dz = cubelmax*2 / REAL(nzpoints-1, 8)

        dx = cubelmax*2 / REAL(nxpoints, 8)
        dy = cubelmax*2 / REAL(nypoints, 8)
        dz = cubelmax*2 / REAL(nzpoints, 8)

        !print*, 'Voxels edges lenght (a.u.)'
        !print*, 'dx =', dx
        !print*, 'dy =', dy
        !print*, 'dz =', dz

   END SUBROUTINE set_cube



   SUBROUTINE find_cube_index(Ri, ix, iy, iz)

        use input_def,    only: nxpoints, nypoints, nzpoints   
        
        IMPLICIT NONE 
        INTEGER, INTENT(OUT) :: ix, iy, iz
        REAL*8, INTENT(IN) :: Ri(3)
         
           !along x
           !IF (Ri(1) >= 0) THEN
           !   ix = CEILING(Ri(1)/dx)
           !   ix = ( nxpoints / 2 ) + ix
           !ELSE
           !   ix = FLOOR(Ri(1)/dx)
           !   ix = ( nxpoints / 2 ) + ix + 1
           !END IF

           ix = NINT(Ri(1) / dx)
           ix = (nxpoints + 1) / 2 + ix
           IF (ix > nxpoints) ix = ix - 1
           IF (ix < 1) ix = ix + 1

           !along y
           !IF (Ri(2) >= 0) THEN
           !   iy = CEILING(Ri(2)/dy)
           !   iy = ( nypoints / 2 ) + iy
           !ELSE
           !   iy = FLOOR(Ri(2)/dy)
           !   iy = ( nypoints / 2 ) + iy + 1
           !END IF

           iy = NINT(Ri(2) / dy)
           iy = (nypoints + 1) / 2 + iy
           IF (iy > nypoints) iy = iy - 1
           IF (iy < 1) iy = iy + 1

           !along z
           !IF (Ri(3)>= 0) THEN
           !   iz = CEILING(Ri(3)/dz)
           !   iz = ( nzpoints / 2 ) + iz 
           !ELSE
           !   iz = FLOOR(Ri(3)/dz)
           !   iz = ( nzpoints / 2 ) + iz + 1
           !END IF

           iz = NINT(Ri(3) / dz)
           iz = (nzpoints + 1) / 2 + iz
           IF (iz > nzpoints) iz = iz - 1
           IF (iz < 1) iz = iz + 1

   END SUBROUTINE find_cube_index



   SUBROUTINE print_cube(density)

        use input_def,    only: x_eq_cart, nat, Z_nuclei
        use input_def,    only: nxpoints, nypoints, nzpoints   
        use io_units_def        
        use constants
 
        implicit none

        REAL*8, intent(in) :: density(nxpoints, nypoints, nzpoints, nat)
        character(25) :: filename
        integer :: i, j, k, iat

      DO iat=1, nat

        !Preparing filename string
        WRITE(filename, "('nucl_',I0,'.cube')") iat
        
        !associate a unit to the cube file
        OPEN (UNIT=unit_cube, FILE=trim(filename), STATUS='replace')

        ! write out cube file introduction
        WRITE(unit_cube,*) "!Cube format assumes densities are printed at vertices of cubes"
        WRITE(unit_cube,*) "!voxels index order: Z, Y, X"
        ! Shifts origin of box (vertex of box) so that it corresponds to the center of 
        ! the first cube in this program convention
        ! TODO: do not shift if the main routine is modified to evaluate densities at vertices
        WRITE(unit_cube,*) nat, -cubelmax + dx*0.5d0 ,&
                                -cubelmax + dy*0.5d0 ,&
                                -cubelmax + dz*0.5d0 
        WRITE(unit_cube,*) nxpoints, dx, zero_real, zero_real
        WRITE(unit_cube,*) nypoints, zero_real, dy, zero_real
        WRITE(unit_cube,*) nzpoints, zero_real, zero_real, dz


!        Print out the equilibrium geometry of the molecule
         DO i=1, nat
            WRITE(unit_cube,*) Z_nuclei(i), real(Z_nuclei(i),8), &
                        x_eq_cart(3*i-2),&
                        x_eq_cart(3*i-1),&
                        x_eq_cart(3*i)  
         END DO

         DO i=1, nxpoints
           DO j=1, nypoints
                 WRITE(unit_cube, '( 6F12.6 )') (density(i, j, k, iat), k=1, nzpoints)
           END DO
         END DO

         CLOSE(unit_cube)

      END DO

   END SUBROUTINE print_cube


   SUBROUTINE print_cube_errorbar(density)

        use input_def,    only: x_eq_cart, nat, Z_nuclei
        use input_def,    only: nxpoints, nypoints, nzpoints   
        use io_units_def        
        use constants
 
        implicit none

        REAL*8, intent(in) :: density(nxpoints, nypoints, nzpoints, nat)
        character(25) :: filename
        integer :: i, j, k, iat

      DO iat=1, nat

        !Preparing filename string
        WRITE(filename, "('err_nucl_',I0,'.cube')") iat
        
        !associate a unit to the cube file
        OPEN (UNIT=unit_cube, FILE=trim(filename), STATUS='replace')

        ! write out cube file "introdction"
        ! Wriete somethig useful in comment lines TBD
        WRITE(unit_cube,*) "!comment"
        WRITE(unit_cube,*) "!comment"
        WRITE(unit_cube,*) nat, -cubelmax,&
                                -cubelmax,&
                                -cubelmax
        WRITE(unit_cube,*) nxpoints, dx, zero_real, zero_real
        WRITE(unit_cube,*) nypoints, zero_real, dy, zero_real
        WRITE(unit_cube,*) nzpoints, zero_real, zero_real, dz


!        Print out the equilibrium geometry of the molecule
         DO i=1, nat
            WRITE(unit_cube,*) Z_nuclei(i), real(Z_nuclei(i),8), &
                        x_eq_cart(3*i-2),&
                        x_eq_cart(3*i-1),&
                        x_eq_cart(3*i)  
         END DO

         DO i=1, nxpoints
           DO j=1, nypoints
                 WRITE(unit_cube, '( 6F12.6 )') (density(i, j, k, iat), k=1, nzpoints)
           END DO
         END DO

         CLOSE(unit_cube)

      END DO

   END SUBROUTINE print_cube_errorbar




   END MODULE manage_cubes
