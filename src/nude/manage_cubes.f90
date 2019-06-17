MODULE manage_cubes
   use input_def
   implicit none

   save

   REAL*8              :: cubelmax, dx, dy, dz,voxel
   REAL*8, ALLOCATABLE :: density(:,:,:,:) !, density_sq(:,:,:,:) , density_errorbar(:,:,:,:)
   REAL*8, ALLOCATABLE :: density_red(:,:,:,:)  !, density_sq_red(:,:,:,:)
   INTEGER             :: density_elem
 
   contains

   SUBROUTINE set_cube

        use input_def,    only: x_eq_cart, nat
        use input_def,    only: nxpoints, nypoints, nzpoints   

        implicit none
       

        REAL*8 :: barx, bary, barz, maxdist
        REAL*8, DIMENSION(nat) :: dist


        INTEGER :: i
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

        voxel=dx*dy*dz

        !print*, 'Voxels edges lenght (a.u.)'
        !print*, 'dx =', dx
        !print*, 'dy =', dy
        !print*, 'dz =', dz

   END SUBROUTINE set_cube

   SUBROUTINE allocate_densities()
      ! Initialize density (a 3D array for each atom)
      allocate(density(nxpoints, nypoints, nzpoints, nat))
      !allocate(density_sq(nxpoints, nypoints, nzpoints, nat))
      allocate(density_red(nxpoints, nypoints, nzpoints, nat))
      !allocate(density_sq_red(nxpoints, nypoints, nzpoints, nat))
      !allocate(density_errorbar(nxpoints, nypoints, nzpoints, nat))
      density_elem = nxpoints*nypoints*nzpoints*nat
      density           = 0.d0
      density_red       = 0.d0
      !density_sq       = 0.d0
      !density_sq_red   = 0.d0
      !density_errorbar = 0.d0
   END SUBROUTINE


    SUBROUTINE update_densities(xx,bar_wfn_sq)
      real*8, intent(in) :: xx(1:ncart)
      real*8, intent(in) :: bar_wfn_sq
      integer            :: i,ix,iy,iz
      real*8             :: Ri(3)

      do i = 1,nat
        !
        ! coordinates in au of i^th atom
        Ri(:) = xx(3*i-2:3*i)
        !
        ! Corresponding grid cube indexes
        call find_cube_index(Ri,ix,iy,iz)
        !
        ! Update density
        density(ix,iy,iz,i) = density(ix,iy,iz,i) + bar_wfn_sq
        !density_sq(ix,iy,iz,i) = density_sq(ix,iy,iz,i) + bar_wfn_sq**2
        !
      enddo

    END SUBROUTINE


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


    SUBROUTINE MPI_REDUCE_DENSITIES()

      use mpimod
      implicit none
      integer :: err_mpi

      CALL MPI_REDUCE(density, density_red, density_elem, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD, err_mpi)
      !CALL MPI_REDUCE(density_sq, density_sq_red, density_elem, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD, err_mpi)

    END SUBROUTINE


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

          WRITE(unit_cube,"(1I5,3(1ES22.15,1x))") nat, -cubelmax + dx*0.5d0 ,&
                                                       -cubelmax + dy*0.5d0 ,&
                                                       -cubelmax + dz*0.5d0
          WRITE(unit_cube,"(1I6,3(1ES22.15,1x))") nxpoints, dx, zero_real, zero_real
          WRITE(unit_cube,"(1I6,3(1ES22.15,1x))") nypoints, zero_real, dy, zero_real
          WRITE(unit_cube,"(1I6,3(1ES22.15,1x))") nzpoints, zero_real, zero_real, dz


  !        Print out the equilibrium geometry of the molecule
           DO i=1, nat
              WRITE(unit_cube,"(1I4,4(1ES22.15,1x))") Z_nuclei(i), real(Z_nuclei(i),8), &
                                                                   x_eq_cart(3*i-2:3*i)
           END DO

           DO i=1, nxpoints
             DO j=1, nypoints
                   WRITE(unit_cube, '( 6F12.6 )') (density(i, j, k, iat), k=1, nzpoints)
             END DO
           END DO

           CLOSE(unit_cube)

        END DO

   END SUBROUTINE print_cube


    SUBROUTINE print_normalized_densities(Nsteps_MC_tot,tot_int_red)
      integer, intent(in) :: Nsteps_MC_tot
      real*8,  intent(in) :: tot_int_red

      integer :: ii
      !
      ! NO
      !density_red = density_red / Nsteps_MC_tot
      !density_sq_red = density_sq_red / Nsteps_MC_tot
      !
      ! NO
      !density = density * (norm_gaussian_normal / Nsteps_MC)
      !
      ! SI'
      density_red = density_red / tot_int_red
      !density_sq_red = density_sq_red / tot_int_red
      !density_errorbar = SQRT( (density_sq_red - density_red**2) / Nsteps_MC_tot )
      !
      !Normalize by the voxel dimension
      density_red = density_red / voxel
      !density_errorbar = density_errorbar / voxel
      !
      !Check the normalization of each nucleus density
      DO ii=1, nat
        print*, 'Check norm of density of atom ',ii,': ', sum(density_red(:,:,:,ii))*voxel
      END DO
      !
      ! Print out nuclear densities on output files
      print*, 'Printing cube files'
      call print_cube(density_red)
      !
      ! print*, 'printing cube files with density error bar'
      ! call print_cube_errorbar(density_errorbar)

    END SUBROUTINE


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

        WRITE(unit_cube,"(1I5,3(1ES22.15,1x))") nat, -cubelmax + dx*0.5d0 ,&
                                                     -cubelmax + dy*0.5d0 ,&
                                                     -cubelmax + dz*0.5d0
        WRITE(unit_cube,"(1I6,3(1ES22.15,1x))") nxpoints, dx, zero_real, zero_real
        WRITE(unit_cube,"(1I6,3(1ES22.15,1x))") nypoints, zero_real, dy, zero_real
        WRITE(unit_cube,"(1I6,3(1ES22.15,1x))") nzpoints, zero_real, zero_real, dz

!        Print out the equilibrium geometry of the molecule
         DO i=1, nat
            WRITE(unit_cube,"(1I4,4(1ES22.15,1x))") Z_nuclei(i), real(Z_nuclei(i),8), &
                                                                 x_eq_cart(3*i-2:3*i)
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


MODULE manage_bonds

    use input_def, only : x_eq_cart, nat, ncart,bondpoints  , nbonds,bond_pair,bond_name
    implicit none
    save

    REAL*8, ALLOCATABLE :: bond_density(:,:), bond_density2(:,:),bond_density_errorbar(:,:)
    REAL*8, ALLOCATABLE :: bonddr(:),bondvol(:)

    INTEGER             :: bond_density_elem
    REAL*8, ALLOCATABLE :: bond_density_red(:,:),bond_density2_red(:,:)

    contains

    SUBROUTINE set_bonds
      use constants
      implicit none
      REAL*8 :: dist

      INTEGER :: i,j,b
      REAL*8, PARAMETER :: maxdist_plus=4.d0 ! Bohr

      ALLOCATE(bonddr(nbonds),bondvol(nbonds))

      !       Compute distances
      DO b=1, nbonds

        i = bond_pair(1,b)
        j = bond_pair(2,b)

        dist = SQRT( (x_eq_cart(3*i - 2)-x_eq_cart(3*j - 2))**2 +&
                     (x_eq_cart(3*i - 1)-x_eq_cart(3*j - 1))**2 +&
                     (x_eq_cart(3*i    )-x_eq_cart(3*j    ))**2  )

        bonddr(b)   = (dist + maxdist_plus) / bondpoints    ! Bohr
        bondvol(b)  = bonddr(b) ! SICCOME NON ABBIAMO ROTAZIONI
        !bondvol(b)  = 4*pi*bonddr(b)**2  ! DA VERIFICARE, SICCOME NON ABBIAMO ROTAZIONI
      END DO

    END SUBROUTINE set_bonds


    SUBROUTINE allocate_bonds()

      allocate(bond_density_red(1:bondpoints,1:nbonds))
      allocate(bond_density2_red(1:bondpoints,1:nbonds))
      allocate(bond_density(1:bondpoints,1:nbonds))
      allocate(bond_density2(1:bondpoints,1:nbonds))
      allocate(bond_density_errorbar(1:bondpoints,1:nbonds))
      bond_density_elem = bondpoints * nbonds
      bond_density_red      = 0.d0
      bond_density2_red     = 0.d0
      bond_density          = 0.d0
      bond_density2         = 0.d0
      bond_density_errorbar = 0.d0

    END SUBROUTINE


    SUBROUTINE update_bonds(xx,bar_wfn_sq)

      real*8, intent(in) :: xx(1:ncart)
      real*8, intent(in) :: bar_wfn_sq
      integer            :: b,ix

      do b=1,nbonds
        call find_bond_index(xx,b,ix)
        bond_density(ix,b)  = bond_density(ix,b)  + bar_wfn_sq
        bond_density2(ix,b) = bond_density2(ix,b) + bar_wfn_sq**2
      enddo

    END SUBROUTINE


    SUBROUTINE find_bond_index(xx,b, ir)
    implicit none

      INTEGER, INTENT(OUT) :: ir
      REAL*8, INTENT(IN)   :: xx(ncart)
      INTEGER, INTENT(IN)  :: b
      INTEGER              :: i,j
      REAL*8               :: dist

      i = bond_pair(1,b)
      j = bond_pair(2,b)
      dist = SQRT( (xx(3*i - 2)-xx(3*j - 2))**2 +&
                   (xx(3*i - 1)-xx(3*j - 1))**2 +&
                   (xx(3*i    )-xx(3*j    ))**2  )

      ir = 1+INT(dist / bonddr(b))
      IF (ir > bondpoints  ) ir = bondpoints

    END SUBROUTINE


    SUBROUTINE MPI_REDUCE_BONDS()
      use mpimod
      implicit none
      integer :: err_mpi
      CALL MPI_REDUCE(bond_density,bond_density_red,bond_density_elem,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,err_mpi)
      CALL MPI_REDUCE(bond_density2, bond_density2_red,bond_density_elem,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,err_mpi)
    END SUBROUTINE


    SUBROUTINE print_bonds(density,density_err)
      use io_units_def
      use constants

      implicit none

      REAL*8, intent(in) :: density(bondpoints  , nbonds)
      REAL*8, intent(in) :: density_err(bondpoints  , nbonds)
      character(25) :: filename
      integer :: i, b

      DO b=1, nbonds

        !Preparing filename string
        WRITE(filename, "('bond_',1A,'.out')") adjustl(trim(bond_name(b)))

        !associate a unit to the cube file
        OPEN (UNIT=unit_bonds_out+b, FILE=trim(filename), STATUS='replace')

        ! write out cube file introduction
        WRITE(unit_bonds_out+b,*) "#Histogram density of bond. Values at the center of intervals (Ang unit)"

        WRITE(unit_bonds_out+b,"(3(1ES22.15,1x))") (bonddr(b)*(i-0.5d0)*FROMauTOang,density(i,b),density_err(i,b), i=1,bondpoints  )

        CLOSE(unit_bonds_out+b)

      ENDDO

    END SUBROUTINE print_bonds

    SUBROUTINE print_normalized_bonds(Nsteps_MC_tot,tot_int_red)
      integer, intent(in) :: Nsteps_MC_tot
      real*8,  intent(in) :: tot_int_red

      integer :: b
      !
      bond_density_red  = bond_density_red  / tot_int_red
      bond_density2_red = bond_density2_red / tot_int_red
      ! Evaluation of Monte Carlo uncertainty of the mean
      bond_density_errorbar = SQRT( abs(bond_density2_red - bond_density2_red**2) / Nsteps_MC_tot )
      !
      do b=1,nbonds
        bond_density_red(:,b)      = bond_density_red(:,b) / bondvol(b)
        bond_density_errorbar(:,b) = bond_density_errorbar(:,b) / bondvol(b)
      enddo

      !Check the normalization of each nucleus density
      DO b=1, nbonds
        print*, 'Check norm of bond distribution ',b,': ', sum(bond_density_red(:,b))*bondvol(b)
      END DO
      !
      ! Print out bond distributions on output files
      print*, 'Printing bond files'
      call print_bonds(bond_density_red,bond_density_errorbar)
      !
    END SUBROUTINE

END MODULE manage_bonds
