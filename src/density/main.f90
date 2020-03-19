! COPYRIGHT (C) 2020 Chiara Donatella Aieta, Marco Micciarelli, Gianluca Bertaina, Michele Ceotto
! See LICENSE for details

   PROGRAM NUCLEARDENSITY

      use constants
      use io_units_def
      use input_def
      use manage_cubes
      use manage_bonds
      use manage_angles
      use manage_dihedrals
      use mpimod

      implicit none

      REAL*8, ALLOCATABLE, DIMENSION(:)       :: q_eq, sigma_gaus 
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: qq, xx 
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: q_expect_value, q_expect_value_h
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: x_expect_value, x_expect_value_h
      !
      integer :: i, istep
      real*8  :: bar_wfn_sq, tot_int, tot_int_sq
      !
      real*8  :: det_omega, norm_gaussian_normal
      !
      integer :: idum
      real*8  :: ran2
      real*8  :: xrand1,xrand2,y1,y2

!!!!!!!MPI variables
      INTEGER :: err_mpi
      REAL*8  :: t_start, t_end

      CHARACTER(len=4)  :: str1
      CHARACTER(len=10) :: str
      CHARACTER(len=20) :: filename

      INTEGER :: reminder, Nsteps_MC_tot
      REAL*8  :: tot_int_red, tot_int_sq_red
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: q_expect_value_red, q_expect_value_h_red
!!!!!!!
      !
      CALL MPI_INIT(err_mpi)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, err_mpi)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, err_mpi)
      t_start = MPI_WTIME()
      !
      idum = -200 + my_rank
      !
      !
      IF (my_rank == 0) THEN
         print*, 'Entering program:  NUclear DEnsity:'
         print*, '  '
         print*, '  '
      END IF
      !
      IF (my_rank == 0) print*, 'Setting parameters from input files'
      call get_input_params
      IF (my_rank == 0) THEN
         print*,'Done'
         print*, '  '
      END IF
      !
      if (do_densities) then
        call set_cube
        IF (my_rank == 0) print*, 'Setting grid and parameters for cube files'
      endif
      !
      if (do_bonds) then
        call set_bonds
        IF (my_rank == 0) print*, "Setting bonds"
      endif
      !
      if (do_angles) then
        call set_angles
        IF (my_rank == 0) print*, "Setting angles"
      endif
      !
      if (do_dihedrals) then
        call set_dihedrals
        IF (my_rank == 0) print*, "Setting dihedrals"
      endif
      !
      IF (my_rank == 0) THEN
         print*,'Done'
         print*, '  '
      END IF
      !
      !  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! CALCULATION OF NUCLEAR DENSITY    !
      ! USING MONTE CARLO                 !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! Equlibrium coordinates in normal modes:
      allocate(q_eq(ncart))
      call to_normal(x_eq_cart,q_eq)
      !
      if (do_densities) call allocate_densities()
      !
      ! Initialize bond density (1D array for each bond specified in input)
      if (do_bonds) call allocate_bonds()
      !
      ! Initialize angle density (1D array for each angle specified in input)
      if (do_angles) call allocate_angles()
      !
      ! Initialize dihedral density (1D array for each dihedral specified in input)
      if (do_dihedrals) call allocate_dihedrals()
      !
      ! Set multivariate gaussian width vector 
      ! for husimi distribution
      allocate(sigma_gaus(nvib))
      do i=1,nvib
        sigma_gaus(i) = sqrt(1.d0/(2.d0 * omega_vh(i)))
      enddo
      !
      ! Initialize structure to generate to the equilibrium
      ! note: this is to leave unchanged rototraslation coordinates !
      allocate(qq(ncart))
      qq = q_eq  
      allocate(q_expect_value(ncart)) 
      q_expect_value = 0.d0  
      allocate(q_expect_value_h(ncart))
      q_expect_value_h = 0.d0  
      allocate(q_expect_value_red(ncart)) 
      q_expect_value_red = 0.d0  
      allocate(q_expect_value_h_red(ncart))
      q_expect_value_h_red = 0.d0  
      !
      ! allocate corresponding structure in cartesian coordinates
      allocate(xx(ncart))
      allocate(x_expect_value(ncart)) 
      allocate(x_expect_value_h(ncart))
      !
      ! open unit for fil with MC traj in xyz
      IF (switch_print_mctraj == 1) THEN
        str1='traj'
        WRITE (str, '(I10)') my_rank
        filename = str1 // TRIM(ADJUSTL(str)) //'.xyz'
        open(unit_trajMC+my_rank, file=filename)
      END IF
      !
      !WORK DIVISION AMONG THREADS:
      !Each rank nkows how many MC steps has to do
      Nsteps_MC_tot = Nsteps_MC
      reminder = MOD(Nsteps_MC, num_procs)
      Nsteps_MC = (Nsteps_MC - reminder) / num_procs
      IF (my_rank == 0) then
        WRITE(*,*) 'All ',num_procs,' processes but rank 0 do ', Nsteps_MC, ' steps.'
        Nsteps_MC = Nsteps_MC + reminder
        WRITE(*,*) 'Rank 0 does ', Nsteps_MC, ' steps.'
        flush(6)
      ENDIF

      ! MC SAMPLING CYCLE
      tot_int=0.d0
      tot_int_sq=0.d0
      DO istep=1,Nsteps_MC

        !
        ! Generate molecular structure in normal coordinates:
        do i = 1, nvib
          xrand1 = ran2(idum)
          xrand2 = ran2(idum)
          call Husimi(xrand1,xrand2,y1,y2,q_eq(i),sigma_gaus(i))
          qq(i) = y1
        enddo
        !
        ! Compute wfn^2/Gaussian at point qq
        call bar_wfn_squared(qq-q_eq,bar_wfn_sq)
        !
        tot_int = tot_int + bar_wfn_sq
        tot_int_sq = tot_int_sq + bar_wfn_sq**2
        !
        ! Bring back structure in cartesian coordinates
        call to_cartesian(xx,qq)
        !
        ! Print out structure in xyz format on file
        IF (Nsteps_MC <= 100000 .and. switch_print_mctraj==1) THEN
          write(unit_trajMC+my_rank,* ) nat
          write(unit_trajMC+my_rank,* )
          do i = 1, nat
            write(unit_trajMC+my_rank,* ) symb_list(i), xx(3*i-2:3*i) * FROMauTOang
          enddo
        END IF
        !
        !
        ! Update value for expectation postition
        q_expect_value   = q_expect_value   + qq * bar_wfn_sq
        q_expect_value_h = q_expect_value_h + qq
        !
        !
        ! Update corresponding densities
        if (do_densities) call update_densities(xx,bar_wfn_sq)
        !
        ! Update bond densities
        if (do_bonds) call update_bonds(xx,bar_wfn_sq)
        !
        ! Update angle densities
        if (do_angles) call update_angles(xx,bar_wfn_sq)
        !
        ! Update dihedral densities
        if (do_dihedrals) call update_dihedrals(xx,bar_wfn_sq)
        !
      ENDDO


!!!!!!!MPI reduce
      CALL MPI_REDUCE(tot_int, tot_int_red, 1, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD, err_mpi)
      CALL MPI_REDUCE(tot_int_sq, tot_int_sq_red, 1, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD, err_mpi)
      CALL MPI_REDUCE(q_expect_value, q_expect_value_red, ncart, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD, err_mpi)
      CALL MPI_REDUCE(q_expect_value_h, q_expect_value_h_red, ncart, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD, err_mpi)
      !
      if (do_densities) CALL MPI_REDUCE_DENSITIES()
      !
      if (do_bonds) CALL MPI_REDUCE_BONDS()
      !
      if (do_angles) CALL MPI_REDUCE_ANGLES()
      !
      if (do_dihedrals) CALL MPI_REDUCE_DIHEDRALS()
!!!!!!!
      ! close unit for fil with MC traj in xyz
      close(unit_trajMC+my_rank)
      !
      IF (my_rank == 0) THEN
        !
        ! Normalize integral for expectation structures:
        q_expect_value_red   = q_expect_value_red   / tot_int_red
        q_expect_value_h_red = q_expect_value_h_red / Nsteps_MC_tot
        !
        q_expect_value_red(nvib+1:ncart)   = q_eq(nvib+1:ncart)
        q_expect_value_h_red(nvib+1:ncart) = q_eq(nvib+1:ncart)

        ! Bring back structures in cartesian coordinates
        call to_cartesian(x_expect_value,   q_expect_value_red)
        call to_cartesian(x_expect_value_h ,q_expect_value_h_red)
        !
        ! Print out structures in xyz format on file
        open(unit_geo_exp,   file='expectation.xyz')
        open(unit_geo_exp_h, file='expectation_H.xyz')
        !
        write(unit_geo_exp,* ) nat
        write(unit_geo_exp,* ) "# Expectation value of geometry in the chosen wavefunction"
        write(unit_geo_exp_h,* ) nat
        write(unit_geo_exp_h,* ) "# Expectation value of geometry in harmonic ground state"
        do i = 1, nat
          write(unit_geo_exp,* )   symb_list(i), x_expect_value(3*i-2:3*i)   * FROMauTOang
          write(unit_geo_exp_h,* ) symb_list(i), x_expect_value_h(3*i-2:3*i) * FROMauTOang
        enddo
        !
        close(unit_geo_exp)
        close(unit_geo_exp_h)
        !
        ! Norm of Gaussian: not used
        det_omega = product(omega_vh(1:nvib))
        norm_gaussian_normal = DSQRT( (2*pi)**nvib * 0.5d0 / det_omega )
        !
        ! Normalize density:
        print*, 'Total integral of wf^2, Nsteps: ', tot_int_red , Nsteps_MC_tot
        !
        if (do_densities) call print_normalized_densities(Nsteps_MC_tot,tot_int_red)
        !
        if (do_bonds) call print_normalized_bonds(Nsteps_MC_tot,tot_int_red)
        !
        if (do_angles) call print_normalized_angles(Nsteps_MC_tot,tot_int_red)
        !
        if (do_dihedrals) call print_normalized_dihedrals(Nsteps_MC_tot,tot_int_red)
        !
        print*,'Done'
        print*, '  '
        !
        ! call cpu_time(finish)
        t_end = MPI_WTIME()
        print '("Time for execution = ",f10.3," seconds.")', t_end-t_start
        print*, '  '
        !
      END IF
       
      call MPI_FINALIZE(err_mpi)

   END PROGRAM NUCLEARDENSITY


!
! From here there are the subroutines used to manipulate coordinates 
! in order to go from normal modes to cartesian and vice versa.
! NOTE 1: the knoledge of the cnorm matrix
! throght the which the tranformations are performed
! is assumed 
! NOTE 2: all coordinates are assumed to be in au
!

      SUBROUTINE to_normal(RR, QQ)

      use input_def,    only: cnorm, xm_sqrt, ncart  

      implicit none

      real*8, dimension(ncart), intent(in)       :: RR
      real*8, dimension(ncart), intent(out)      :: QQ

      real*8, dimension(ncart)                   :: RR_mas
      integer :: i,j

      RR_mas = RR * xm_sqrt

      QQ = 0.d0
      do i = 1, ncart
          do j = 1, ncart
                 QQ(i) = QQ(i) + cnorm(j,i) * RR_mas(j)
          enddo
      enddo

      
      END SUBROUTINE to_normal



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      SUBROUTINE to_cartesian(RR, QQ)

      use input_def,    only: cnorm, xm_sqrt, ncart    

      implicit none

      real*8, dimension(ncart), intent(in)       :: QQ
      real*8, dimension(ncart), intent(out)      :: RR

      real*8, dimension(ncart)                   :: RR_mas
      integer :: i,j


      RR_mas = 0.d0
      do i = 1, ncart
          do j = 1, ncart
                 RR_mas(i) = RR_mas(i) + cnorm(i,j) * QQ(j)
          enddo
      enddo

      RR = RR_mas / xm_sqrt


      END SUBROUTINE to_cartesian


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine Husimi(x1,x2,y1,y2,u,sigma)
      use constants
           implicit none
      
           real(8) :: x1,x2,y1,y2,u,sigma
      
           y1 = sqrt(-2.0d0*log(x1)) * cos( 2.0d0*pi*x2 )
           y2 = sqrt(-2.0d0*log(x1)) * sin( 2.0d0*pi*x2 )
      
           y1 = y1*sigma + u
           y2 = y2*sigma + u
      
      end subroutine Husimi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


















