   PROGRAM NUDE

      use constants
      use io_units_def
      use input_def
      use manage_cubes
      use manage_bonds

      implicit none

      include "mpif.h"

      REAL*8, ALLOCATABLE, DIMENSION(:)       :: q_eq, sigma_gaus 
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: qq, xx 
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: q_expect_value, q_expect_value_h
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: x_expect_value, x_expect_value_h
      !
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: density !, density_sq , density_errorbar
      REAL*8, ALLOCATABLE, DIMENSION(:,:)     :: bond_density, bond_density2,bond_density_errorbar

      integer :: i, ix,iy,iz, istep, ii,b
      real*8  :: Ri(3)
      real*8  :: bar_wfn_sq, tot_int, tot_int_sq
      !
      real*8  :: det_omega, norm_gaussian_normal
      !
      integer :: idum
      real*8  :: ran2
      real*8  :: xrand1,xrand2,y1,y2

!!!!!!!MPI variables
      INTEGER :: my_rank, num_procs, err_mpi
      REAL*8  :: t_start, t_end

      CHARACTER(len=4)  :: str1
      CHARACTER(len=10) :: str
      CHARACTER(len=20) :: filename

      INTEGER :: density_elem, reminder, Nsteps_MC_tot,bond_density_elem
      REAL*8  :: tot_int_red, tot_int_sq_red
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: q_expect_value_red, q_expect_value_h_red
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: density_red  !, density_sq_red
      REAL*8, ALLOCATABLE, DIMENSION(:,:)     :: bond_density_red,bond_density2_red
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
      !
      IF (my_rank == 0) print*, 'setting parameters from input files'
      call get_input_params
      IF (my_rank == 0) THEN
         print*,'done'
         print*, '  '
         print*, '  '
      END IF
      !
      !
      IF (my_rank == 0) print*, 'setting grid and parameters for cube files'
      call set_cube
      !
      IF (my_rank == 0) print*, "setting bonds"
      if (do_bonds) call set_bonds
      !
      IF (my_rank == 0) THEN
         print*,'done'
         print*, '  '
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
      if (do_densities) then
        ! Initialize density (a 3D array for each atom)
        allocate(density(nxpoints, nypoints, nzpoints, nat))
        !allocate(density_sq(nxpoints, nypoints, nzpoints, nat))
        allocate(density_red(nxpoints, nypoints, nzpoints, nat))
        !allocate(density_sq_red(nxpoints, nypoints, nzpoints, nat))
        !allocate(density_errorbar(nxpoints, nypoints, nzpoints, nat))
        density_elem = nxpoints*nypoints*nzpoints*nat
      else
        allocate(density(1,1,1,1))
        allocate(density_red(1,1,1,1))
        density_elem = 1
      endif
      density           = 0.d0
      density_red       = 0.d0
      !density_sq       = 0.d0
      !density_sq_red   = 0.d0
      !density_errorbar = 0.d0
      !
      ! Initialize bond density (1D array for each bond specified in input)
      if (do_bonds) then
        allocate(bond_density_red(1:bondpoints,1:nbonds))
        allocate(bond_density2_red(1:bondpoints,1:nbonds))
        allocate(bond_density(1:bondpoints,1:nbonds))
        allocate(bond_density2(1:bondpoints,1:nbonds))
        allocate(bond_density_errorbar(1:bondpoints,1:nbonds))
        bond_density_elem = bondpoints * nbonds
      else
        allocate(bond_density_red(1,1))
        allocate(bond_density2_red(1,1))
        allocate(bond_density(1,1))
        allocate(bond_density2(1,1))
        allocate(bond_density_errorbar(1,1))
        bond_density_elem = 1
      endif
      bond_density_red      = 0.d0
      bond_density2_red     = 0.d0
      bond_density          = 0.d0
      bond_density2         = 0.d0
      bond_density_errorbar = 0.d0
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
        WRITE(*,*) 'All processes but rank 0 do ', Nsteps_MC, ' steps.'
        Nsteps_MC = Nsteps_MC + reminder
        WRITE(*,*) 'Rank 0 does ', Nsteps_MC, ' steps.'
      ENDIF

      ! MC SAMPLING CYCLE
      tot_int=0.d0
      tot_int_sq=0.d0
      DO istep=1,Nsteps_MC

        !TODO: valutare se far stamapare questo a ogni thread
        !IF ( MOD(istep,1000)==0 ) THEN
        !   WRITE(*,*) istep, ' steps of ', Nsteps_MC, ' done!'
        !END IF

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
        IF (switch_harm /= 1) THEN
          call bar_wfn_squared(qq-q_eq,bar_wfn_sq)
        ELSE
          bar_wfn_sq = 1.d0
        END IF
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
        if (do_densities) then
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
        endif
        !
        ! Update bond densities
        if (do_bonds) then
          do b=1,nbonds
            call find_bond_index(xx,b,ix)
            bond_density(ix,b)  = bond_density(ix,b)  + bar_wfn_sq
            bond_density2(ix,b) = bond_density2(ix,b) + bar_wfn_sq**2
          enddo
        endif
        !
      ENDDO


!!!!!!!MPI reduce
      CALL MPI_REDUCE(tot_int, tot_int_red, 1, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD, err_mpi)
      CALL MPI_REDUCE(tot_int_sq, tot_int_sq_red, 1, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD, err_mpi)
      CALL MPI_REDUCE(q_expect_value, q_expect_value_red, ncart, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD, err_mpi)
      CALL MPI_REDUCE(q_expect_value_h, q_expect_value_h_red, ncart, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD, err_mpi)
      !
      if (do_densities) then
        CALL MPI_REDUCE(density, density_red, density_elem, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD, err_mpi)
        !CALL MPI_REDUCE(density_sq, density_sq_red, density_elem, MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD, err_mpi)
      endif
      !
      if (do_bonds) then
        CALL MPI_REDUCE(bond_density,bond_density_red,bond_density_elem,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,err_mpi)
        CALL MPI_REDUCE(bond_density2, bond_density2_red,bond_density_elem,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,err_mpi)
      endif
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
        q_expect_value_red(4:9)   = q_eq(4:9)
        q_expect_value_h_red(4:9) = q_eq(4:9)
        !print*, 'q_eq    = ', q_eq
        !print*, 'q_exp_h = ', q_expect_value_h_red
        !print*, 'q_exp   = ', q_expect_value_red

        ! Bring back structures in cartesian coordinates
        call to_cartesian(x_expect_value,   q_expect_value_red)
        call to_cartesian(x_expect_value_h ,q_expect_value_h_red)
        !
        ! Print out structures in xyz format on file
        open(unit_geo_exp,   file='expectation.xyz')
        open(unit_geo_exp_h, file='expectation_H.xyz')
        !
        write(unit_geo_exp,* ) nat
        write(unit_geo_exp,* )
        write(unit_geo_exp_h,* ) nat
        write(unit_geo_exp_h,* )
        do i = 1, nat
          write(unit_geo_exp,* )   symb_list(i), x_expect_value(3*i-2:3*i)   * FROMauTOang
          write(unit_geo_exp_h,* ) symb_list(i), x_expect_value_h(3*i-2:3*i) * FROMauTOang
        enddo
        !
        close(unit_geo_exp)
        close(unit_geo_exp_h)
        !
        det_omega = product(omega_vh(1:nvib))
        !
        norm_gaussian_normal = DSQRT( (2*pi)**nvib * 0.5d0 / det_omega )
        !
        ! Normalize density:
        print*, 'Total integral of wf^2, Nsteps, Gauss norm: ', tot_int_red , Nsteps_MC_tot ,norm_gaussian_normal
        !
        if (do_densities) then
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
          !Check the normalization of each nucleus density
          DO ii=1, nat
            print*, ii, ' unnormalized nucleus norm: ', sum(density_red(:,:,:,ii))
          END DO
          !
          !Normalize by the voxel dimension
          density_red = density_red / (dx*dy*dz)
          !density_errorbar = density_errorbar / (dx*dy*dz)
          !
          ! Print out nuclear densities on output files
          print*, 'Printing cube files'
          call print_cube(density_red)
          !
          ! print*, 'printing cube files with density error bar'
          ! call print_cube_errorbar(density_errorbar)
        endif
        !
        if (do_bonds) then
          !
          bond_density_red  = bond_density_red  / tot_int_red
          bond_density2_red = bond_density2_red / tot_int_red
          bond_density_errorbar = SQRT( abs(bond_density2_red - bond_density2_red**2) / Nsteps_MC_tot )
          !
          do b=1,nbonds
            bond_density_red(:,b)      = bond_density_red(:,b) / bondvol(b)
            bond_density_errorbar(:,b) = bond_density_errorbar(:,b) / bondvol(b)
          enddo
          !
          ! Print out bond distributions on output files
          print*, 'Printing bond files'
          call print_bonds(bond_density_red,bond_density_errorbar)
          !
        endif
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

   END PROGRAM NUDE


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
           implicit none
      
           real(8) :: x1,x2,y1,y2,u,sigma
           REAL(8),PARAMETER :: pi = 3.141592653589793238462643383279503d0
      
           y1 = sqrt(-2.0d0*log(x1)) * cos( 2.0d0*pi*x2 )
           y2 = sqrt(-2.0d0*log(x1)) * sin( 2.0d0*pi*x2 )
      
           y1 = y1*sigma + u
           y2 = y2*sigma + u
      
      end subroutine Husimi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   
      FUNCTION ran2(idum)
        INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
        REAL*8 ran2,AM,EPS,RNMX
        PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
              IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211, &
              IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)

        INTEGER idum2,j,k,iv(NTAB),iy
        SAVE iv,iy,idum2
        DATA idum2/123456789/, iv/NTAB*0/, iy/0/
        if (idum.le.0) then
           idum=max(-idum,1)
           idum2=idum
           do j=NTAB+8 , 1 ,-1
              k = idum / IQ1
              idum=IA1*(idum-k*IQ1)-k*IR1
              if (idum.lt.0) idum=idum+IM1
              if (j.le.NTAB) iv(j)=idum
           enddo
           iy=iv(1)
         endif
         k=idum/IQ1
         idum=IA1*(idum-k*IQ1)-k*IR1
         if (idum.lt.0) idum=idum+IM1
         k=idum2/IQ2
         idum2=IA2*(idum2-k*IQ2)-k*IR2
         if (idum2.lt.0) idum2=idum2+IM2
         j=1+iy/NDIV
         iy=iv(j)-idum2
         iv(j)=idum
         if(iy.lt.1)iy=iy+IMM1
         ran2=min(AM*iy,RNMX)
         return
      END function ran2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


















