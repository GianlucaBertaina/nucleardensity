! COPYRIGHT (C) 2020 Chiara Donatella Aieta, Marco Micciarelli, Gianluca Bertaina, Michele Ceotto
! See LICENSE for details

!     This module reads input parameter from file and sets appropriate variables.
      MODULE input_def

      use constants
      use io_units_def 

      implicit none

      save
      !
      !
      !
!!!!!! parameters from input file !!!!!
      !
      integer :: nat, nrotrasl
      integer :: ncart, nvib
      integer :: dim_K
      integer :: Nsteps_MC
      integer :: nxpoints, nypoints, nzpoints
      integer :: print_mctraj
      character (len=80) :: file_geo, file_cnorm
      character (len=80) :: file_wfn, file_omega,file_bonds_input
      character (len=80) :: file_angles_input,file_dihedrals_input
      !
!!!!!! from file with eq geometry frame !!!!      
      !
      real*8, allocatable, dimension(:) :: xm_sqrt, x_eq_cart
      character (len=2), allocatable, dimension (:) :: symb_list
      integer, allocatable, dimension(:) :: Z_nuclei
      !
!!!!!! from file with harmonic frequencies !!!!      
      !
      real*8, allocatable, dimension(:) :: omega_vh
      !
!!!!!! from file with cnorm matrix !!!!      
      !
      real*8, allocatable, dimension(:,:) :: cnorm
      !
!!!!!! from file with wavefunction !!!!      
      !
      integer, allocatable, dimension(:,:) :: h_vec
      real*8, allocatable, dimension(:)   :: coef 
      !
!!!!!! Settings for bonds
      integer                       :: nbonds,bondpoints
      integer, allocatable          :: bond_pair(:,:)
      character(len=30),allocatable :: bond_name(:)
!!!!!! Settings for dihedrals
      integer                       :: ndihedrals,dihedralpoints
      integer, allocatable          :: dihedral_atoms(:,:)
      character(len=30),allocatable :: dihedral_name(:)
!!!!!! Settings for angles
      integer                       :: nangles,anglepoints
      integer, allocatable          :: angle_atoms(:,:)
      character(len=30),allocatable :: angle_name(:)
!!!!!! Options
      logical                       :: do_bonds=.false.
      logical                       :: do_densities=.false.
      logical                       :: do_dihedrals=.false.
      logical                       :: do_angles=.false.
      !
      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      subroutine get_input_params
      ! main subroutine calling the others in this module in the right order 
      ! Print status sentences commented out for parrallel execution
      !
      !
      implicit none
      !
      !print*, '    reading parameters from file input_nucleardensity.dat'
      call read_input
      !
      !print*, '    getting strcture masses from eq geometry frame:' 
      allocate(symb_list(nat))
      allocate(xm_sqrt(ncart))
      allocate(x_eq_cart(ncart))
      call get_geo_symb_mas(x_eq_cart,symb_list,xm_sqrt)
      !
      !print*, '    getting omegas:' 
      call get_omega
      !
      !print*, '    getting cnorm matrix:' 
      call get_cnorm
      !
      !print*, '    getting vibrational wavefunction:' 
      call get_wfn
      !
      if (do_bonds) call get_bonds()
      !
      if (do_angles) call get_angles()
      !
      if (do_dihedrals) call get_dihedrals()
      !
      !
      end subroutine get_input_params


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      subroutine read_input 

      implicit none

      character(30),parameter:: inputfile='input_nucleardensity.dat'

      logical :: file_exists=.true.

      !! check input file existence
      inquire(file=trim(inputfile), exist=file_exists)
      if (.not.file_exists) then
        print*, "Missing input file. Stopping program."
        STOP
      endif


      ! Read options from input file
      open(unit_input,file=trim(inputfile),status='old', action="read")
      !
      read(unit_input,*) ! Here a comment line
      read(unit_input,*) ! Here a comment line
      read(unit_input,*) nat,nrotrasl
      ncart = nat * 3
      nvib = ncart - nrotrasl
      read(unit_input,*) ! Here a comment line
      read(unit_input,*) dim_K ! eff harm basis size
      read(unit_input,*) ! Here a comment line
      read(unit_input,*) nxpoints, nypoints, nzpoints 
      read(unit_input,*) ! Here a comment line
      read(unit_input,*) Nsteps_MC
!!!!!!!!!!!!!!!  list of different filenames:  !!!!!!!!!!!!!!!!
      read(unit_input,*) ! Here a comment line
      read(unit_input,*) file_geo ! file with eq geometry 
      read(unit_input,*) ! Here a comment line
      read(unit_input,*) file_cnorm ! file with eq geometry 
      read(unit_input,*) ! Here a comment line
      read(unit_input,*) file_omega ! file with eq geometry 
      read(unit_input,*) ! Here a comment line
      read(unit_input,*) file_wfn ! file with SC wavefunction 
      read(unit_input,*) ! Here a comment line
      read(unit_input,*) print_mctraj ! Number of printed MC samples
      read(unit_input,*) ! Here a comment line
      read(unit_input,*) do_densities ! option
      read(unit_input,*) ! Here a comment line
      read(unit_input,*) file_bonds_input ! file with bonds settings
      read(unit_input,*) ! Here a comment line
      read(unit_input,*) file_angles_input ! file with angles settings
      read(unit_input,*) ! Here a comment line
      read(unit_input,*) file_dihedrals_input ! file with dihedrals settings
      !
      close(unit_input)

      !
      !
      if ((mod(nxpoints,2) == 0) .OR. &
          (mod(nypoints,2) == 0) .OR. &
          (mod(nzpoints,2) == 0)      ) then
        print*, "Number of grid points must be odd!"
        STOP
      end if
      !
      !

!!    Check all other input files existence
      !
      inquire(file=trim(file_geo), exist=file_exists)
      if (.not.file_exists) then
        print*, "Missing file with molecular geometry at equilibrium."
        print*, " Stopping program"
        STOP
      endif
      !
      inquire(file=trim(file_cnorm), exist=file_exists)
      if (.not.file_exists) then
        print*, "Missing file with cnorm matrix. Stopping program"
        STOP
      endif
      !
      inquire(file=trim(file_omega), exist=file_exists)
      if (.not.file_exists) then
        print*, "Missing file with harmonic frequencies." 
        print*, "Stopping program"
        STOP
      endif
      !
      !
      inquire(file=trim(file_wfn), exist=file_exists)
      if (.not.file_exists) then
        print*, "Missing file with wavefunction. Stopping program"
        STOP
      endif
      !
      if (trim(file_bonds_input)=="F") then
        do_bonds=.false.
      else
        do_bonds=.true.
        inquire(file=trim(file_bonds_input), exist=file_exists)
        if (.not.file_exists) then
          print*, "Missing file with bond definitions. Stopping program"
          STOP
        endif
      endif
      !
      if (trim(file_angles_input)=="F") then
        do_angles=.false.
      else
        do_angles=.true.
        inquire(file=trim(file_angles_input), exist=file_exists)
        if (.not.file_exists) then
          print*, "Missing file with angle definitions. Stopping program"
          STOP
        endif
      endif
      !
      if (trim(file_dihedrals_input)=="F") then
        do_dihedrals=.false.
      else
        do_dihedrals=.true.
        inquire(file=trim(file_dihedrals_input), exist=file_exists)
        if (.not.file_exists) then
          print*, "Missing file with dihedral definitions. Stopping program"
          STOP
        endif
      endif
      !
      !
      end subroutine read_input


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      subroutine get_geo_symb_mas(x_eq_cart,symb_list,xm_sqrt)

      use io_units_def,    only : unit_eq_geo
      use constants

      implicit none

      real*8, dimension(nat)   :: amass
      real*8, dimension(ncart) :: xm

      integer :: i, j, k

      character (len=2), dimension (nat), intent(out) :: symb_list
      real*8, dimension(ncart), intent(out) :: x_eq_cart, xm_sqrt

      ! file with equilibrium geometry in xyz format:
      open(unit_eq_geo, file= file_geo, status='old', action="read") 

      ! allocate nuclear mass vector for output cube files
      allocate(Z_nuclei(nat))

      read (unit_eq_geo,*)
      read (unit_eq_geo,*)

      do i = 1, nat

          read (unit_eq_geo,*) symb_list(i),  x_eq_cart(3*i-2:3*i)

          select case (symb_list(i))

          case ('H')
              amass(i) = 1837.15d0
              Z_nuclei(i) = 1
          case ('D')  !Deuterium
              amass(i) = 3671.48d0
              Z_nuclei(i) = 1
          case ('O')
              amass(i) = 29156.96d0
              Z_nuclei(i) = 8
          case ('C')
              amass(i) = 21874.66d0
              Z_nuclei(i) = 6
          case ('N')
              amass(i) = 25526.06d0
              Z_nuclei(i) = 7
          case ('Ti')
              amass(i) = 87256.20d0
              Z_nuclei(i) = 22
          case default
               write(*,*) 'One or more atoms in your molecule is not', &
                          ' present in our database'
               stop

          end select

      enddo

      close(unit_eq_geo)

      ! convert Angstrom (as in xyz format) to au
      x_eq_cart(:) = x_eq_cart(:) / FROMauTOang 


!     cartesian masses
      k=0
      do i = 1, nat
         do  j = 1, 3  ! over three dimensional x,y,z

            k = k + 1
            xm( k ) = amass( i )

         enddo
      enddo
      xm_sqrt=dsqrt(xm)


      end subroutine get_geo_symb_mas


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      subroutine get_omega

      use io_units_def,    only : unit_omega

      implicit none

      integer :: i
      !
      !
      allocate(omega_vh(nvib))
      !
      !
      open(unit_omega,file=trim(file_omega),status='old', action="read")
      !
      !
      read(unit_omega,*) ! A comment line
      do i = 1, nvib
        read(unit_omega,*) omega_vh(i)
      enddo
      !
      omega_vh = omega_vh * FROMcmTOau
      !
      close(unit_omega)

      end subroutine get_omega



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


      subroutine get_cnorm

      use io_units_def,    only : unit_cnorm

      implicit none

      integer :: i
      !
      !
      allocate(cnorm(ncart,ncart))
      !
      !
      open(unit_cnorm,file=trim(file_cnorm),status='old', action="read")
      !
      !
      do i= 1, ncart
          !write(unit,"(9g15.6)") cnorm(i,:)
          read(unit_cnorm,*) cnorm(:,i)
      enddo
      !
      !
      close(unit_cnorm)

      end subroutine get_cnorm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


      subroutine get_wfn

      use io_units_def,    only : unit_wfn  

      implicit none

      integer :: K
      real*8  :: indx
      !
      !
      allocate(h_vec(nvib,dim_K))
      allocate(coef(dim_K))
      !
      !
      open(unit_wfn,file=trim(file_wfn),status='old', action="read")
      !
      !
      read(unit_wfn,*) ! seven lines for comments
      read(unit_wfn,*)
      read(unit_wfn,*)
      read(unit_wfn,*)
      read(unit_wfn,*)
      read(unit_wfn,*)
      read(unit_wfn,*)
      do K = 1, dim_K
          read(unit_wfn,*) indx, h_vec(:,K), coef(K)
      enddo
      !
      !
      close(unit_wfn)
      !
      !
      end subroutine get_wfn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      subroutine get_bonds()
        use io_units_def,    only : unit_bonds_in
        implicit none
        integer :: b
        !
        open(unit_bonds_in,file=trim(file_bonds_input),status='old', action="read")
        read(unit_bonds_in,*) ! Info
        read(unit_bonds_in,*) bondpoints
        read(unit_bonds_in,*) ! Info
        read(unit_bonds_in,*) nbonds
        allocate(bond_pair(1:2,1:nbonds))
        allocate(bond_name(1:nbonds))

        do b=1,nbonds
          read(unit_bonds_in,*) ! Info
          read(unit_bonds_in,*) bond_name(b)
          read(unit_bonds_in,*) ! Info
          read(unit_bonds_in,*) bond_pair(1:2,b)
        enddo

        close(unit_bonds_in)

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_angles()
        use io_units_def,    only : unit_angles_in
        implicit none
        integer :: b
        !
        open(unit_angles_in,file=trim(file_angles_input),status='old', action="read")
        read(unit_angles_in,*) ! Info
        read(unit_angles_in,*) anglepoints
        read(unit_angles_in,*) ! Info
        read(unit_angles_in,*) nangles
        allocate(angle_atoms(1:3,1:nangles))
        allocate(angle_name(1:nangles))

        do b=1,nangles
          read(unit_angles_in,*) ! Info
          read(unit_angles_in,*) angle_name(b)
          read(unit_angles_in,*) ! Info
          read(unit_angles_in,*) angle_atoms(1:3,b)
        enddo

        close(unit_angles_in)

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_dihedrals()
        use io_units_def,    only : unit_dihedrals_in
        implicit none
        integer :: b
        !
        open(unit_dihedrals_in,file=trim(file_dihedrals_input),status='old', action="read")
        read(unit_dihedrals_in,*) ! Info
        read(unit_dihedrals_in,*) dihedralpoints
        read(unit_dihedrals_in,*) ! Info
        read(unit_dihedrals_in,*) ndihedrals
        allocate(dihedral_atoms(1:4,1:ndihedrals))
        allocate(dihedral_name(1:ndihedrals))

        do b=1,ndihedrals
          read(unit_dihedrals_in,*) ! Info
          read(unit_dihedrals_in,*) dihedral_name(b)
          read(unit_dihedrals_in,*) ! Info
          read(unit_dihedrals_in,*) dihedral_atoms(1:4,b)
        enddo

        close(unit_dihedrals_in)

      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  END MODULE input_def






