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
      integer :: switch_harm
      character (len=80) :: file_geo, file_cnorm
      character (len=80) :: file_wfn, file_omega
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
      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      subroutine get_input_params
      ! main subroutine calling the others in this module in the right order 
      ! Print status sentences commented out for parrallel execution
      !
      !
      implicit none
      !
      !print*, '    reading parameters from file input_nude.dat' 
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
      !
      end subroutine get_input_params


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      subroutine read_input 


      implicit none

      character(30),parameter:: inputfile='input_nude.dat'
      integer :: i,j,k

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
      read(unit_input,*) switch_harm! Harmonic 1 / Hanaramonic /=1 
      !
      close(unit_input)

      !
      !

!!    Check that the number of grid points in all directions is even
      !if (                            &
	!  (mod(nxpoints,2) /= 0) .OR. &
        !  (mod(nypoints,2) /= 0) .OR. &
        !  (mod(nzpoints,2) /= 0)      & 
	!  ) then
        !print*, "number of grid points must be even!"
	!STOP
      !endif

      if (                            &
	  (mod(nxpoints,2) == 0) .OR. &
          (mod(nypoints,2) == 0) .OR. &
          (mod(nzpoints,2) == 0)      & 
	  ) then
        print*, "number of grid points must be odd!"
	STOP
      end if
      !
      !

!!    Check all other input files existence
      !
      inquire(file=trim(file_geo), exist=file_exists)
      if (.not.file_exists) then
        print*, "Missing file with molecular geometry at equilibrium."
        print*,  "Stopping program"
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
      !
      !
      allocate(h_vec(nvib,dim_K))
      allocate(coef(dim_K))
      !
      !
      open(unit_wfn,file=trim(file_wfn),status='old', action="read")
      !
      !
      read(unit_wfn,*) ! two lines for comments
      read(unit_wfn,*)
      do K = 1, dim_K
          read(unit_wfn,*) h_vec(:,K), coef(K)
      enddo
      !
      !
      close(unit_cnorm)
      !
      !
      end subroutine get_wfn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      END MODULE input_def






