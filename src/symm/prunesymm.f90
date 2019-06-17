PROGRAM prunesymm
  implicit none

  integer, parameter   :: unit_input=10,unit_wfn=20
  character (len=80)   :: file_wfn, file_char, file_wfn_pruned

  integer              :: nvib,nat,ncart,nrotrasl

  integer              :: dim_K_in,dim_K_out
  integer, allocatable :: h_vec_in(:,:),h_vec_out(:,:)
  real*8,  allocatable :: coef_in(:),coef_out(:)
  real*8,  allocatable :: indx_in(:),indx_out(:)

  integer, allocatable :: characters(:),reference(:)
  integer              :: refcharacter
  real*8               :: norm

  call read_input
  !
  call get_wfn
  !
  call get_characters
  !
  refcharacter = calc_character(reference)
  !
  call prune_wavefunction
  !
  ! Normalize coefficients
  norm = sqrt(sum(coef_out(1:dim_K_out)**2))
  coef_out(1:dim_K_out) = coef_out(1:dim_K_out) / norm
  !
  call print_pruned_wavefunction()
  !
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_input
    implicit none

    character(30),parameter:: inputfile='input_prunesymm.dat'
    logical :: file_exists=.true.

    !! check input file existence
    inquire(file=trim(inputfile), exist=file_exists)
    if (.not.file_exists) then
      print*, "Missing prunesymm input file. Stopping program."
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
    allocate(reference(1:nvib))
    read(unit_input,*) ! Here a comment line
    read(unit_input,*) dim_K_in ! eff initial harm basis size
    read(unit_input,*) ! Here a comment line
    read(unit_input,*) file_wfn ! file with SC wavefunction
    read(unit_input,*) ! Here a comment line
    read(unit_input,*) file_wfn_pruned ! file with symmetry pruned wavefunction
    read(unit_input,*) ! Here a comment line
    read(unit_input,*) reference(1:nvib) ! quantum harmonic numbers of the state of interest
    read(unit_input,*) ! Here a comment line
    read(unit_input,*) file_char ! file with the characters of each normal mode
    !
    close(unit_input)
    !
    ! Check all other input files existence
    inquire(file=trim(file_wfn), exist=file_exists)
    if (.not.file_exists) then
      print*, "Missing file with initial wavefunction. Stopping program"
      STOP
    endif
    !
    inquire(file=trim(file_char), exist=file_exists)
    if (.not.file_exists) then
      print*, "Missing file with characters of normal modes. Stopping program"
      STOP
    endif
    !
  end subroutine read_input
  !
  !
  subroutine get_wfn
    implicit none
    !
    integer :: K
    !
    allocate(h_vec_in(nvib,dim_K_in),h_vec_out(nvib,dim_K_in))
    allocate(coef_in(dim_K_in),coef_out(dim_K_in))
    allocate(indx_in(dim_K_in),indx_out(dim_K_in))
    !
    open(unit_wfn,file=trim(file_wfn),status='old', action="read")
    !
    read(unit_wfn,*) ! seven lines for comments
    read(unit_wfn,*)
    read(unit_wfn,*)
    read(unit_wfn,*)
    read(unit_wfn,*)
    read(unit_wfn,*)
    read(unit_wfn,*)
    do K = 1, dim_K_in
      read(unit_wfn,*) indx_in(K), h_vec_in(1:nvib,K), coef_in(K)
    enddo
    !
    close(unit_wfn)
  end subroutine
  !
  !
  subroutine get_characters
    implicit none
    !
    integer :: i
    !
    allocate(characters(nvib))
    !
    open(unit_wfn,file=trim(file_char),status='old', action="read")
    !
    read(unit_wfn,*) ! one line for comments
    do i = 1, nvib
      read(unit_wfn,*) characters(i)
    enddo
    !
    close(unit_wfn)
  end subroutine
  !
  !
  function calc_character(state)
    integer             :: calc_character
    integer, intent(in) :: state(1:nvib)
    integer             :: i
    !
    calc_character = 1
    do i=1,nvib
      calc_character = calc_character * characters(i)**state(i)
    enddo
    !
  end function
  !
  !
  subroutine prune_wavefunction
    implicit none
    !
    integer :: cha,K
    !
    dim_K_out = 0
    do K = 1, dim_K_in
      cha = calc_character(h_vec_in(1:nvib,K))
      !
      if (cha == refcharacter) then
        !
        dim_K_out = dim_K_out + 1
        h_vec_out(1:nvib,dim_K_out) = h_vec_in(1:nvib,K)
        indx_out(dim_K_out)         = indx_in(K)
        coef_out(dim_K_out)         = coef_in(K)
        !
      endif
    enddo
    !
  end subroutine
  !
  !
  subroutine print_pruned_wavefunction
    implicit none
    !
    integer :: K
    !
    open(unit_wfn,file=trim(file_wfn_pruned),status='unknown', action="write")
    !
    Write(unit_wfn,'("# Pruned wavefunction, with character: ",1I3)') refcharacter
    Write(unit_wfn,'("# Number of resulting basis elements: ",1I8)') dim_K_out
    Write(unit_wfn,'("# ")')
    Write(unit_wfn,'("# ")')
    Write(unit_wfn,'("# ")')
    Write(unit_wfn,'("# ")')
    Write(unit_wfn,'("# ")')
    do K = 1, dim_K_out
      write(unit_wfn,*) indx_out(K), h_vec_out(1:nvib,K), coef_out(K)
    enddo
    !
    close(unit_wfn)
  end subroutine
  !
  !
END PROGRAM prunesymm
