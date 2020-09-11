PROGRAM gramschmidt
  implicit none

  integer, parameter   :: unit_input=10,unit_wfn=20
  character (len=80),allocatable :: file_wfn_in(:), file_wfn_out(:)

  integer              :: nvib,nat,ncart,nrotrasl,i

  integer              :: dim_K, nwf
  integer, allocatable :: h_vec(:,:)
  real*8,  allocatable :: coef(:,:)
  real*8,  allocatable :: indx(:)
  real*8               :: norm

  character(30),parameter    :: inputfile='input_gramschmidt.dat'
  character(len=7),parameter :: app="_gs.dat"

  call read_input
  !
  call get_wfn
  !
  write(*,*) "Initial overlap matrix"
  call calc_overlaps
  !
  call gs
  !
  write(*,*) "Overlap matrix after Gram Schmidt"
  call calc_overlaps
  !
  call print_wavefunction

  !
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_input
    implicit none

    logical :: file_exists=.true.
    integer :: i
    integer :: ppos

    !! check input file existence
    inquire(file=trim(inputfile), exist=file_exists)
    if (.not.file_exists) then
      print*, "Missing gramschmidt input file. Stopping program."
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
    read(unit_input,*) nwf ! Number of wavefunctions
    allocate(file_wfn_in(1:nwf),file_wfn_out(1:nwf))
    !
    read(unit_input,*) ! Here a comment line
    do i = 1, nwf
      read(unit_input,*) file_wfn_in(i) ! file with SC wavefunction
      ! Check all other input files existence
      inquire(file=trim(file_wfn_in(i)), exist=file_exists)
      if (.not.file_exists) then
        print*, "Missing file ",file_wfn_in(i) ,". Stopping program"
        STOP
      endif
      ppos = scan(trim(file_wfn_in(i)),".", BACK= .true.)
      if ( ppos > 0 ) file_wfn_out(i) = file_wfn_in(i)(1:ppos-1)//app
    enddo
    !
    close(unit_input)
    !
  end subroutine read_input
  !
  !
  subroutine get_wfn
    implicit none
    !
    integer :: K,i
    !
    allocate(h_vec(nvib,dim_K))
    allocate(coef(dim_K,1:nwf))
    allocate(indx(dim_K))
    !
    do i = 1,nwf
      open(unit_wfn,file=trim(file_wfn_in(i)),status='old', action="read")
      !
      read(unit_wfn,*) ! seven lines for comments
      read(unit_wfn,*)
      read(unit_wfn,*)
      read(unit_wfn,*)
      read(unit_wfn,*)
      read(unit_wfn,*)
      read(unit_wfn,*)
      do K = 1, dim_K
        read(unit_wfn,*) indx(K), h_vec(1:nvib,K), coef(K,i)
      enddo
      !
      close(unit_wfn)
    enddo
  end subroutine
  !
  !
  subroutine calc_overlaps()
    implicit none
    integer           :: i,j
    character(len=80) :: form
    real*8            :: overlaps(1:nwf,1:nwf)
    !
    write(form,'("(",1I4,"(1E13.6,2x))")') nwf
    !
    do i = 1,nwf
      do j = 1,nwf
        overlaps(i,j) = dot_product(coef(:,i),coef(:,j))
      enddo
    enddo
    write(*,form) overlaps
    !
  end subroutine
  !
  !
  subroutine gs()
    implicit none
    integer :: i,j
    real*8  :: overlap,norm
    !
    do i = 1,nwf-1
      ! Normalize coefficients of state to be removed
      norm = sqrt(sum(coef(:,i)**2))
      coef(:,i) = coef(:,i) / norm

      do j=i+1,nwf
        ! projection of j onto i
        ! Notice real coefficients are assumed
        overlap = dot_product(coef(:,i),coef(:,j))
        ! remove projection of i from j
        coef(:,j) = coef(:,j) - overlap * coef(:,i)
      enddo
    enddo
    ! Normalize last state
    norm = sqrt(sum(coef(:,nwf)**2))
    coef(:,nwf) = coef(:,nwf) / norm
    !
  end subroutine
  !
  !
  subroutine print_wavefunction
    implicit none
    !
    integer           :: K,i
    character(len=80) :: form
    !
    write(form,'("(1F8.2,1x,",1I4,"(I3,1x),1x,1E17.10)")') nvib

    do i = 1,nwf
      open(unit_wfn,file=trim(file_wfn_out(i)),status='unknown', action="write")
      !
      Write(unit_wfn,'("# Gram Schmidt")')
      Write(unit_wfn,'("# Number of basis elements: ",1I8)') dim_K
      Write(unit_wfn,'("# ")')
      Write(unit_wfn,'("# ")')
      Write(unit_wfn,'("# ")')
      Write(unit_wfn,'("# ")')
      Write(unit_wfn,'("# ")')
      do K = 1, dim_K
        write(unit_wfn,form) indx(K), h_vec(1:nvib,K), coef(K,i)
      enddo
      !
      close(unit_wfn)
    enddo
  end subroutine
  !
  !
END PROGRAM gramschmidt
