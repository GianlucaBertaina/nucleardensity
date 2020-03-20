! COPYRIGHT (C) 2020 Chiara Donatella Aieta, Marco Micciarelli, Gianluca Bertaina, Michele Ceotto
! See LICENSE for details

        PROGRAM crude
        !
        !-----------------!
        !  C-ustomized    !
        !  R-endering of  !
        !  U-nified       !
        !  DE-nsities     !
        !-----------------!
        ! 
        IMPLICIT NONE

        INTEGER :: i, j, k, ii
        INTEGER :: atom1, atom2, ix, iy, iz
        CHARACTER(len=500) :: arg, filecube1, filecube2, task
        REAL*8 :: coef1, coef2
        REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: density1, density2, &
                                                 density
        INTEGER :: nat, nxpoints, nypoints, nzpoints
        REAL*8 :: dx, dy, dz, zero_real, cubelmax
        REAL*8 :: dummy, rdx, rdy, rdz, xold, yold, zold, dnucleai
        REAL*8 :: xnew, ynew, znew, step, xbar, zbar, ybar
        REAL*8 :: x0, y0, z0, dp, cutdens
        INTEGER, ALLOCATABLE, DIMENSION(:) :: Z_nuclei
        REAL*8, ALLOCATABLE, DIMENSION(:) :: x_eq_cart

!       Input selection       
        i=0
        CALL get_command_argument(i, arg)
        i=1
        CALL get_command_argument(i, arg)
        task = TRIM(arg)

        IF (task=="comb") THEN
 
           DO i=2, 5
   
             CALL get_command_argument(i, arg)
             
             SELECT CASE (i)
                CASE(2) 
                    filecube1 = TRIM(arg)           
                CASE(3)             
                    filecube2 = TRIM(arg)           
                CASE(4)             
                    READ(arg,*) coef1
                CASE(5)
                    READ(arg,*) coef2
                CASE DEFAULT
                    WRITE(*,*) "Ignored argument: ", arg
             END SELECT             
   
           END DO
 
        WRITE(*,*) "Input par: ", trim(filecube1), trim(filecube2), coef1, coef2

        ELSEIF (task=="cut") THEN

           DO i=2, 4
   
             CALL get_command_argument(i, arg)
             
             SELECT CASE (i)
                CASE(2) 
                    filecube1 = TRIM(arg)           
                CASE(3)             
                    READ(arg,*) atom1
                CASE(4)
                    READ(arg,*) atom2
                CASE DEFAULT
                    WRITE(*,*) "Ignored argument: ", arg
             END SELECT             
   
            END DO          
 
          WRITE(*,*) "Input par: ", filecube1, atom1, atom2

        ELSEIF (task=="--help") THEN

          WRITE(*,*)
          WRITE(*,*)  "CRUDE - Customized Rendering of Unified &
                                                              DEnsities"
          WRITE(*,*)
          WRITE(*,*) "Calling instructions:"
          WRITE(*,*)
          WRITE(*,*) "=> To combine two cube files" 
          WRITE(*,*) "   ./crude.x comb <file1.cube> <file2.cube> &
                                                                <a> <b>"
          WRITE(*,*) "   <file1.cube>  : 1st cube file name"
          WRITE(*,*) "   <file2.cube>  : 2nd cube file name"
          WRITE(*,*) "   <a>           : 1st linear comb. coefficient"
          WRITE(*,*) "   <b>           : 2nd linear comb. coefficient"
          WRITE(*,*) "   Output: a * file1.cube + b * file2.cube"
          WRITE(*,*)
          WRITE(*,*) "=> To generate cube file cut along a two &
                                                              atom axis"
          WRITE(*,*) "   ./crude.x cut <file.cube> <a> <b>"
          WRITE(*,*) "   <file.cube>  : cube file name"
          WRITE(*,*) "   <a>          &
                         : 1st atom index (as listed in xyz file)"
          WRITE(*,*) "   <b>          &
                         : 2nd atom index (as listed in xyz file)"
  
          STOP

        ELSE
 
            WRITE(*,*) "Input command ", arg, "is unknown."
            WRITE(*,*) "Try: ./crude --help"
            WRITE(*,*) "CRUDE stops."
            STOP

        END IF

        IF (task=="comb") THEN
        
           !associate the first cube file to unit 30
           OPEN (UNIT=30, FILE=filecube1, STATUS='old')
   
           ! write out cube file "introdction"
           ! Write somethig useful in comment lines TBD
           READ(30,*) !"!comment"
           READ(30,*) !"!comment"
           READ(30,*) nat, cubelmax,cubelmax,cubelmax
           READ(30,*) nxpoints, dx, zero_real, zero_real
           READ(30,*) nypoints, zero_real, dy, zero_real
           READ(30,*) nzpoints, zero_real, zero_real, dz
   
           ALLOCATE(Z_nuclei(nat), x_eq_cart(3*nat))
           ALLOCATE(density1(nxpoints,nypoints,nzpoints))
           ALLOCATE(density2(nxpoints,nypoints,nzpoints))
   
   
           !Print out the equilibrium geometry of the molecule
           DO i=1, nat
              READ(30,*) Z_nuclei(i), dummy, &
                          x_eq_cart(3*i-2),&
                          x_eq_cart(3*i-1),&
                          x_eq_cart(3*i)  
           END DO
   
           DO i=1, nxpoints
             DO j=1, nypoints
                   READ(30,'( 6F12.6 )') (density1(i, j, k), k=1, nzpoints)
             END DO
           END DO
   
           CLOSE(30)
   
           !associate the second cube file to unit 31
           OPEN (UNIT=31, FILE=filecube2, STATUS='old')
   
           ! write out cube file "introdction"
           ! Write somethig useful in comment lines TBD
           READ(31,*) !"!comment"
           READ(31,*) !"!comment"
           READ(31,*) nat!, -cubelmax,&
                      !             -cubelmax,&
                      !             -cubelmax
           READ(31,*) nxpoints!, dx, zero_real, zero_real
           READ(31,*) nypoints!, zero_real, dy, zero_real
           READ(31,*) nzpoints!, zero_real, zero_real, dz
   
   
           !Print out the equilibrium geometry of the molecule
           DO i=1, nat
              READ(31,*) !Z_nuclei(i), real(Z_nuclei(i),8), &
                          !x_eq_cart(3*i-2),&
                          !x_eq_cart(3*i-1),&
                          !x_eq_cart(3*i)  
           END DO
   
           DO i=1, nxpoints
             DO j=1, nypoints
                   READ(31,'( 6F12.6 )') (density2(i, j, k), k=1, nzpoints)
             END DO
           END DO
   
           CLOSE(31)
   
           OPEN (UNIT=32, FILE='cuberes.cube', STATUS='replace')
   
           ! write out cube file "introdction"
           ! Write somethig useful in comment lines TBD
           WRITE(32,*) "Combined cubefile from input files:"
           WRITE(32,*) coef1, trim(filecube1), coef2, trim(filecube2)
           WRITE(32,*) nat, cubelmax,&
                                   cubelmax,&
                                   cubelmax
           WRITE(32,*) nxpoints, dx, zero_real, zero_real
           WRITE(32,*) nypoints, zero_real, dy, zero_real
           WRITE(32,*) nzpoints, zero_real, zero_real, dz
   
   
           !Print out the equilibrium geometry of the molecule
           DO i=1, nat
              WRITE(32,*) Z_nuclei(i), real(Z_nuclei(i),8), &
                          x_eq_cart(3*i-2),&
                          x_eq_cart(3*i-1),&
                          x_eq_cart(3*i)  
           END DO
   
           DO i=1, nxpoints
             DO j=1, nypoints
                WRITE(32,'( 6F12.6 )')& 
                     (coef1*density1(i, j, k) +& 
                      coef2*density2(i, j, k), k=1, nzpoints)
             END DO
           END DO
   
           CLOSE(32)

       ELSE IF (task == "cut") THEN
        
           OPEN (UNIT=33, FILE=filecube1, STATUS='old')
           OPEN (UNIT=34, FILE="denscut.dat", STATUS='replace')

           !Read the cubefile
           READ(33,*) !"!comment"
           READ(33,*) !"!comment"
           READ(33,*) nat, cubelmax,&
                                   cubelmax,&
                                   cubelmax
           READ(33,*) nxpoints, dx, zero_real, zero_real
           READ(33,*) nypoints, zero_real, dy, zero_real
           READ(33,*) nzpoints, zero_real, zero_real, dz

           ALLOCATE(Z_nuclei(nat), x_eq_cart(3*nat))
           ALLOCATE(density(nxpoints,nypoints,nzpoints))

           ! x_eq_cart components are expressed in Bohr
           DO i=1, nat
              READ(33,*) Z_nuclei(i), dummy, &
                          x_eq_cart(3*i-2),&
                          x_eq_cart(3*i-1),&
                          x_eq_cart(3*i)  
           END DO

           DO i=1, nxpoints
             DO j=1, nypoints
                   READ(33,'( 6F12.6 )') (density(i, j, k), k=1,& 
                                                          nzpoints)
             END DO
           END DO

           dp = SQRT(dx**2 + dy**2 + dz**2)
           WRITE(*,*) "dp ", dp
 
           !vec = rdx i + rdy j + rdz k
           rdx = x_eq_cart(3*atom2-2) - x_eq_cart(3*atom1-2)
           rdy = x_eq_cart(3*atom2-1) - x_eq_cart(3*atom1-1)
           rdz = x_eq_cart(3*atom2) - x_eq_cart(3*atom1)

           WRITE(*,*) "Bond length", SQRT(rdx**2 + rdy**2 + rdz**2)

           x0 = x_eq_cart(3*atom1-2) 
           y0 = x_eq_cart(3*atom1-1)
           z0 = x_eq_cart(3*atom1)

           DO i=1, nat
              WRITE(36,*) x_eq_cart(3*i-2),&
                          x_eq_cart(3*i-1),&
                          x_eq_cart(3*i)  
           END DO

           rdx = rdx / 101
           rdy = rdy / 101
           rdz = rdz / 101

           step = SQRT(rdx**2 + rdy**2 + rdz**2)

           WRITE(*,*) "step", step

!             first part before atom1
              xold = x0 - rdx*401
              yold = y0 - rdy*401
              zold = z0 - rdz*401

              DO i = 1, 400

                 xnew = xold + rdx
                 ynew = yold + rdy
                 znew = zold + rdz
                 WRITE(35,*) xnew, ynew, znew

                 CALL find_cube_index(nxpoints, nypoints, nzpoints,&
                                      dx, dy, dz,&
                                      xnew, ynew, znew,& 
                                      ix, iy, iz)

                 IF( (ix+1 <= nxpoints .and. ix-1 > 0) .and. & 
                     (iy+1 <= nypoints .and. iy-1 > 0) .and. & 
                     (iz+1 <= nzpoints .and. iz-1 > 0) ) THEN
                 
                 cutdens =   density(ix, iy, iz) +&
                             density(ix+1, iy, iz) +&
                             density(ix, iy+1, iz) +&
                             density(ix, iy, iz+1) +&
                             density(ix-1, iy, iz) +&
                             density(ix, iy-1, iz) +&
                             density(ix, iy, iz-1) +&
                             density(ix+1, iy+1, iz) +&
                             density(ix, iy+1, iz+1) +&
                             density(ix+1, iy, iz+1) +&
                             density(ix-1, iy+1, iz) +&
                             density(ix+1, iy-1, iz) +&
                             density(ix, iy-1, iz+1) +&
                             density(ix, iy+1, iz-1) +&
                             density(ix-1, iy, iz+1) +&
                             density(ix+1, iy, iz-1) +&
                             density(ix-1, iy-1, iz) +&
                             density(ix, iy-1, iz-1) +&
                             density(ix-1, iy, iz-1) +&
                             density(ix+1, iy+1, iz+1) +&
                             density(ix-1, iy-1, iz-1) +&
                             density(ix+1, iy-1, iz-1) +&
                             density(ix-1, iy+1, iz+1) +&
                             density(ix-1, iy+1, iz-1) +&
                             density(ix+1, iy-1, iz+1) +&
                             density(ix-1, iy-1, iz+1) +&
                             density(ix+1, iy+1, iz-1) 
                             
                 cutdens = cutdens / 27.d0

                 WRITE(34,*) (-401+i)*step, cutdens

                 ELSE
                
                 WRITE(34,*) (-401+i)*step, 0.d0

                 END IF
 
                 xold = xnew
                 yold = ynew
                 zold = znew

              END DO

!             starting from 0
              xold = x0
              yold = y0
              zold = z0

              DO i = 1, 400

                 xnew = xold + rdx
                 ynew = yold + rdy
                 znew = zold + rdz
                 WRITE(35,*) xnew, ynew, znew

                 CALL find_cube_index(nxpoints, nypoints, nzpoints,&
                                      dx, dy, dz,&
                                      xnew, ynew, znew,& 
                                      ix, iy, iz)

                 IF( (ix+1 <= nxpoints .and. ix-1 > 0) .and. & 
                     (iy+1 <= nypoints .and. iy-1 > 0) .and. & 
                     (iz+1 <= nzpoints .and. iz-1 > 0) ) THEN

                 cutdens =   density(ix, iy, iz) +&
                             density(ix+1, iy, iz) +&
                             density(ix, iy+1, iz) +&
                             density(ix, iy, iz+1) +&
                             density(ix-1, iy, iz) +&
                             density(ix, iy-1, iz) +&
                             density(ix, iy, iz-1) +&
                             density(ix+1, iy+1, iz) +&
                             density(ix, iy+1, iz+1) +&
                             density(ix+1, iy, iz+1) +&
                             density(ix-1, iy+1, iz) +&
                             density(ix+1, iy-1, iz) +&
                             density(ix, iy-1, iz+1) +&
                             density(ix, iy+1, iz-1) +&
                             density(ix-1, iy, iz+1) +&
                             density(ix+1, iy, iz-1) +&
                             density(ix-1, iy-1, iz) +&
                             density(ix, iy-1, iz-1) +&
                             density(ix-1, iy, iz-1) +&
                             density(ix+1, iy+1, iz+1) +&
                             density(ix-1, iy-1, iz-1) +&
                             density(ix+1, iy-1, iz-1) +&
                             density(ix-1, iy+1, iz+1) +&
                             density(ix-1, iy+1, iz-1) +&
                             density(ix+1, iy-1, iz+1) +&
                             density(ix-1, iy-1, iz+1) +&
                             density(ix+1, iy+1, iz-1) 

                 cutdens = cutdens / 27.d0
                 
                 WRITE(34,*) i*step, cutdens

                 ELSE
                
                 WRITE(34,*) i*step, 0.d0

                 END IF

                 xold = xnew
                 yold = ynew
                 zold = znew

              END DO

       END IF

        END PROGRAM crude

!==================================================================
!==================================================================
!==================================================================

      SUBROUTINE find_cube_index(nxpoints, nypoints, nzpoints,&
                                 dx, dy, dz,&
                                 x, y, z,& 
                                 ix, iy, iz)

        IMPLICIT NONE
        INTEGER, INTENT(OUT) :: ix, iy, iz
        INTEGER, INTENT(IN) :: nxpoints, nypoints, nzpoints
        REAL*8, INTENT(IN) :: x, y, z, dx, dy, dz

           ix = NINT(x / dx)
           ix = (nxpoints + 1) / 2 + ix
           IF (ix > nxpoints) ix = ix - 1
           IF (ix < 1) ix = ix + 1

           iy = NINT(y / dy)
           iy = (nypoints + 1) / 2 + iy
           IF (iy > nypoints) iy = iy - 1
           IF (iy < 1) iy = iy + 1

           iz = NINT(z / dz)
           iz = (nzpoints + 1) / 2 + iz
           IF (iz > nzpoints) iz = iz - 1
           IF (iz < 1) iz = iz + 1

      END SUBROUTINE find_cube_index

