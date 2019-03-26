   PROGRAM NUDE

      include "mpif.h"
      use constants
      use io_units_def
      use input_def
      use manage_cubes


      implicit none

      REAL*8, ALLOCATABLE, DIMENSION(:)       :: q_eq, sigma_gaus 
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: qq, xx 
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: q_expect_value, q_expect_value_h
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: x_expect_value, x_expect_value_h
      REAL*8, ALLOCATABLE, DIMENSION(:)       :: norm 
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: density, density_sq, density_errorbar 


      integer :: i, ix,iy,iz, istep, ii, j, k
      real*8, dimension(3) :: Ri 
      real*8 :: bar_wfn_sq, tot_int, tot_int_sq
      !
      real*8  :: det_omega, norm_gaussian_normal
      !
      integer :: idum
      real*8  :: ran2
      real*8  :: xrand1,xrand2,y1,y2
      real    :: start, finish
      !
      idum = -200
      !
      !
      call cpu_time(start)
      print*, 'Entering program:  NUclear DEnsity:'
      print*, '  '
      print*, '  '
      !
      !
      print*, 'setting parameters from input files'
      call get_input_params
      print*,'done'
      print*, '  '
      print*, '  '
      !
      !
      print*, 'setting grid and parameters for cube files'
      call set_cube
      print*,'done'
      print*, '  '
      print*, '  '
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
      ! Initialize density (a 3D array for each atom)
      allocate(density(nxpoints, nypoints, nzpoints, nat))
      density = 0.d0
      allocate(density_sq(nxpoints, nypoints, nzpoints, nat))
      density_sq = 0.d0
      allocate(density_errorbar(nxpoints, nypoints, nzpoints, nat))
      density_errorbar = 0.d0
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
      !
      ! allocate corresponding structure in cartesian coordinates
      allocate(xx(ncart))
      allocate(x_expect_value(ncart)) 
      allocate(x_expect_value_h(ncart))
      !
      ! open unit for fil with MC traj in xyz
      open(unit_trajMC, file='traj.xyz')
      !
      ! MC SAMPLING CYCLE
      tot_int=0.d0
      tot_int_sq=0.d0
      DO istep=1,Nsteps_MC

          IF ( MOD(istep,1000)==0 ) THEN
             WRITE(*,*) istep, ' steps of ', Nsteps_MC, ' done!' 
          END IF
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

          tot_int = tot_int + bar_wfn_sq
          tot_int_sq = tot_int_sq + bar_wfn_sq**2
          !
	  ! Bring back structure in cartesian coordinates
          call to_cartesian(xx,qq)
          !
	  ! Print out structure in xyz format on file
          IF (Nsteps_MC <= 100000) THEN
           
	  write(unit_trajMC,* ) nat
	  write(unit_trajMC,* ) 
	  do i = 1, nat
	      write(unit_trajMC,* ) symb_list(i), xx(3*i-2:3*i) * FROMauTOang
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
	      density_sq(ix,iy,iz,i) = density_sq(ix,iy,iz,i) + bar_wfn_sq**2
	      !
	  enddo
          !
          !
      ENDDO
      !
      ! close unit for fil with MC traj in xyz
      close(unit_trajMC)

      !
      ! Normalize integral for expectation structures:
      q_expect_value   = q_expect_value   /tot_int
      q_expect_value_h = q_expect_value_h /Nsteps_MC
      !
      q_expect_value(4:9) = q_eq(4:9)
      q_expect_value_h(4:9) = q_eq(4:9)
      print*, 'q_eq    = ', q_eq
      print*, 'q_exp_h = ', q_expect_value_h
      print*, 'q_exp   = ', q_expect_value

      ! Bring back structures in cartesian coordinates
      call to_cartesian(x_expect_value,  q_expect_value)
      call to_cartesian(x_expect_value_h,q_expect_value_h)
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
      !
      !
      ! Normalize density:
      print*, 'total integ: ', tot_int / Nsteps_MC
      !
      det_omega = 1.d0
      do i = 1, nvib
         det_omega = det_omega * omega_vh(i)
      enddo
      !
      norm_gaussian_normal = DSQRT( (2*pi)**nvib * 0.5d0 / det_omega )

      !print*, 'norm gauss',  norm_gaussian_normal
      !
      !density = density / Nsteps_MC
      !CHIARA: density della riga sopra diviso per total integ 
      density = density / Nsteps_MC  
      !density = density / tot_int 
      !!!!density_sq = density_sq / tot_int_sq
      density_sq = density_sq / Nsteps_MC
      density_errorbar = SQRT( (density_sq - density**2) / Nsteps_MC )   
        
      !density = density * (norm_gaussian_normal / Nsteps_MC)
      !
      !
      !
      ! Print out nuclear densities on output files **SPOSTATO SOTTO**
      !print*, 'printing cube files'
      !call print_cube(density)
      
      !Check the normalization of each nucleus density
      ALLOCATE(norm(nat))
      norm = 0.d0
      DO ii=1, nat
         DO i=1, nxpoints
            DO j=1, nypoints
               DO k=1, nxpoints
                  norm(ii) = norm(ii) + density(i,j,k,ii)
               END DO
            END DO
          END DO
          print*, ii, ' nucleus norm: ', norm(ii) 
      END DO
      !
      !
      !
      !Normalize by the voxel dimension
      density = density / (dx*dy*dz) 
      !
      !
      !
      ! Print out nuclear densities on output files
      print*, 'printing cube files'
      call print_cube(density)
      !
      !
      !
      print*, 'printing cube files with density error bar'
      call print_cube_errorbar(density_errorbar)
      !
      !
      !
      print*,'done'
      print*, '  '
      print*, '  '
      !
      !
      !
      call cpu_time(finish)
      print '("Time for execution = ",f6.3," seconds.")',finish-start
      print*, '  '
      !
      !
      print*, 'ciao'



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
           REAL,PARAMETER :: pi = 3.1415926536
      
           y1 = sqrt(-2.0*log(x1)) * cos( 2.0*pi*x2 )
           y2 = sqrt(-2.0*log(x1)) * sin( 2.0*pi*x2 )
      
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


















