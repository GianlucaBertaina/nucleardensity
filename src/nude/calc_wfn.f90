      subroutine bar_wfn_squared(dq,wfn_b2)

      use input_def,  only: dim_K
      use input_def,  only: h_vec, coef
      use input_def,  only: ncart

      implicit none 

      real*8, dimension(ncart), intent(in)  :: dq
      real*8, intent(out)                   :: wfn_b2

      integer :: k
      real*8  :: psi

      wfn_b2 = 0.d0

      do k = 1, dim_K
	  call harm_wfn_bar(dq,h_vec(:,k),psi)
          wfn_b2 = wfn_b2 + coef(k) * psi
      enddo


      wfn_b2 = wfn_b2 **2


      end subroutine bar_wfn_squared


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine harm_wfn_bar(dq,vec,wfn)

      use input_def,  only: omega_vh, ncart,nvib


      implicit none 

      real*8, dimension(ncart), intent(in) :: dq
      integer, dimension(nvib ), intent(in) :: vec
      real*8, intent(out) :: wfn



      real*8  :: zz, norm, pol, fharm
      integer :: i
      


      wfn = 1.d0
      do i = 1, nvib
          zz = dq(i) * dsqrt(omega_vh(i))
          norm = 1.d0 ! hermite pols already normalized in routine
	  call hn_polynomial_value ( vec(i), zz, pol)
          fharm = pol !* EXP(-zz**2/2.d0)
          wfn = wfn * fharm * norm
      enddo

      !IF (ABS(wfn) < 1d-24) wfn = 0.d0


      end subroutine harm_wfn_bar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
