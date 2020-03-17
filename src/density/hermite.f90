subroutine hn_polynomial_value ( n, x, p_out )

!*****************************************************************************80
!
!! HN_POLYNOMIAL_VALUE evaluates Hn(n,x).
!
!    Hn(n,x) is the normalized physicist's Hermite polynomial of degree n.
!
!    These polynomials satisfy the orthonormality condition:
!
!      Integral ( -oo < X < +oo ) 
!        exp ( - X^2 ) * Hn(M,X) Hn(N,X) dX = delta ( N, M )
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
!    NIST Handbook of Mathematical Functions,
!    Cambridge University Press, 2010,
!    ISBN: 978-0521192255,
!    LC: QA331.N57.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order polynomial to compute.
!
!    Input, real ( kind = 8 ) X, the evaluation points.
!
!    Output, real ( kind = 8 ) P(N), the values of the polynomials in X  
!    of index N.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) fact
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(0:n)
  real ( kind = 8 ) p_out
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) two_power
  real ( kind = 8 ) x

  p(0) = 1.0D+00

  if ( n == 0 ) then
    p_out = p(n)
    return
  end if

  p(1) = 2.0D+00 * x
 
  do j = 2, n
    p(j) = 2.0D+00 * x * p(j-1) &
      - 2.0D+00 * real ( j - 1, kind = 8 ) * p(j-2)
  end do
!
!  Normalize.
!
  fact = 1.0D+00
  two_power = 1.0D+00
  do j = 0, n
    !p(j) = p(j) / sqrt ( fact * two_power * sqrt ( r8_pi ) ) ORIGINALE
    p(j) = p(j) / sqrt ( fact * two_power )
    fact = fact * real ( j + 1, kind = 8 )
    two_power = two_power * 2.0D+00
  end do

! marco: 
! in output only n^th component
  p_out = p(n)

  return
end subroutine
