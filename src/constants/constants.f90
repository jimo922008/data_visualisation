module constants

  implicit none

  integer, public, parameter       :: dp=8 ! Single = 4 Double = 8 Quad = 16

  !integer, public, parameter       :: sp=4 (dgemm to sgmm)
  ! (apple silicon )

  integer, public, parameter       :: stdin=5
  integer, public, parameter       :: stdout=6
  integer, public, parameter       :: stderr=0

  real(kind=dp), public, parameter :: pi=3.141592653589793238462643383279502884197_dp
  real(kind=dp), public, parameter :: tpi=2.0_dp*pi
  real(kind=dp), public, parameter :: gr=(sqrt(5.0_dp)+1.0_dp)/2.0_dp
  real(kind=dp), public, parameter :: dgrd = pi/180.0_dp
  real(kind=dp), public, parameter :: evbyang3=160.2176487_dp
  real(kind=dp), public, parameter :: bohr2ang = 0.529177210903_dp
  real(kind=dp), public, parameter :: delta = 1e-13_dp

  real(kind=dp), public, parameter, dimension(3,3) :: &
       ident=reshape((/1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp/),(/3,3/))

end module constants