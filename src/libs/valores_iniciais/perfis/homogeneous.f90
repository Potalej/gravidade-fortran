! ************************************************************
!! Homogeneous Sphere profile
!
!> Objectives
!  This module implements functions relative to the Homogeneous
!  sphere profile. Because of the non continuity of the density
!  function the ergodic distribution function seems to be not 
!  physical.
!
!> Modified
!  2026.09.28
!
!> Created
!  2026.09.28
!
!> Author
!  oap
!
MODULE iv_homogeneous_sphere_mod
  IMPLICIT NONE
  PUBLIC
  PRIVATE PI, M0

  REAL(8), PARAMETER :: PI = 4.0d0 * ATAN(1.0d0)
  REAL(8), PARAMETER :: M0 = 1.0d0
CONTAINS

FUNCTION homogeneous_sphere_density (r, pars) RESULT(density)
! Homogeneous sphere mass density profile:
!   rho(r) = rho_0 1_[0,R](r),
! where
!   rho_0 = 3 / (4 pi R^3),
! in order to have M0 = 1.
  REAL(8), INTENT(IN) :: r, pars(2)
  REAL(8) :: density
  REAL(8) :: Rpar, rho0

  Rpar = pars(2)
  rho0 = 3.0d0 / (4.0d0 * PI * Rpar**3)

  IF (r > Rpar) THEN
    density = 0.0d0
  ELSE
    density = rho0
  ENDIF
END FUNCTION

FUNCTION homogeneous_sphere_inverse_P (x, pars) RESULT(invP)
! The cumulated density is P(r) = M(r) / M_0, this solves P^-1(X) = r
! In this case, the cumulative for the Homogeneous Sphere is
!   M(r) = r^3 / R^3 if r < R else 1
! so we have
!   r = R (x M_0)^3
  REAL(8), INTENT(IN) :: x, pars(2)
  REAL(8) :: invP
  REAL(8) :: Rpar

  Rpar = pars(2)

  invP = Rpar * (x*M0)**(1.0d0/3.0d0)
END FUNCTION

FUNCTION homogeneous_sphere_potential (r, pars) RESULT(potential)
! Homogeneous Sphere potential, given by
!   Phi(r) = G (r^2 - 3 R^2)/2
! if r < R, and
!   Phi(r) = G / r
! else.
! The function Phi solves the Poisson equation for the Hernquist
! density profile.
  REAL(8), INTENT(IN) :: r, pars(2)
  REAL(8) :: potential
  REAL(8) :: G, Rpar

  G = pars(1)
  Rpar = pars(2)

  IF (r > Rpar) THEN
    potential = G / r
  ELSE
    potential = 0.5d0 * G * (r**2 - 3.0d0 * Rpar**2)
  ENDIF
END FUNCTION

END MODULE