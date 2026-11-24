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
!  2026.11.24
!
!> Created
!  2026.09.28
!
!> Author
!  oap
!
MODULE iv_homogeneous_sphere_mod
  USE tipos
  IMPLICIT NONE
  PUBLIC
  PRIVATE PI, M0

  REAL(pf), PARAMETER :: PI = 4.0_pf * ATAN(1.0_pf)
  REAL(pf), PARAMETER :: M0 = 1.0_pf
CONTAINS

FUNCTION homogeneous_sphere_density (r, pars) RESULT(density)
! Homogeneous sphere mass density profile:
!   rho(r) = rho_0 1_[0,R](r),
! where
!   rho_0 = 3 / (4 pi R^3),
! in order to have M0 = 1.
  REAL(pf), INTENT(IN) :: r, pars(2)
  REAL(pf) :: density
  REAL(pf) :: Rpar, rho0

  Rpar = pars(2)
  rho0 = 3.0_pf / (4.0_pf * PI * Rpar**3)

  IF (r > Rpar) THEN
    density = 0.0_pf
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
  REAL(pf), INTENT(IN) :: x, pars(2)
  REAL(pf) :: invP
  REAL(pf) :: Rpar

  Rpar = pars(2)

  invP = Rpar * (x*M0)**(1.0_pf/3.0_pf)
END FUNCTION

FUNCTION homogeneous_sphere_potential (r, pars) RESULT(potential)
! Homogeneous Sphere potential, given by
!   Phi(r) = G (r^2 - 3 R^2)/2
! if r < R, and
!   Phi(r) = G / r
! else.
! The function Phi solves the Poisson equation for the Hernquist
! density profile.
  REAL(pf), INTENT(IN) :: r, pars(2)
  REAL(pf) :: potential
  REAL(pf) :: G, Rpar

  G = pars(1)
  Rpar = pars(2)

  IF (r > Rpar) THEN
    potential = G / r
  ELSE
    potential = 0.5_pf * G * (r**2 - 3.0_pf * Rpar**2)
  ENDIF
END FUNCTION

END MODULE