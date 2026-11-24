! ************************************************************
!! Plummer profile (Plummer 1911)
!
!> Objectives
!  This module implements functions relative to the Plummer
!  sphere profile, both the mass density and for isotropic
!  velocities.
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
MODULE iv_plummer_mod
  USE tipos
  IMPLICIT NONE
  PUBLIC
  PRIVATE PI, M0

  REAL(pf), PARAMETER :: PI = 4.0_pf * ATAN(1.0_pf)
  REAL(pf), PARAMETER :: M0 = 1.0_pf
CONTAINS

FUNCTION plummer_density (r, pars) RESULT(density)
! Plummer mass density profile:
! rho(r) = 3 M_0 b^2 / (4 pi (r^2 + b^2)^(5/2))
  REAL(pf), INTENT(IN) :: r, pars(2)
  REAL(pf) :: density
  REAL(pf) :: G, b

  G = pars(1)
  b = pars(2)

  density = 3.0_pf * M0 * b**2 / (4.0_pf * PI * SQRT(r**2 + b**2)**5.0_pf)
END FUNCTION

FUNCTION plummer_radius_pdf (r, pars) RESULT(pdf)
! Plummer pdf:
! dP / dr = 4 \pi r^2 \rho(r) / M_0
  REAL(pf), INTENT(IN) :: r, pars(2)
  REAL(pf) :: pdf
  REAL(pf) :: G, b

  G = pars(1)
  b = pars(2)

  pdf = 3.0_pf * b**2 * r**2 / (r**2 + b**2)**2.5_pf
END FUNCTION

FUNCTION plummer_inverse_P (x, pars) RESULT(invP)
! The cumulated density is P(r) = M(r) / M_0, this solves P^-1(X) = r
  REAL(pf), INTENT(IN) :: x, pars(2)
  REAL(pf) :: invP
  REAL(pf) :: b

  b = pars(2)

  invP = x**(1.0_pf/3.0_pf) * b / SQRT(1.0_pf - x**(2.0_pf/3.0_pf))
END FUNCTION

FUNCTION plummer_potential (r, pars) RESULT(potential)
! Plummer potential, given by
! Phi(r) = - G M_0 / (r^2 + b^2)^(1/2)
! The function Phi solves the Poisson equation for the Plummer
! density profile.
  REAL(pf), INTENT(IN) :: r, pars(2)
  REAL(pf) :: potential
  REAL(pf) :: G, b

  G = pars(1)
  b = pars(2)

  potential = - G * M0 / SQRT(r**2 + b**2)
END FUNCTION

FUNCTION plummer_dPhi (r, pars) RESULT(dPhi)
! Derivative of the Plummer potential wrt the radius
! dPhi / dt = G M_0 r / (r^2 + b^2)^(3/2) 
  REAL(pf), INTENT(IN) :: r, pars(2)
  REAL(pf) :: dPhi
  REAL(pf) :: G, b

  G = pars(1)
  b = pars(2)

  dPhi = G * M0 * r / SQRT(r**2 + b**2)**3
END FUNCTION

FUNCTION plummer_edf_norm (x) RESULT(edf)
! Ergodic distribution function without the constants.
  REAL(pf), INTENT(IN) :: x
  REAL(pf) :: edf

  edf = (1.0_pf - x**2)**3.5_pf
END FUNCTION

FUNCTION plummer_maximum_pdf () RESULT(maximum)
! Maximum value of the probability distribution function
  REAL(pf) :: maximum
  maximum = 0.0922108113291_pf
END FUNCTION

END MODULE