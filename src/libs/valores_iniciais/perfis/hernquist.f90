! ************************************************************
!! Hernquist profile (Hernquist 1990)
!
!> Objectives
!  This module implements functions relative to the Hernquist
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
MODULE iv_hernquist_mod
  USE tipos
  IMPLICIT NONE
  PUBLIC
  PRIVATE PI, M0

  REAL(pf), PARAMETER :: PI = 4.0_pf * ATAN(1.0_pf)
  REAL(pf), PARAMETER :: M0 = 1.0_pf
CONTAINS

FUNCTION hernquist_density (r, pars) RESULT(density)
! Hernquist mass density profile:
!   rho(r) = M_0 a / (2 pi r (r+a)^3)
  REAL(pf), INTENT(IN) :: r, pars(2)
  REAL(pf) :: density
  REAL(pf) :: G, a

  G = pars(1)
  a = pars(2)

  density = M0 * a / (2.0_pf * PI * r * (r + a)**3.0_pf)
END FUNCTION

FUNCTION hernquist_radius_pdf (r, pars) RESULT(pdf)
! Hernquist pdf:
! dP / dr = 4 \pi r^2 \rho(r) / M_0
  REAL(pf), INTENT(IN) :: r, pars(2)
  REAL(pf) :: pdf
  REAL(pf) :: G, a

  G = pars(1)
  a = pars(2)

  pdf = 2.0_pf * r * a / (r + a)**3.0_pf
END FUNCTION

FUNCTION hernquist_inverse_P (x, pars) RESULT(invP)
! The cumulated density is P(r) = M(r) / M_0, this solves P^-1(X) = r
! In this case, the cumulative for the Hernquist is
!   M(r) = M_0 r^2 / (a + r)^2
! so we have
!   r = a SQRT(x) / (1 - SQRT(x))
  REAL(pf), INTENT(IN) :: x, pars(2)
  REAL(pf) :: invP
  REAL(pf) :: a

  a = pars(2)

  invP = a * SQRT(x) / (1.0_pf - SQRT(x))
END FUNCTION

FUNCTION hernquist_potential (r, pars) RESULT(potential)
! Hernquist potential, given by
!   Phi(r) = - G M_0 / (a + r)
! The function Phi solves the Poisson equation for the Hernquist
! density profile.
  REAL(pf), INTENT(IN) :: r, pars(2)
  REAL(pf) :: potential
  REAL(pf) :: G, a

  G = pars(1)
  a = pars(2)

  potential = - G * M0 / (a + r)
END FUNCTION

FUNCTION hernquist_dPhi (r, pars) RESULT(dphi)
! Hernquist potential derivative
  REAL(pf), INTENT(IN) :: r, pars(2)
  REAL(pf) :: dphi
  REAL(pf) :: G, a

  G = pars(1)
  a = pars(2)

  dphi = G * M0 / (a + r)**2.0_pf
END FUNCTION

FUNCTION hernquist_edf_norm (r, x, pars) RESULT(edf)
! Ergodic distribution function without the constants.
  REAL(pf), INTENT(IN) :: r, x, pars(2)
  REAL(pf) :: edf, a, E
  a = pars(2)

  E = SQRT(-hernquist_potential(r, pars) * (1.0_pf - x**2))

  edf = 1.0_pf / (8.0_pf * SQRT(2.0_pf) * a * PI**3)
  edf = edf / SQRT(1.0_pf - E)**5.0_pf
  edf = edf * 3.0_pf * ASIN(SQRT(E))
  edf = edf + edf * SQRT(E*(1.0_pf - E))*(1.0_pf-2.0_pf*E)*(8.0_pf*E**2-8.0_pf*E-3.0_pf)
END FUNCTION

FUNCTION hernquist_pdf_derivative (r, x, h, pars) RESULT(df)
! derivative of the probability function p(v|r) \propto v^2 f(E)
! for this, we use the centered finite differences
  REAL(pf), INTENT(IN) :: r, x, h, pars(2)
  REAL(pf) :: df

  df = (x+h)**2 * hernquist_edf_norm(r, x+h, pars) 
  df = df - (x-h)**2 * hernquist_edf_norm(r, x-h, pars)
  df = df / (2.0_pf * h)
END FUNCTION

FUNCTION hernquist_maximum_pdf (r, pars) RESULT(maximum)
! Maximum value of the probability distribution function
! For this, we look for the root of g'(q) using finite differences
  REAL(pf), INTENT(IN) :: r, pars(2)
  REAL(pf) :: x, h, der_befo, der_next, maximum
  
  h = 0.001_pf
  x = h
  der_befo = hernquist_pdf_derivative(r, x, h, pars)
  maximum = 1.0_pf
  DO WHILE (x + h < 1.0_pf)
    x = x + h
    der_next = hernquist_pdf_derivative(r, x, h, pars)

    ! if the signal changes, we crossed g'(q)=0
    IF (der_befo * der_next <= 0) THEN
      maximum = x*x*hernquist_edf_norm(r, x, pars)
      EXIT
    ENDIF

    der_befo = der_next
  END DO
END FUNCTION

END MODULE