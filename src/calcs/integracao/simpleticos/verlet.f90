! ************************************************************
!! METODO NUMERICO: Velocity-Verlet
!
! Objetivos:
!   Aplicacao do metodo simpletico de Velocity-Verlet.
!
! Modificado:
!   15 de marco de 2024 (criado)
!   18 de janeiro de 2026 (modificado)
!
! Autoria:
!   oap
!
MODULE verlet
  USE tipos
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_verlet

  TYPE, EXTENDS(integracao) :: integracao_verlet

    CONTAINS
      PROCEDURE :: metodo, metodo_mi

  END TYPE

CONTAINS

! ************************************************************
!! Metodo numerico
!
! Objetivos:
!   Aplicacao do metodo em si.
!
! Modificado:
!   18 de janeiro de 2026
!
! Autoria:
!   oap
!
SUBROUTINE metodo (self, R, P, FSomas)

  IMPLICIT NONE
  class(integracao_verlet), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas

  ! Integrando as posicoes
  R = R + self % h * (P + 0.5_pf * self%h * FSomas) * self%massasInvertidas

  ! Velocidades
  P = P + 0.5_pf * self % h * FSomas

  ! Calcula as novas forcas
  FSomas = self%forcas(R)

  ! Integrando as velocidades
  P = P + 0.5_pf * self % h * FSomas

END SUBROUTINE metodo

! ************************************************************
!! Metodo numerico
!
! Objetivos:
!   Aplicacao do metodo em si.
!
! Modificado:
!   18 de janeiro de 2026
!
! Autoria:
!   oap
!
SUBROUTINE metodo_mi (self, R, P, FSomas)

  IMPLICIT NONE
  class(integracao_verlet), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  REAL(pf) :: const
  
  const = 0.5_pf * self % h * self % m_esc ! Constante auxiliar

  ! Integrando as posicoes
  R = R + self % h * (P * self % m_inv + const * FSomas)

  ! Velocidades
  const = const * self % m_esc
  P = P + const*Fsomas

  ! Calcula as novas forcas
  FSomas = self%forcas(R)

  ! Integrando as velocidades
  P = P + const*Fsomas

END SUBROUTINE metodo_mi

END MODULE verlet