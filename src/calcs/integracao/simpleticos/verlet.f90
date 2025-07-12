! ************************************************************
!! METODO NUMERICO: Velocity-Verlet
!
! Objetivos:
!   Aplicacao do metodo simpletico de Velocity-Verlet.
!
! Modificado:
!   15 de marco de 2024 (criado)
!   12 de julho de 2025 (modificado)
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
!   12 de julho de 2025
!
! Autoria:
!   oap
!
FUNCTION metodo (self, R, P, FSomas_ant)

  IMPLICIT NONE
  class(integracao_verlet), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, FSomas_prox
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo

  ! Integrando as posicoes
  R1 = (P + 0.5_pf * self%h * FSomas_ant) * self%massasInvertidas
  R1 = R + self % h * R1

  ! Calcula as novas forcas
  FSomas_prox = self%forcas(R1)

  ! Integrando as velocidades
  P1 = P + 0.5_pf*self%h*(FSomas_ant + FSomas_prox)

  metodo(1,:,:) = R1
  metodo(2,:,:) = P1
  metodo(3,:,:) = FSomas_prox

END FUNCTION metodo

! ************************************************************
!! Metodo numerico
!
! Objetivos:
!   Aplicacao do metodo em si.
!
! Modificado:
!   12 de julho de 2025
!
! Autoria:
!   oap
!
FUNCTION metodo_mi (self, R, P, FSomas_ant)

  IMPLICIT NONE
  class(integracao_verlet), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, FSomas_prox
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo_mi
  REAL(pf) :: const

  ! Integrando as posicoes
  const = 0.5_pf * self % h * self % m_esc
  R1 = P * self % m_inv + const * FSomas_ant
  R1 = R + self % h * R1

  ! Calcula as novas forcas
  FSomas_prox = self%forcas(R1)

  ! Integrando as velocidades
  const = const * self % m_esc
  P1 = P + const*(Fsomas_ant + FSomas_prox)

  metodo_mi(1,:,:) = R1
  metodo_mi(2,:,:) = P1
  metodo_mi(3,:,:) = FSomas_prox

END FUNCTION metodo_mi

END MODULE verlet