! ************************************************************
!! METODO NUMERICO: Euler explicito
!
! Objetivos:
!   Aplicacao do metodo de Euler explicito.
!
! Modificado:
!   14 de setembro de 2024 (criado)
!   18 de janeiro de 2026 (modificado)
!
! Autoria:
!   oap
!
MODULE euler_exp
  USE tipos
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_euler_exp

  TYPE, EXTENDS(integracao) :: integracao_euler_exp

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
  class(integracao_euler_exp), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  INTEGER :: a

  ! Integrando as posicoes
  DO a = 1, self % N
    R(a,:) = R(a,:) + (self % h / self % m(a)) * P(a,:)
  END DO

  ! Integrando as velocidades
  P = P + self%h * FSomas

  ! Atualizando as forcas
  FSomas = self % forcas(R)

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
  class(integracao_euler_exp), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas

  ! Integrando as posicoes
  R = R + self % h * self % m_inv * P

  ! Integrando as velocidades
  P = P + self%h * (self%m2 * FSomas)

  ! Atualizando as forcas
  FSomas = self%forcas(R)

END SUBROUTINE metodo_mi

END MODULE euler_exp