! ************************************************************
!! METODO NUMERICO: Euler Implicito
!
! Objetivos:
!   Aplicacao do metodo de Euler Implicito.
!
! Modificado:
!   14 de setembro de 2024 (criado)
!   18 de janeiro de 2026 (modificado)
!
! Autoria:
!   oap
!
MODULE euler_imp
  USE tipos
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_euler_imp

  TYPE, EXTENDS(integracao) :: integracao_euler_imp

    CONTAINS
      PROCEDURE :: metodo, metodo_mi

  END TYPE
  
CONTAINS

! ************************************************************
!! Metodo numerico
!
! Objetivos:
!   Aplicacao do metodo em si usando o metodo de Picard.
!
! Modificado:
!   18 de janeiro de 2026
!
! Autoria:
!   oap
!
SUBROUTINE metodo (self, R, P, FSomas)

  IMPLICIT NONE
  class(integracao_euler_imp), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, Rchute
  INTEGER :: i, a

  R1 = R
  P1 = P

  DO i = 1, 100
    ! Integrando as posicoes
    DO a = 1, self % N
      Rchute(a,:) = R(a,:) + self % h * P1(a,:) / self % m(a)
    END DO

    ! Integrando as velocidades
    P1 = P + self%h * FSomas

    ! Calcula as novas forcas
    FSomas = self%forcas(Rchute)

    IF (NORM2(R1 - Rchute) < self % h * self % h) EXIT
    R1 = Rchute
  END DO

  R = Rchute
  P = P1

END SUBROUTINE metodo

! ************************************************************
!! Metodo numerico (massas iguais)
!
! Objetivos:
!   Aplicacao do metodo em si usando o metodo de Picard.
!
! Modificado:
!   18 de janeiro de 2026
!
! Autoria:
!   oap
!
SUBROUTINE metodo_mi (self, R, P, FSomas)

  IMPLICIT NONE
  class(integracao_euler_imp), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, Rchute
  INTEGER :: i

  R1 = R
  P1 = P

  DO i = 1, 100
    ! Integrando as posicoes
    Rchute = R + self % h * self % m_inv * P1

    ! Integrando as velocidades
    P1 = P + self%h * (self % m2 * FSomas)

    ! Calcula as novas forcas
    FSomas = self%forcas(Rchute)

    IF (NORM2(R1 - Rchute) < self % h * self % h) EXIT
    R1 = Rchute
  END DO

  R = Rchute
  P = P1

END SUBROUTINE metodo_mi

END MODULE euler_imp