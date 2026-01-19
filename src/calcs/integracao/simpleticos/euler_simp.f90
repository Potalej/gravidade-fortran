! ************************************************************
!! METODO NUMERICO: Euler Simpletico
!
! Objetivos:
!   Aplicacao do metodo simpletico de Euler.
!
! Modificado:
!   14 de setembro de 2024 (criado)
!   18 de janeiro de 2026 (modificado)
!
! Autoria:
!   oap
!
MODULE euler_simp
  USE tipos
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_euler_simp

  TYPE, EXTENDS(integracao) :: integracao_euler_simp

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
  class(integracao_euler_simp), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  INTEGER :: a

  ! Integrando as velocidades
  P = P + self%h * FSomas

  ! Integrando as posicoes
  DO a = 1, self % N
    R(a,:) = R(a,:) + self % h * P(a,:) / self % m(a)
  END DO

  FSomas = self % forcas(R)
END SUBROUTINE metodo

! ************************************************************
!! Metodo numerico (massas iguais)
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
  class(integracao_euler_simp), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  
  ! Integrando as velocidades
  P = P + self%h * (self%m2 * FSomas)

  ! Integrando as posicoes
  R = R + self % h * P * self % m_inv

  FSomas = self % forcas(R)
END SUBROUTINE metodo_mi

END MODULE euler_simp