! ************************************************************
!! METODO NUMERICO: Euler Simpletico
!
! Objetivos:
!   Aplicacao do metodo simpletico de Euler.
!
! Modificado:
!   14 de setembro de 2024 (criado)
!   12 de julho de 2025 (modificado)
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
!   12 de julho de 2025
!
! Autoria:
!   oap
!
FUNCTION metodo (self, R, P, FSomas_ant)

  IMPLICIT NONE
  class(integracao_euler_simp), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo
  INTEGER :: a

  ! Integrando as velocidades
  P1 = P + self%h * FSomas_ant

  ! Integrando as posicoes
  DO a = 1, self % N
    R1(a,:) = R(a,:) + self % h * P1(a,:) / self % m(a)
  END DO

  metodo(1,:,:) = R1
  metodo(2,:,:) = P1
  metodo(3,:,:) = self % forcas(R1)

END FUNCTION metodo

! ************************************************************
!! Metodo numerico (massas iguais)
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
  class(integracao_euler_simp), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo_mi
  
  ! Integrando as velocidades
  P1 = P + self%h * (self%m2 * FSomas_ant)

  ! Integrando as posicoes
  R1 = R + self % h * P1 * self % m_inv

  metodo_mi(1,:,:) = R1
  metodo_mi(2,:,:) = P1
  metodo_mi(3,:,:) = self % forcas(R1)

END FUNCTION metodo_mi

END MODULE euler_simp