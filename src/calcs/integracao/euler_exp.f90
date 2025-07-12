! ************************************************************
!! METODO NUMERICO: Euler explicito
!
! Objetivos:
!   Aplicacao do metodo de Euler explicito.
!
! Modificado:
!   17 de novembro de 2024
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
!   14 de setembro de 2024
!
! Autoria:
!   oap
!
FUNCTION metodo (self, R, P, FSomas_ant)

  IMPLICIT NONE
  class(integracao_euler_exp), INTENT(IN) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, FSomas
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo
  INTEGER :: a

  ! Integrando as posicoes
  DO a = 1, self % N
    R1(a,:) = R(a,:) + self % h * P(a,:) / self % m(a)
  END DO

  ! Integrando as velocidades
  P1 = P + self%h * FSomas_ant

  metodo(1,:,:) = R1
  metodo(2,:,:) = P1
  metodo(3,:,:) = self % forcas(R1)

END FUNCTION metodo

! ************************************************************
!! Metodo numerico
!
! Objetivos:
!   Aplicacao do metodo em si.
!
! Modificado:
!   14 de setembro de 2024
!
! Autoria:
!   oap
!
FUNCTION metodo_mi (self, R, P, FSomas_ant)

  IMPLICIT NONE
  class(integracao_euler_exp), INTENT(IN) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, FSomas
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo_mi

  ! Integrando as posicoes
  R1 = R + self % h * self % m_inv * P

  ! Integrando as velocidades
  P1 = P + self%h * (self%m2 * FSomas)

  metodo_mi(1,:,:) = R1
  metodo_mi(2,:,:) = P1
  metodo_mi(3,:,:) = self%forcas(R)

END FUNCTION metodo_mi

END MODULE euler_exp