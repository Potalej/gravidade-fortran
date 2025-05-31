! ************************************************************
!! METODO NUMERICO: Euler Implicito
!
! Objetivos:
!   Aplicacao do metodo de Euler Implicito.
!
! Modificado:
!   17 de novembro de 2024
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
      PROCEDURE :: metodo

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
  class(integracao_euler_imp), INTENT(IN) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, FSomas
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo
  INTEGER :: a, i

  R1 = R
  P1 = P

  DO i = 1, 100
    ! Calcula as forcas
    FSomas = self%forcas(R1)

    ! Integrando as posicoes
    DO a = 1, self % N
      R1(a,:) = R(a,:) + self % h * P1(a,:) / self % m(a)
    END DO

    ! Integrando as velocidades
    P1 = P + self%h * FSomas

  END DO

  metodo(1,:,:) = R1
  metodo(2,:,:) = P1

END FUNCTION metodo

END MODULE euler_imp