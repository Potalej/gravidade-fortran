! ************************************************************
!! METODO NUMERICO: Ruth 3a Ordem
!
! Objetivos:
!   Aplicacao do metodo simpletico de Ruth com 3a Ordem
!
! Modificado:
!   14 de setembro de 2024
!
! Autoria:
!   oap
!
MODULE ruth3
  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_ruth3

  TYPE, EXTENDS(integracao) :: integracao_ruth3

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
  class(integracao_ruth3), INTENT(IN) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, FSomas_prox
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo
  INTEGER :: a
  REAL(pf) :: d3 = 0.29166666666666666_pf, c3 = 0.66666666666666666_pf
  REAL(pf) :: d2 = 0.75_pf, c2 = - 0.66666666666666666_pf
  REAL(pf) :: d1 = -0.041666666666666666_pf, c1 = 1.0_pf

  ! i = 3

  ! Momentos
  P1 = P + d3 * self%h * FSomas_ant
  
  ! R1 = R + c3 * self % h * P1 * self % massasInvertidas

  ! Posicoes
  DO a = 1, self % N
    R1(a,:) = R(a,:) + c3 * self % h * P1(a,:) / self % m(a)
  END DO


  ! i = 2
  FSomas_prox = self%forcas(R1)

  ! Momentos
  P1 = P1 + d2 * FSomas_prox * self % h

  ! R1 = R1 + c2 * self % h * P1 * self % massasInvertidas

  ! Posicoes
  DO a = 1, self % N
    R1(a,:) = R1(a,:) + c2 * self % h * P1(a,:) / self % m(a)
  END DO


  ! i = 1
  FSomas_prox = self%forcas(R1)

  ! Momentos
  P1 = P1 + d1 * FSomas_prox * self % h

  ! R1 = R1 + c1 * self % h * P1 * self % massasInvertidas

  ! Posicoes
  DO a = 1, self % N
    R1(a,:) = R1(a,:) + c1 * self % h * P1(a,:) / self % m(a)
  END DO


  ! Calcula as novas forcas
  FSomas_prox = self%forcas(R1)


  metodo(1,:,:) = R1
  metodo(2,:,:) = P1
  metodo(3,:,:) = FSomas_prox

END FUNCTION metodo

END MODULE ruth3