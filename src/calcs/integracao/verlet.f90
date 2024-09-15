! ************************************************************
!! METODO NUMERICO: Velocity-Verlet
!
! Objetivos:
!   Aplicacao do metodo simpletico de Velocity-Verlet.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
MODULE verlet
  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_verlet

  TYPE, EXTENDS(integracao) :: integracao_verlet

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
!   15 de marco de 2024
!
! Autoria:
!   oap
!
FUNCTION metodo (self, R, P, FSomas_ant)

  IMPLICIT NONE
  class(integracao_verlet), INTENT(IN) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, FSomas_prox
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo
  INTEGER :: a

  ! Integrando as posicoes
  DO a = 1, self % N
    R1(a,:) = R(a,:) + self % h * P(a,:) / self % m(a) + 0.5*(self%h**2)*FSomas_ant(a,:)/self%m(a)
  END DO

  ! Calcula as novas forcas
  FSomas_prox = self%forcas(R1)

  ! Integrando as velocidades
  P1 = P + 0.5*self%h*(FSomas_ant + FSomas_prox)

  metodo(1,:,:) = R1
  metodo(2,:,:) = P1
  metodo(3,:,:) = FSomas_prox

END FUNCTION metodo

END MODULE verlet