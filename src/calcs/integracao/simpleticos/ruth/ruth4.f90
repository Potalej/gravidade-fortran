! ************************************************************
!! METODO NUMERICO: Ruth 4a Ordem
!
! Objetivos:
!   Aplicacao do metodo simpletico de Ruth com 4a Ordem
!
! Modificado:
!   14 de setembro de 2024
!
! Autoria:
!   oap
!
MODULE ruth4
  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_ruth4

  TYPE, EXTENDS(integracao) :: integracao_ruth4

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
  class(integracao_ruth4), INTENT(IN) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, FSomas_prox
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo
  INTEGER :: a
  REAL(pf) :: d4 = 0.0_pf, c4 = 0.67560359597982889_pf
  REAL(pf) :: d3 = 1.3512071919596578_pf, c3 = -0.17560359597982883_pf
  REAL(pf) :: d2 = -1.7024143839193153_pf, c2 = -0.17560359597982883_pf
  REAL(pf) :: d1 = 1.3512071919596578_pf, c1 = 0.67560359597982889_pf

  ! i = 4

  ! Momentos: p_ = p + d4 F(q) h
  P1 = P
  
  ! Posicoes: q_ = q + c4 p_/m h
  DO a = 1, self % N
    R1(a,:) = R(a,:) + c4 * self % h * P1(a,:) / self % m(a)
  END DO

  ! i = 3
  FSomas_prox = self%forcas(R1)

  ! Momentos: p_ = p + d3 F(q) h
  P1 = P1 + d3*self%h*FSomas_prox
  
  ! Posicoes: q_ = q + c3 p_/m h
  DO a = 1, self % N
    R1(a,:) = R1(a,:) + c3 * self % h * P1(a,:) / self % m(a)
  END DO


  ! i = 2
  FSomas_prox = self%forcas(R1)

  ! Momentos: p_ = p_ + d2 F(q_) h
  P1 = P1 + d2 * FSomas_prox * self % h

  ! Posicoes: q_ = q_ + c2 p_/m h
  DO a = 1, self % N
    R1(a,:) = R1(a,:) + c2 * self % h * P1(a,:) / self % m(a)
  END DO


  ! i = 1
  FSomas_prox = self%forcas(R1)

  ! Momentos: p_ = p_ + d1 F(q_) h
  P1 = P1 + d1 * FSomas_prox * self % h

  ! Posicoes: q_ = q_ + c1 p_/m h
  DO a = 1, self % N
    R1(a,:) = R1(a,:) + c1 * self % h * P1(a,:) / self % m(a)
  END DO


  ! Calcula as novas forcas
  FSomas_prox = self%forcas(R1)


  metodo(1,:,:) = R1
  metodo(2,:,:) = P1
  metodo(3,:,:) = FSomas_prox

END FUNCTION metodo

END MODULE ruth4