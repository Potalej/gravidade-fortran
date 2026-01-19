! ************************************************************
!! METODO NUMERICO: Ruth 3a Ordem
!
! Objetivos:
!   Aplicacao do metodo simpletico de Ruth com 3a Ordem
!
! Modificado:
!   14 de setembro de 2024 (criado)
!   18 de janeiro de 2026 (atualizado)
!
! Autoria:
!   oap
!
MODULE ruth3
  USE tipos
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_ruth3

  ! Parametros
  REAL(pf) :: d1 = -0.041666666666666666_pf, c1 =  1.0_pf
  REAL(pf) :: d2 = 0.75_pf,                  c2 = -0.66666666666666666_pf
  REAL(pf) :: d3 = 0.29166666666666666_pf,   c3 =  0.66666666666666666_pf

  TYPE, EXTENDS(integracao) :: integracao_ruth3

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
  class(integracao_ruth3), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas

  ! i = 3
  P = P + d3 * self % h * FSomas
  R = R + c3 * self % h * P * self % massasInvertidas

  ! i = 2
  FSomas = self%forcas(R)
  P = P + d2 * self % h * FSomas
  R = R + c2 * self % h * P * self % massasInvertidas

  ! i = 1
  FSomas = self%forcas(R)
  P = P + d1 * self % h * FSomas
  R = R + c1 * self % h * P * self % massasInvertidas

  ! Calcula as novas forcas
  FSomas = self%forcas(R)

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
  class(integracao_ruth3), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  
  ! i = 3
  P = P + d3 * self % h * (self % m2 * FSomas)
  R = R + c3 * self % h * self % m_inv * P

  ! i = 2
  FSomas = self%forcas(R)
  P = P + d2 * self % h * (self % m2 * FSomas)
  R = R + c2 * self % h * self % m_inv * P

  ! i = 1
  FSomas = self%forcas(R)
  P = P + d1 * self % h * (self % m2 * FSomas)
  R = R + c1 * self % h * self % m_inv * P

  ! Calcula as novas forcas
  FSomas = self%forcas(R)

END SUBROUTINE metodo_mi

END MODULE ruth3