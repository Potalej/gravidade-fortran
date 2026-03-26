! ************************************************************
!! METODO NUMERICO: Ruth 4a Ordem
!
! Objetivos:
!   Aplicacao do metodo simpletico de Ruth com 4a Ordem
!
! Modificado:
!   14 de setembro de 2024 (criado)
!   26 de marco de 2026 (atualizado)
!
! Autoria:
!   oap
!
MODULE ruth4
  USE tipos
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_ruth4

  ! Parametros
  REAL(pf) :: d1, d2, d3, d4
  REAL(pf) :: c1, c2, c3, c4

  TYPE, EXTENDS(integracao) :: integracao_ruth4

    CONTAINS
      PROCEDURE :: metodo, metodo_mi, atualizar_constantes

  END TYPE
  
CONTAINS

SUBROUTINE atualizar_constantes (self)
  IMPLICIT NONE
  class(integracao_ruth4), INTENT(IN) :: self
  REAL(pf) :: dois3, dif2, x

  dois3 = 2.0_pf ** (1.0_pf / 3.0_pf)
  dif2 = 2.0_pf - dois3
  x = (dois3 + 1.0_pf/dois3 - 1.0_pf)/6.0_pf
  
  c1 = 0
  c2 = 2.0_pf * x + 1.0_pf
  c3 = -4.0_pf * x - 1.0_pf
  c4 = 2.0_pf * x + 1.0_pf

  d1 = x + 0.5_pf
  d2 = -x
  d3 = -x
  d4 = x + 0.5_pf
END SUBROUTINE

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
  class(integracao_ruth4), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas

  ! i = 4
  P = P + d4 * self % h * FSomas
  R = R + c4 * self % h * P * self % massasInvertidas

  ! i = 3
  FSomas = self%forcas(R)
  P = P + d3 * self % h * FSomas
  R = R + c3 * self % h * P * self % massasInvertidas

  ! i = 2
  FSomas = self%forcas(R)
  P = P + d2 * self % h * FSomas
  R = R + c2 * self % h * P * self % massasInvertidas

  ! i = 1
  FSomas = self%forcas(R)
  P = P + d1 * self % h * FSomas
  ! R = R + c1 * self % h * P * self % massasInvertidas ! c1=0

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
  class(integracao_ruth4), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  
  ! i = 4
  P = P
  R = R + c4 * self % h * self % m_inv * P

  ! i = 3
  FSomas = self%forcas(R)
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

END SUBROUTINE metodo_mi

END MODULE ruth4