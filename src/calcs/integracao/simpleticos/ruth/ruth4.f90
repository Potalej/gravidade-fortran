! ************************************************************
!! METODO NUMERICO: Ruth 4a Ordem
!
! Objetivos:
!   Aplicacao do metodo simpletico de Ruth com 4a Ordem
!
! Modificado:
!   14 de setembro de 2024 (criado)
!   18 de janeiro de 2026 (atualizado)
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
  REAL(pf) :: d1 =  1.3512071919596578_pf, c1 =  0.67560359597982889_pf
  REAL(pf) :: d2 = -1.7024143839193153_pf, c2 = -0.17560359597982883_pf
  REAL(pf) :: d3 =  1.3512071919596578_pf, c3 = -0.17560359597982883_pf
  REAL(pf) :: d4 =  0.0_pf,                c4 =  0.67560359597982889_pf

  TYPE, EXTENDS(integracao) :: integracao_ruth4

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
  class(integracao_ruth4), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas

  ! i = 4
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
  R = R + c1 * self % h * P * self % massasInvertidas

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