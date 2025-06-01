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

  ! i = 4
  P1 = P
  R1 = R + c4 * self % h * P1 * self % massasInvertidas

  ! i = 3
  FSomas_prox = self%forcas(R1)
  P1 = P1 + d3 * self % h * FSomas_prox
  R1 = R1 + c3 * self % h * P1 * self % massasInvertidas

  ! i = 2
  FSomas_prox = self%forcas(R1)
  P1 = P1 + d2 * self % h * FSomas_prox
  R1 = R1 + c2 * self % h * P1 * self % massasInvertidas

  ! i = 1
  FSomas_prox = self%forcas(R1)
  P1 = P1 + d1 * self % h * FSomas_prox
  R1 = R1 + c1 * self % h * P1 * self % massasInvertidas

  ! Calcula as novas forcas
  FSomas_prox = self%forcas(R1)

  metodo(1,:,:) = R1
  metodo(2,:,:) = P1
  metodo(3,:,:) = FSomas_prox

END FUNCTION metodo

! ************************************************************
!! Metodo numerico (massas iguais)
!
! Objetivos:
!   Aplicacao do metodo em si.
!
! Modificado:
!   01 de junho de 2025
!
! Autoria:
!   oap
!
FUNCTION metodo_mi (self, R, P, FSomas_ant)
  IMPLICIT NONE
  class(integracao_ruth4), INTENT(IN) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, FSomas_prox
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo_mi
  
  ! i = 4
  P1 = P
  R1 = R + c4 * self % h * self % m_inv * P1

  ! i = 3
  FSomas_prox = self%forcas(R1)
  P1 = P1 + d3 * self % h * (self % m2 * FSomas_prox)
  R1 = R1 + c3 * self % h * self % m_inv * P1

  ! i = 2
  FSomas_prox = self%forcas(R1)
  P1 = P1 + d2 * self % h * (self % m2 * FSomas_prox)
  R1 = R1 + c2 * self % h * self % m_inv * P1

  ! i = 1
  FSomas_prox = self%forcas(R1)
  P1 = P1 + d1 * self % h * (self % m2 * FSomas_prox)
  R1 = R1 + c1 * self % h * self % m_inv * P1

  ! Calcula as novas forcas
  FSomas_prox = self%forcas(R1)

  metodo_mi(1,:,:) = R1
  metodo_mi(2,:,:) = P1
  metodo_mi(3,:,:) = FSomas_prox

END FUNCTION metodo_mi

END MODULE ruth4