! ************************************************************
!! METODO NUMERICO RK4
!
! Objetivos:
!   Aplicacao do metodo de Runge-Kutta de Quarta Ordem em Quatro
!   Estagios (RK4 ou RK44)
!
! Modificado:
!   08 de agosto de 2025
!
! Autoria:
!   oap
! 
MODULE rungekutta4
  USE tipos
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_rungekutta4

  TYPE, EXTENDS(integracao) :: integracao_rungekutta4

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
  class(integracao_rungekutta4), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self % N, self % dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self % N, self % dim) :: R1, P1
  REAL(pf), DIMENSION(3, self % N, self % dim) :: metodo

  ! componentes da integracao (kappas)
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_q, k2_q, k3_q, k4_q
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_p, k2_p, k3_p, k4_p

  k1_q = P * self % massasInvertidas
  k1_p = FSomas_ant

  k2_q = (P + 0.5_pf * self % h * k1_p) * self % massasInvertidas
  k2_p = self % forcas (R + 0.5_pf * self % h * k1_q)

  k3_q = (P + 0.5_pf * self % h * k2_p) * self % massasInvertidas
  k3_p = self % forcas (R + 0.5_pf * self % h * k2_q)

  k4_q = (P + self % h * k3_p) * self % massasInvertidas
  k4_p = self % forcas (R + self % h * k3_q)

  ! fator para integracao
  R1 = R + (self % h / 6.0_pf) * (k1_q + 2.0_pf * k2_q + 2.0_pf * k3_q + k4_q)
  P1 = P + (self % h / 6.0_pf) * (k1_p + 2.0_pf * k2_p + 2.0_pf * k3_p + k4_p)

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
  class(integracao_rungekutta4), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self % N, self % dim) :: R1, P1
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo_mi
  
  ! componentes da integracao (kappas)
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_q, k2_q, k3_q, k4_q
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_p, k2_p, k3_p, k4_p

  k1_q = P * self % m_inv
  k1_p = FSomas_ant * self % m2

  k2_q = (P + 0.5_pf * self % h * k1_p) * self % m_inv
  k2_p = self % forcas (R + 0.5_pf * self % h * k1_q)  * self % m2

  k3_q = (P + 0.5_pf * self % h * k2_p) * self % m_inv
  k3_p = self % forcas (R + 0.5_pf * self % h * k2_q)  * self % m2

  k4_q = (P + self % h * k3_p) * self % m_inv
  k4_p = self % forcas (R + self % h * k3_q) * self % m2

  ! fator para integracao
  R1 = R + (self % h / 6.0_pf) * (k1_q + 2.0_pf * k2_q + 2.0_pf * k3_q + k4_q)
  P1 = P + (self % h / 6.0_pf) * (k1_p + 2.0_pf * k2_p + 2.0_pf * k3_p + k4_p)

  metodo_mi(1,:,:) = R1
  metodo_mi(2,:,:) = P1
  metodo_mi(3,:,:) = self % forcas(R1)
END FUNCTION metodo_mi

END MODULE rungekutta4