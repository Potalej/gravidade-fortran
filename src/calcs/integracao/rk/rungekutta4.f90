! ************************************************************
!! METODO NUMERICO RK4
!
! Objetivos:
!   Aplicacao do metodo de Runge-Kutta de Quarta Ordem em Quatro
!   Estagios (RK4 ou RK44)
!
! Modificado:
!   18 de janeiro de 2026
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
!   18 de janeiro de 2026
!
! Autoria:
!   oap
!
SUBROUTINE metodo (self, R, P, FSomas)
  class(integracao_rungekutta4), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self % N, self % dim), INTENT(INOUT) :: R, P, FSomas

  ! componentes da integracao (kappas)
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_q, k2_q, k3_q, k4_q
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_p, k2_p, k3_p, k4_p

  k1_q = P * self % massasInvertidas
  k1_p = FSomas

  k2_q = (P + 0.5_pf * self % h * k1_p) * self % massasInvertidas
  k2_p = self % forcas (R + 0.5_pf * self % h * k1_q)

  k3_q = (P + 0.5_pf * self % h * k2_p) * self % massasInvertidas
  k3_p = self % forcas (R + 0.5_pf * self % h * k2_q)

  k4_q = (P + self % h * k3_p) * self % massasInvertidas
  k4_p = self % forcas (R + self % h * k3_q)

  ! fator para integracao
  R = R + (self % h / 6.0_pf) * (k1_q + 2.0_pf * k2_q + 2.0_pf * k3_q + k4_q)
  P = P + (self % h / 6.0_pf) * (k1_p + 2.0_pf * k2_p + 2.0_pf * k3_p + k4_p)

  FSomas = self % forcas(R)
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
  class(integracao_rungekutta4), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  
  ! componentes da integracao (kappas)
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_q, k2_q, k3_q, k4_q
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_p, k2_p, k3_p, k4_p

  k1_q = P * self % m_inv
  k1_p = FSomas * self % m2

  k2_q = (P + 0.5_pf * self % h * k1_p) * self % m_inv
  k2_p = self % forcas (R + 0.5_pf * self % h * k1_q)  * self % m2

  k3_q = (P + 0.5_pf * self % h * k2_p) * self % m_inv
  k3_p = self % forcas (R + 0.5_pf * self % h * k2_q)  * self % m2

  k4_q = (P + self % h * k3_p) * self % m_inv
  k4_p = self % forcas (R + self % h * k3_q) * self % m2

  ! fator para integracao
  R = R + (self % h / 6.0_pf) * (k1_q + 2.0_pf * k2_q + 2.0_pf * k3_q + k4_q)
  P = P + (self % h / 6.0_pf) * (k1_p + 2.0_pf * k2_p + 2.0_pf * k3_p + k4_p)

  FSomas = self % forcas(R)
END SUBROUTINE metodo_mi

END MODULE rungekutta4