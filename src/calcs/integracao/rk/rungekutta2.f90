! ************************************************************
!! METODO NUMERICO RK2
!
! Objetivos:
!   Aplicacao do metodo de Runge-Kutta de Segunda Ordem em Dois
!   Estagios (RK2 ou RK22)
!
! Modificado:
!   18 de janeiro de 2026
!
! Autoria:
!   oap
! 
MODULE rungekutta2
  USE tipos
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_rungekutta2

  TYPE, EXTENDS(integracao) :: integracao_rungekutta2

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
  class(integracao_rungekutta2), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self % N, self % dim), INTENT(INOUT) :: R, P, FSomas

  ! componentes da integracao (kappas)
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_q, k2_q
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_p, k2_p

  k1_q = P * self % massasInvertidas
  k1_p = FSomas

  k2_q = (P + 0.5_pf * self % h * k1_p) * self % massasInvertidas
  k2_p = self % forcas(R + 0.5_pf * self % h * k1_q)

  R = R + self % h * k2_q
  P = P + self % h * k2_p

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
  IMPLICIT NONE
  class(integracao_rungekutta2), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  
  ! componentes da integracao (kappas)
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_q, k2_q
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_p, k2_p
  
  k1_q = P * self % m_inv
  k1_p = FSomas * self % m2

  k2_q = (p + 0.5_pf * self % h * k1_p) * self % m_inv
  k2_p = self % forcas(R + 0.5_pf * self % h * k1_q) * self % m2

  R = R + self % h * k2_q
  P = P + self % h * k2_p

  FSomas = self % forcas(R)

END SUBROUTINE metodo_mi

END MODULE rungekutta2