! ************************************************************
!! METODO NUMERICO: Velocity-Verlet (DKD)
!
! Objetivos:
!   Aplicacao do metodo simpletico de Velocity-Verlet usando o
!   Drift-Kick-Drift (DKD). Eh o padrao usado pelo Rebound.
!
! Modificado:
!   15 de setembro de 2026 (criado)
!   18 de setembro de 2026 (modificado)
!
! Autoria:
!   oap
!
MODULE verlet_dkd
  USE tipos
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_verlet_dkd

  TYPE, EXTENDS(integracao) :: integracao_verlet_dkd

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
  class(integracao_verlet_dkd), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas

  R = R + 0.5_pf * self % h * P * self%massasInvertidas
  FSomas = self%forcas(R)
  P = P + self % h * FSomas
  R = R + 0.5_pf * self % h * P * self%massasInvertidas

END SUBROUTINE metodo

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
SUBROUTINE metodo_mi (self, R, P, FSomas)

  IMPLICIT NONE
  class(integracao_verlet_dkd), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  
  ! integrando as posicoes
  R = R + 0.5_pf * self%h * self%m_inv * P
  
  ! calcula as novas forcas
  FSomas = self%forcas(R)
  
  ! velocidades
  P = P + self%h * self%m_esc * self%m_esc * FSomas
  
  ! novas posicoes
  R = R + 0.5_pf * self%h * self%m_inv * P

END SUBROUTINE metodo_mi

END MODULE verlet_dkd