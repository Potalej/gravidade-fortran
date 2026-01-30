! ************************************************************
!! METODO NUMERICO MULTIPASSO: Adams-Bashforth 2
!
! Objetivos:
!   Aplicacao do metodo multipasso Adams-Bashforth 2
!
! Modificado:
!   29 de janeiro de 2026 (criado)
!   29 de janeiro de 2026 (modificado)
!
! Autoria:
!   oap
!
MODULE ab2
  USE tipos
  USE OMP_LIB
  USE integrador
  USE json_utils_mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_ab2

  TYPE, EXTENDS(integracao) :: integracao_ab2

    CONTAINS
      PROCEDURE :: metodo, metodo_mi
      PROCEDURE :: iniciar => iniciar_ab2

  END TYPE

CONTAINS

SUBROUTINE iniciar_ab2 (self, infos, timestep, massas)
  IMPLICIT NONE
  class(integracao_ab2), INTENT(INOUT) :: self
  TYPE(json_value), POINTER, INTENT(IN) :: infos
  REAL(pf), INTENT(IN), allocatable :: massas(:)
  REAL(pf), INTENT(IN)              :: timestep
  
  CALL iniciar_base(self, infos, timestep, massas)
  self % multipasso = 2

  ALLOCATE(self % ps_ant(self%multipasso, self % N, self % dim))
  ALLOCATE(self % fs_ant(self%multipasso, self % N, self % dim))

END SUBROUTINE

! ************************************************************
!! Metodo numerico
!
! Objetivos:
!   Aplicacao do metodo em si.
!
! Modificado:
!   29 de janeiro de 2026
!
! Autoria:
!   oap
!
SUBROUTINE metodo (self, R, P, FSomas)

  IMPLICIT NONE
  class(integracao_ab2), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas

  ! Integra as posicoes
  R = R + self % h * (1.5_pf * P - 0.5_pf * self % Ps_ant(1,:,:)) * self % massasInvertidas
  
  ! Atualiza o vetor de momentos
  self % Ps_ant(1,:,:) = P

  ! Integra as velocidades
  P = P + self % h * (1.5_pf * FSomas - 0.5_pf * self % Fs_ant(1,:,:))

  ! Atualiza o vetor de forcas
  self % Fs_ant(1,:,:) = FSomas

  ! Calcula as novas forcas
  FSomas = self%forcas(R)

END SUBROUTINE metodo

! ************************************************************
!! Metodo numerico
!
! Objetivos:
!   Aplicacao do metodo em si.
!
! Modificado:
!   29 de janeiro de 2026
!
! Autoria:
!   oap
!
SUBROUTINE metodo_mi (self, R, P, FSomas)

  IMPLICIT NONE
  class(integracao_ab2), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  REAL(pf) :: cnst

  ! Integra as posicoes
  R = R + self % h * self % m_inv * (1.5_pf * P - 0.5_pf * self % Ps_ant(1,:,:))
  
  ! Atualiza o vetor de momentos
  self % Ps_ant(1,:,:) = P

  ! Integra as velocidades
  cnst = self % h * self % m_esc * self % m_esc
  P = P + cnst * (1.5_pf * FSomas - 0.5_pf * self % Fs_ant(1,:,:))

  ! Atualiza o vetor de forcas
  self % Fs_ant(1,:,:) = FSomas

  ! Calcula as novas forcas
  FSomas = self%forcas(R)

END SUBROUTINE metodo_mi

END MODULE ab2