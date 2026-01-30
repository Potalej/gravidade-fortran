! ************************************************************
!! METODO NUMERICO MULTIPASSO: Adams-Bashforth 4
!
! Objetivos:
!   Aplicacao do metodo multipasso Adams-Bashforth 4
!
! Modificado:
!   29 de janeiro de 2026 (criado)
!   29 de janeiro de 2026 (modificado)
!
! Autoria:
!   oap
!
MODULE ab4
  USE tipos
  USE OMP_LIB
  USE integrador
  USE json_utils_mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_ab4

  TYPE, EXTENDS(integracao) :: integracao_ab4

    CONTAINS
      PROCEDURE :: metodo, metodo_mi
      PROCEDURE :: iniciar => iniciar_ab4

  END TYPE

CONTAINS

SUBROUTINE iniciar_ab4 (self, infos, timestep, massas)
  IMPLICIT NONE
  class(integracao_ab4), INTENT(INOUT) :: self
  TYPE(json_value), POINTER, INTENT(IN) :: infos
  REAL(pf), INTENT(IN), ALLOCATABLE :: massas(:)
  REAL(pf), INTENT(IN)              :: timestep
  
  CALL iniciar_base(self, infos, timestep, massas)
  self % multipasso = 4

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
  class(integracao_ab4), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  REAL(pf), DIMENSION(self%N, self%dim) :: aux

  ! Integra as posicoes
  aux =       55.0_pf/24.0_pf * P
  aux = aux - 59.0_pf/24.0_pf * self%Ps_ant(self % mp_rb(1),:,:)
  aux = aux + 37.0_pf/24.0_pf * self%Ps_ant(self % mp_rb(2),:,:)
  aux = aux -  9.0_pf/24.0_pf * self%Ps_ant(self % mp_rb(3),:,:)
  R = R + self % h * aux * self % massasInvertidas
  
  ! Atualiza o vetor de momentos
  self % Ps_ant(self % mp_rb(), :, :) = P

  ! Integra as velocidades
  aux =       55.0_pf/24.0_pf * FSomas
  aux = aux - 59.0_pf/24.0_pf * self%Fs_ant(self % mp_rb(1),:,:)
  aux = aux + 37.0_pf/24.0_pf * self%Fs_ant(self % mp_rb(2),:,:)
  aux = aux -  9.0_pf/24.0_pf * self%Fs_ant(self % mp_rb(3),:,:)
  P = P + self % h * aux

  ! Atualiza o vetor de forcas
  self % Fs_ant(self % mp_rb(), :, :) = FSomas

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
  class(integracao_ab4), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  REAL(pf), DIMENSION(self%N, self%dim) :: aux
  REAL(pf) :: cnst

  cnst = self % h * self % m_esc * self % m_esc

  ! Integra as posicoes
  aux =       55.0_pf/24.0_pf * P
  aux = aux - 59.0_pf/24.0_pf * self%Ps_ant(self % mp_rb(1),:,:)
  aux = aux + 37.0_pf/24.0_pf * self%Ps_ant(self % mp_rb(2),:,:)
  aux = aux -  9.0_pf/24.0_pf * self%Ps_ant(self % mp_rb(3),:,:)
  R = R + self % h * self % m_inv * aux

  ! Atualiza o vetor de momentos
  self % Ps_ant(self % mp_rb(), :, :) = P

  ! Integra as velocidades
  aux =       55.0_pf/24.0_pf * FSomas
  aux = aux - 59.0_pf/24.0_pf * self%Fs_ant(self % mp_rb(1),:,:)
  aux = aux + 37.0_pf/24.0_pf * self%Fs_ant(self % mp_rb(2),:,:)
  aux = aux -  9.0_pf/24.0_pf * self%Fs_ant(self % mp_rb(3),:,:)
  P = P + cnst * aux

  ! Atualiza o vetor de forcas
  self % Fs_ant(self % mp_rb(), :, :) = FSomas

  ! Calcula as novas forcas
  FSomas = self%forcas(R)

END SUBROUTINE metodo_mi

END MODULE ab4