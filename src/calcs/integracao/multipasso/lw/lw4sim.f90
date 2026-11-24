! ************************************************************
!! METODO NUMERICO: Lambert-Watson Simetrico de 4ª Ordem
!
! Objetivos:
!   Aplicacao do metodo simetrico de Lambert e Watson de 4ª 
!   Ordem. Publicado em Lambert e Watson (1976), metodo IV com
!   a = 0. Acesso em Hairer, Lubich, Wanner, p. 570 (C).
!
! Modificado:
!   23 de novembro de 2026 (criado)
!   23 de novembro de 2026 (modificado)
!
! Autoria:
!   oap
!
MODULE lw4sim
  USE tipos
  USE OMP_LIB
  USE integrador
  USE json_utils_mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_lw4sim

  TYPE, EXTENDS(integracao) :: integracao_lw4sim

    CONTAINS
      PROCEDURE :: metodo, metodo_mi
      PROCEDURE :: iniciar => iniciar_lw4sim

  END TYPE

CONTAINS

SUBROUTINE iniciar_lw4sim (self, infos, timestep, massas, auxiliar)
  IMPLICIT NONE
  class(integracao_lw4sim), INTENT(INOUT) :: self
  TYPE(json_value), POINTER, INTENT(IN) :: infos
  REAL(pf), INTENT(IN), ALLOCATABLE :: massas(:)
  LOGICAL,  INTENT(IN), OPTIONAL :: auxiliar
  REAL(pf), INTENT(IN)              :: timestep
  
  CALL iniciar_base(self, infos, timestep, massas, .FALSE.)
  self % multipasso = 4

  ALLOCATE(self % rs_ant(self%multipasso, self % N, self % dim))
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
  class(integracao_lw4sim), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  REAL(pf), DIMENSION(self%N, self%dim) :: aux, ace
  REAL(pf) :: old(self % N, self % dim), v(self % N, self % dim)
  INTEGER :: a
  
  aux =       2.0_pf * R
  aux = aux - 2.0_pf * self % Rs_ant(self % mp_rb(1),:,:)
  aux = aux + 2.0_pf * self % Rs_ant(self % mp_rb(2),:,:)
  aux = aux -          self % Rs_ant(self % mp_rb(3),:,:)

  ace =       (7.0_pf / 6.0_pf) * FSomas
  ace = ace - (1.0_pf / 3.0_pf) * self % Fs_ant(self % mp_rb(1),:,:)
  ace = ace + (7.0_pf / 6.0_pf) * self % Fs_ant(self % mp_rb(2),:,:)
  ace = ace * self % massasInvertidas

  old = self % Rs_ant(self % mp_rb(3),:,:)
  self % Rs_ant(self % mp_rb(), :, :) = R
  self % Fs_ant(self % mp_rb(), :, :) = FSomas

  R = aux + self % h * self % h * ace
  FSomas = self % forcas(R)

  self % Ps_ant(self % mp_rb(), :, :) = P
  v =     25.0_pf * R
  v = v - 48.0_pf * self % Rs_ant(self % mp_rb(3),  :,:)
  v = v + 36.0_pf * self % Rs_ant(self % mp_rb(1), :,:)
  v = v - 16.0_pf * self % Rs_ant(self % mp_rb(2), :,:)
  v = v +  3.0_pf * old
  v = v / (12.0_pf * self % h)
  DO a = 1, self % N
    P(a,:) = v(a,:) * self % m(a)
  END DO


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
  class(integracao_lw4sim), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  REAL(pf), DIMENSION(self%N, self%dim) :: aux, ace
  REAL(pf) :: old(self % N, self % dim), v(self % N, self % dim)
  REAL(pf) :: cnst

  cnst = self % h * self % h * self % m_esc
  
  aux =       2.0_pf * R
  aux = aux - 2.0_pf * self % Rs_ant(self % mp_rb(1),:,:)
  aux = aux + 2.0_pf * self % Rs_ant(self % mp_rb(2),:,:)
  aux = aux -          self % Rs_ant(self % mp_rb(3),:,:)

  ace =       (7.0_pf / 6.0_pf) * FSomas
  ace = ace - (1.0_pf / 3.0_pf) * self % Fs_ant(self % mp_rb(1),:,:)
  ace = ace + (7.0_pf / 6.0_pf) * self % Fs_ant(self % mp_rb(2),:,:)

  old = self % Rs_ant(self % mp_rb(3),:,:)
  self % Rs_ant(self % mp_rb(), :, :) = R
  self % Fs_ant(self % mp_rb(), :, :) = FSomas

  R = aux + cnst * ace
  FSomas = self % forcas(R)

  v =     25.0_pf * R
  v = v - 48.0_pf * self % Rs_ant(self % mp_rb(3),  :,:)
  v = v + 36.0_pf * self % Rs_ant(self % mp_rb(1), :,:)
  v = v - 16.0_pf * self % Rs_ant(self % mp_rb(2), :,:)
  v = v +  3.0_pf * old
  v = v / (12.0_pf * self % h)
  P = v * self % m_esc

END SUBROUTINE metodo_mi

END MODULE lw4sim