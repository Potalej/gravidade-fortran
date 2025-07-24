! ************************************************************
!! METODO NUMERICO RK2
!
! Objetivos:
!   Aplicacao do metodo de Runge-Kutta de Segunda Ordem em Dois
!   Estagios (RK2 ou RK22)
!
! Modificado:
!   24 de julho de 2025
!
! Autoria:
!   oap
! 
MODULE rungekutta2
  USE tipos
  USE rungekutta
  USE funcoes_forca
  USE funcoes_forca_mi
  USE integrador
  USE json_utils_mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_rk2

  TYPE, EXTENDS(integracao) :: integracao_rk2
    
    ! Base do Runge-Kutta
    TYPE(RK) :: baseRK
    ! Modulos adicionais
    CONTAINS
      PROCEDURE :: Iniciar, metodo, metodo_mi

  END TYPE integracao_rk2

CONTAINS

! ************************************************************
!! Construtor da classe
!
! Objetivos:
!   Define o principal, salvando os valores e inicializando o
!   metodo.
!
! Modificado:
!   24 de julho de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE Iniciar (self, infos, timestep, massas, E0, J0)
  IMPLICIT NONE
  CLASS(integracao_rk2), INTENT(INOUT) :: self
  TYPE(json_value), POINTER :: infos
  REAL(pf) :: timestep
  REAL(pf), allocatable :: massas(:)
  REAL(pf)              :: E0, J0(3)
  INTEGER :: a, i

  LOGICAL :: encontrado
  REAL(pf) :: colmd, densidade
  REAL(pf) :: PI = 4.D0*DATAN(1.D0)

  ! quantidade de partículas
  self % N = SIZE(massas)
  
  ! massas
  CALL json % get(infos, 'massas_iguais', self % mi, encontrado)
  IF (.NOT. encontrado) self % mi = .FALSE.
  ALLOCATE(self % m(self % N))
  self % m = massas
  IF (self % mi) THEN
    self % m_esc = massas(1)
    self % m_inv = 1/self % m_esc
    self % m2 = self % m_esc * self % m_esc
  ENDIF

  ! gravidade
  self % G = json_get_float(infos, 'G')
  ! Passo
  self % h = timestep
  ! Softening do potencial
  self % potsoft = json_get_float(infos, 'integracao.amortecedor')
  self % potsoft2 = self%potsoft * self%potsoft

  ! Valores iniciais
  self % E0 = E0
  self % J0 = J0

  ! Distancias
  ALLOCATE(self % distancias(INT(self%N * (self%N-1)/2)))

  ! Se vai ou nao corrigir
  CALL json % get(infos, 'correcao.corrigir', self % corrigir)
  self % corme = json_get_float(infos, 'correcao.margem_erro')
  CALL json % get(infos, 'correcao.max_num_tentativas', self % cormnt)

  ! Colisoes
  CALL json % get(infos, 'colisoes.colidir', self % colidir)
  self % colisoes_modo = json_get_string(infos, 'colisoes.metodo')
  densidade = json_get_float(infos, 'colisoes.densidade')
  colmd = (0.75_pf / (PI * densidade))**(1.0_pf/3.0_pf)
  self % colmd = colmd

  ! Deixando os raios calculados de antemao
  ALLOCATE(self % raios(self % N))
  DO a = 1, self % N
    self % raios(a) = self % colmd * massas(a)**(1.0_pf / 3.0_pf)
  END DO
  
  ! Inicia o base do RK 
  CALL self % baseRK % Iniciar(self % n, self % m, self % G, self % h, self % potsoft)
  
  ! Codigo paralelo
  CALL json % get(infos, 'paralelo', self % paralelo)
  IF (self % paralelo) THEN
    IF (self % mi) THEN
      self % forcas_mi_funcao => forcas_mi_par
    ELSE
      self % forcas_funcao => forcas_par
    ENDIF
  ELSE
    IF (self % mi) THEN
      self % forcas_mi_funcao => forcas_mi_seq
    ELSE
      self % forcas_funcao => forcas_seq
    ENDIF
  ENDIF

END SUBROUTINE Iniciar


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

  IMPLICIT NONE
  class(integracao_rk2), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self % N, self % dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self % N, self % dim) :: R1, P1
  REAL(pf), DIMENSION(3, self % N, self % dim) :: metodo

  ! componentes da integracao (kappas)
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_q
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_p

  k1_p = P * self % baseRK % massasInvertidas + 0.5_pf * self % h * FSomas_ant
  R1 = R + self % h * k1_p 
  
  k1_q = R + 0.5_pf * self % h * P * self % baseRK % massasInvertidas
  P1 = P + self % h * self % forcas(k1_q)

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
  IMPLICIT NONE
  class(integracao_rk2), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self % N, self % dim) :: R1, P1
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo_mi
  
  ! componentes da integracao (kappas)
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_q
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_p

  k1_p = P * self % m_inv + 0.5_pf * self % h * FSomas_ant * self % m2
  R1 = R + self % h * k1_p
  
  k1_q = R + 0.5_pf * self % h * P * self % m_inv
  P1 = P + self % h * self % forcas(k1_q) * self % m2

  metodo_mi(1,:,:) = R1
  metodo_mi(2,:,:) = P1
  metodo_mi(3,:,:) = self % forcas(R1)

END FUNCTION metodo_mi

END module rungekutta2