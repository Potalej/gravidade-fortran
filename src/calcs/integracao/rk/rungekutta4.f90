! ************************************************************
!! METODO NUMERICO RK4
!
! Objetivos:
!   Aplicacao do metodo de Runge-Kutta de Quarta Ordem em Quatro
!   Estagios (RK4 ou RK44)
!
! Modificado:
!   12 de julho de 2025
!
! Autoria:
!   oap
! 
MODULE rungekutta4
  USE tipos
  USE rungekutta
  USE funcoes_forca
  USE funcoes_forca_mi
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_rk4

  TYPE, EXTENDS(integracao) :: integracao_rk4
    
    ! Base do Runge-Kutta
    TYPE(RK) :: baseRK
    ! Modulos adicionais
    CONTAINS
      PROCEDURE :: Iniciar, metodo, metodo_mi

  END TYPE integracao_rk4

CONTAINS

! ************************************************************
!! Construtor da classe
!
! Objetivos:
!   Define o principal, salvando os valores e inicializando o
!   metodo.
!
! Modificado:
!   12 de julho de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE Iniciar (self, massas, G, h, potsoft, E0, J0, corrigir, corme, cormnt, colidir, colmodo, colmd, paralelo, mi)
  IMPLICIT NONE
  CLASS(integracao_rk4), INTENT(INOUT) :: self
  LOGICAL,INTENT(IN) :: corrigir, paralelo
  LOGICAL, INTENT(IN)          :: colidir
  CHARACTER(LEN=*), INTENT(IN) :: colmodo
  REAL(pf), allocatable :: massas(:)
  REAL(pf)              :: G, h, potsoft, colmd, E0, J0(3)
  REAL(pf) :: corme
  INTEGER :: cormnt
  LOGICAL :: mi
  INTEGER :: a

  ! quantidade de partículas
  self % N = SIZE(massas)
  
  ! massas
  self % mi = mi
  ALLOCATE(self % m (self % N))
  self % m = massas
  IF (mi) THEN
    self % m_esc = massas(1)
    self % m_inv = 1/self % m_esc
    self % m2 = self % m_esc * self % m_esc
  ENDIF

  ! gravidade
  self % G = G
  ! passo
  self % h = h
  ! Softening do potencial
  self % potsoft = potsoft

  ! Distancias
  ALLOCATE(self % distancias(INT(self%N * (self%N-1)/2)))

  ! Valores iniciais
  self % E0 = E0
  self % J0 = J0

  ! Se vai ou nao corrigir
  self % corrigir = corrigir
  self % corme = corme
  self % cormnt = cormnt

  ! Se vai ou nao colidir
  self % colidir = colidir
  self % colisoes_modo = colmodo
  self % colmd = colmd
  ALLOCATE(self % raios(self % N))
  DO a = 1, self % N
    self % raios(a) = self % colmd * massas(a)**(1.0_pf / 3.0_pf)
  END DO
  
  ! Inicia o base do RK 
  CALL self % baseRK % Iniciar(self % n, self % m, self % G, self % h, self % potsoft)
  
  ! Codigo paralelo
  self % paralelo = paralelo
  IF (paralelo) THEN
    IF (mi) THEN
      self % forcas_mi_funcao => forcas_mi_par
    ELSE
      self % forcas_funcao => forcas_par
    ENDIF
  ELSE
    IF (mi) THEN
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
  class(integracao_rk4), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self % N, self % dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self % N, self % dim) :: R1, P1
  REAL(pf), DIMENSION(3, self % N, self % dim) :: metodo

  ! componentes da integracao (kappas)
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_q, k2_q, k3_q, k4_q
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_p, k2_p, k3_p, k4_p

  k1_q = P * self % baseRK % massasInvertidas
  k1_p = FSomas_ant

  k2_q = (P + 0.5_pf * self % h * k1_p) * self % baseRK % massasInvertidas
  k2_p = self % forcas (R + 0.5_pf * self % h * k1_q)

  k3_q = (P + 0.5_pf * self % h * k2_p) * self % baseRK % massasInvertidas
  k3_p = self % forcas (R + 0.5_pf * self % h * k2_q)

  k4_q = (P + self % h * k3_p) * self % baseRK % massasInvertidas
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
  IMPLICIT NONE
  class(integracao_rk4), INTENT(INOUT) :: self
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

END module rungekutta4