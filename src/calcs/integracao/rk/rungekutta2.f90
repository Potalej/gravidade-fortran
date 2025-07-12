! ************************************************************
!! METODO NUMERICO RK2
!
! Objetivos:
!   Aplicacao do metodo de Runge-Kutta de Segunda Ordem em Dois
!   Estagios (RK2 ou RK22)
!
! Modificado:
!   12 de julho de 2025
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
!   12 de julho de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE Iniciar (self, massas, G, h, potsoft, E0, J0, corrigir, corme, cormnt, colidir, colmodo, colmd, paralelo, mi)
  IMPLICIT NONE
  CLASS(integracao_rk2), INTENT(INOUT) :: self
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