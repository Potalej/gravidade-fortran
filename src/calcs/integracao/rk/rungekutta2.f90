! ************************************************************
!! METODO NUMERICO RK2
!
! Objetivos:
!   Aplicacao do metodo de Runge-Kutta de Segunda Ordem em Dois
!   Estagios (RK2 ou RK22)
!
! Modificado:
!   17 de novembro de 2024
!
! Autoria:
!   oap
! 
MODULE rungekutta2
  USE tipos
  USE rungekutta
  USE funcoes_forca
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_rk2

  TYPE, EXTENDS(integracao) :: integracao_rk2
    
    ! Base do Runge-Kutta
    TYPE(RK) :: baseRK
    ! Modulos adicionais
    CONTAINS
      PROCEDURE :: Iniciar, metodo

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
!   17 de novembro de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE Iniciar (self, massas, G, h, potsoft, E0, J0, corrigir, corme, cormnt, colidir, colmodo, colmd, paralelo)
  IMPLICIT NONE
  CLASS(integracao_rk2), INTENT(INOUT) :: self
  LOGICAL,INTENT(IN) :: corrigir, paralelo
  LOGICAL, INTENT(IN)          :: colidir
  CHARACTER(LEN=*), INTENT(IN) :: colmodo
  REAL(pf), allocatable :: massas(:)
  REAL(pf)              :: G, h, potsoft, colmd, E0, J0(3)
  INTEGER :: a, i
  REAL(pf) :: corme
  INTEGER :: cormnt

  ! quantidade de partículas
  self % N = SIZE(massas)
  ! massas
  ALLOCATE(self % m (self % N))
  self % m = massas

  ! gravidade
  self % G = G
  ! passo
  self % h = h
  ! Softening do potencial
  self % potsoft = potsoft

  ! Se vai ou nao corrigir
  self % corrigir = corrigir
  self % corme = corme
  self % cormnt = cormnt

  ! Valores iniciais
  self % E0 = E0
  self % J0 = J0

  ! Se vai ou nao colidir
  self % colidir = colidir
  self % colisoes_modo = colmodo
  self % colmd = colmd
  
  ! Inicia o base do RK 
  CALL self % baseRK % Iniciar(self % n, self % m, self % G, self % h, self % potsoft)
  
  ! Codigo paralelo
  self % paralelo = paralelo
  IF (paralelo) THEN
    self % forcas_funcao => forcas_par
  ELSE
    self % forcas_funcao => forcas_seq
  ENDIF

END SUBROUTINE Iniciar


! ************************************************************
!! Metodo numerico
!
! Objetivos:
!   Aplicacao do metodo em si.
!
! Modificado:
!   17 de novembro de 2024
!
! Autoria:
!   oap
!
FUNCTION metodo (self, R, P, FSomas_ant)

  IMPLICIT NONE
  class(integracao_rk2), INTENT(IN) :: self
  REAL(pf), DIMENSION(self % N, self % dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self % N, self % dim) :: R1, P1
  REAL(pf), DIMENSION(3, self % N, self % dim) :: metodo

  ! componentes da integracao (kappas)
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_q, k2_q, k3_q, k4_q
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_p, k2_p, k3_p, k4_p

  k1_q = R + 0.5_pf * self % h * P * self % baseRK % massasInvertidas
  R1 = R + self % h * (P + 0.5_pf * self % h * self % forcas (R))
  P1 = P + self % h * self % forcas (k1_q)

  metodo(1,:,:) = R1
  metodo(2,:,:) = P1

END FUNCTION metodo

END module rungekutta2