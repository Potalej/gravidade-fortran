! ************************************************************
!! METODO NUMERICO RK2
!
! Objetivos:
!   Aplicacao do metodo de Runge-Kutta de Terceira Ordem em Tres
!   Estagios (RK3 ou RK33)
!
! Modificado:
!   17 de novembro de 2024
!
! Autoria:
!   oap
! 
MODULE rungekutta3
  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  USE rungekutta
  USE funcoes_forca
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_rk3

  TYPE, EXTENDS(integracao) :: integracao_rk3
    
    ! Base do Runge-Kutta
    TYPE(RK) :: baseRK
    ! Modulos adicionais
    CONTAINS
      PROCEDURE :: Iniciar, metodo

  END TYPE integracao_rk3

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
SUBROUTINE Iniciar (self, massas, G, h, potsoft, corrigir, corme, cormnt, colidir, colmd, paralelo)
  IMPLICIT NONE
  CLASS(integracao_rk3), INTENT(INOUT) :: self
  LOGICAL,INTENT(IN) :: corrigir, colidir, paralelo
  REAL(pf), allocatable :: massas(:)
  REAL(pf)              :: G, h, potsoft, colmd
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

  ! Se vai ou nao colidir
  self % colidir = colidir
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
  class(integracao_rk3), INTENT(IN) :: self
  REAL(pf), DIMENSION(self % N, self % dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self % N, self % dim) :: R1, P1
  REAL(pf), DIMENSION(3, self % N, self % dim) :: metodo

  ! componentes da integracao (kappas)
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_q, k2_q, k3_q
  REAL(pf), DIMENSION(self % N, self % dim) :: k1_p, k2_p, k3_p

  k1_q = P * self % baseRK % massasInvertidas
  k1_p = self % forcas (R)

  k2_q = (P + 0.5 * self % h * k1_p) * self % baseRK % massasInvertidas
  k2_p = self % forcas (R + 0.5 * self % h * k1_q)

  k3_q = (P - self % h * k1_p + 2.0_pf * self % h * k2_p) * self % baseRK % massasInvertidas
  k3_p = self % forcas (R - self % h * k1_q + 2.0_pf * self % h * k2_q)

  ! fator para integracao
  R1 = R + (self % h / 6) * (k1_q + 4.0_pf * k2_q + k3_q)
  P1 = P + (self % h / 6) * (k1_p + 4.0_pf * k2_p + k3_p)

  metodo(1,:,:) = R1
  metodo(2,:,:) = P1

END FUNCTION metodo

END module rungekutta3