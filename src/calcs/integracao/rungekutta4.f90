! ************************************************************
!! METODO NUMERICO RK4
!
! Objetivos:
!   Aplicacao do metodo de Runge-Kutta de Quarta Ordem em Quatro
!   Estagios (RK4 ou RK44)
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
MODULE rungekutta4
  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  USE rungekutta
  USE correcao
  USE colisao
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_rk4

  TYPE, EXTENDS(integracao) :: integracao_rk4
    
    ! Base do Runge-Kutta
    TYPE(RK) :: baseRK
    ! Modulos adicionais
    CONTAINS
      PROCEDURE :: Iniciar, metodo, aplicarNVezes

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
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE Iniciar (self, massas, G, h, corrigir, colidir)
  IMPLICIT NONE
  CLASS(integracao_rk4), INTENT(INOUT) :: self
  LOGICAL,INTENT(IN) :: corrigir, colidir
  REAL(pf), allocatable :: massas(:)
  REAL(pf)              :: G, h
  INTEGER :: a, i

  ! quantidade de partículas
  self % N = SIZE(massas)
  ! massas
  ALLOCATE(self % m (self % N))
  self % m = massas

  ! gravidade
  self % G = G
  ! passo
  self % h = h

  ! Se vai ou nao corrigir
  self % corrigir = corrigir

  ! Se vai ou nao colidir
  self % colidir = colidir

  ! alocando variaveis de correcao
  ALLOCATE(self%grads(10, 6*self%N))
  ALLOCATE(self%gradsT(6*self%N,10))
  ALLOCATE(self%vetorCorrecao(1:6*self%N))
  
  ! Inicia o base do RK 
  CALL self % baseRK % Iniciar(self % n, self % m, self % G, self % h)

END SUBROUTINE Iniciar


! ************************************************************
!! Metodo numerico
!
! Objetivos:
!   Aplicacao do metodo em si.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
FUNCTION metodo (self, R, P, FSomas)

  IMPLICIT NONE
  class(integracao_rk4), INTENT(IN) :: self
  REAL(pf), DIMENSION(self % N, self % dim), INTENT(IN) :: R, P, FSomas
  REAL(pf), DIMENSION(self % N, self % dim) :: R1, P1
  REAL(pf), DIMENSION(2, self % N, self % dim) :: metodo

  ! componentes da integracao (kappas)
  REAL(pf), DIMENSION(self % N, self % dim) :: k1, k2, k3, k4, fator

  ! faz a integracao sobre as equacoes x'
  k1 = P * self % baseRK % massasInvertidas
  k2 = k1 * self % baseRK % massasInvertidas
  k3 = k2 * self % baseRK % massasInvertidas
  k4 = k3 * self % baseRK % massasInvertidas

  ! fator para integracao
  fator = (self % h / 6) * (6*k1 + 3*self % h*k2 + self % h**2 * k3 + 0.25 * self % h**3 * k4)

  ! integra as posicoes
  R1 = R + fator

  ! integra os momentos
  P1 = P + self % h * FSomas

  metodo(1,:,:) = R1
  metodo(2,:,:) = P1

END FUNCTION metodo

! ************************************************************
!! Aplicacao iterada do metodo
!
! Objetivos:
!   Aplica o metodo iterativamente N vezes.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE aplicarNVezes (self, R, P, passos_antes_salvar, E0, J0)

  IMPLICIT NONE
  class (integracao_rk4), INTENT(INOUT)                    :: self
  REAL(pf), DIMENSION(self % N, self % dim), INTENT(INOUT) :: R, P
  INTEGER, INTENT(IN)                                   :: passos_antes_salvar
  REAL(pf), INTENT(IN)                                  :: E0
  REAL(pf), DIMENSION(3), INTENT(IN)                    :: J0
  ! para cada passo
  INTEGER :: i, qntd = 10
  ! para verificar se corrigiu
  LOGICAL :: corrigiu = .FALSE.
  ! para as forcas e passos pos-integracao
  REAL(pf), DIMENSION (self % N, self % dim) :: F, R1, P1
  REAL(pf), DIMENSION (2, self % N, self % dim) :: resultado  
  R1 = R
  P1 = P

  DO i = 1, passos_antes_salvar
    ! calcula as forcas
    F = self % baseRK % forcas (R1)
    ! aplicada o metodo
    resultado = self % metodo (R1, P1, F)
    
    R1 = resultado(1,:,:)
    P1 = resultado(2,:,:)

    ! aplica a correcao geral, se solicitado
    IF (self % corrigir) THEN
      CALL corrigir(self % G, self % m, R1, P1,self%grads,self%gradsT,self%vetorCorrecao, corrigiu, E0, J0)
    ENDIF

    ! se tiver colisoes, aplica
    IF (self % colidir .AND. .NOT. corrigiu) THEN
      CALL verificar_e_colidir(self % m, R1, P1)
    ENDIF

  END DO

  R = R1
  P = P1

END SUBROUTINE aplicarNVezes


END module rungekutta4