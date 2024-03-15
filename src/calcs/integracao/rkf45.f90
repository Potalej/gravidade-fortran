! ************************************************************
!! METODO NUMERICO: RKF45
!(!!! No momento nao esta em operacao !!!)
!
! Objetivos:
!   Aplicacao do metodo de Runge-Kutta-Fehlberg 45 (RKF45).
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
MODULE rkf45
  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  USE rungekutta
  USE correcao
  USE colisao
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_rkf45

  TYPE, EXTENDS(integracao) :: integracao_rkf45

    ! base do runge-kutta
    TYPE(RK) :: baseRK
  
    ! constantes do metodo
    REAL, DIMENSION(6) :: c = (/0, 1/4, 3/8, 12/13, 1, 1/2/)
    REAL, DIMENSION(6) :: a2 = (/1/4, 0, 0, 0, 0, 0/)    
    REAL, DIMENSION(6) :: a3 = (/3/32, 9/32, 0, 0, 0, 0/)
    REAL, DIMENSION(6) :: a4 = (/1932/2197, -7200/2197, 7296/2197, 0, 0, 0/)
    REAL, DIMENSION(6) :: a5 = (/ 439/216, -8, 3680/513, -845/4104, 0, 0 /)
    REAL, DIMENSION(6) :: a6 = (/ -8/27, 2, -3544/2565, 1859/4104, -11/40, 0 /)
    REAL, DIMENSION(6) :: b1 = (/16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55 /) 
    REAL, DIMENSION(6) :: b2 = (/25/216, 0, 1408/2565, 2197/4104, -1/5, 0 /)
    ! REAL, DIMENSION(6) :: CT = (/ 1.0/360.0, 0.0, -128.0/4275.0, -2197.0/75240.0, 1.0/50.0, 2.0/55.0 /)
    REAL(pf),DIMENSION(6)::CT=(/1.0_pf/360.0_pf,0.0_pf,-128.0_pf/4275.0_pf,-2197.0_pf/75240.0_pf,1.0_pf/50.0_pf,2.0_pf/55.0_pf/)
    REAL(pf) :: epsilon = 0.00001_pf

    CONTAINS
      PROCEDURE :: Iniciar, metodo, aplicarNVezes, tolerancia, aplicarNVezesControleAutomatico

  END TYPE

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
  CLASS(integracao_rkf45), INTENT(INOUT) :: self
  REAL(pf), ALLOCATABLE :: massas(:)
  REAL(pf)              :: G, h
  LOGICAL,INTENT(IN) :: corrigir, colidir
  INTEGER :: a, i

  ! quantidade de particulas
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

  ! inicia o metodo
  CALL self % baseRK % Iniciar(self % n, self % m, self % G, self % h)

  ! alocando variaveis de correcao
  ALLOCATE(self%grads(10, 6*self%N))
  ALLOCATE(self%gradsT(6*self%N,10))
  ALLOCATE(self%vetorCorrecao(1:6*self%N))
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
FUNCTION metodo (self, R, P, FSomas, controle)

  IMPLICIT NONE
  CLASS(integracao_rkf45), INTENT(INOUT)                      :: self
  REAL(pf), DIMENSION(self % N, self % dim), INTENT(IN) :: R, P, FSomas
  REAL(pf), DIMENSION(self % N, self % dim)             :: R1, P1
  REAL(pf), DIMENSION(2, self % N, self % dim)          :: metodo
  LOGICAL                                               :: controle

  ! componentes da integracao (kappas)
  REAL(pf), DIMENSION(self % N, self % dim) :: k1, k2, k3, k4, k5, k6, fator
  REAL(pf) :: TE

  DO WHILE (.TRUE.)

    ! metodo RKF45
    k1 = self % h * (P * self % baseRK % massasInvertidas)
    k2 = self % h * (k1 + (self%a2(1)*k1) * self%baseRK%massasInvertidas)
    k3 = self % h * (k1 + (self%a3(1)*k1 + self%a3(2)*k2) * self%baseRK%massasInvertidas)
    k4 = self % h * (k1 + (self%a4(1)*k1 + self%a4(2)*k2 + self%a4(3)*k3) * self%baseRK%massasInvertidas)
    k5 = self % h * (k1 + (self%a5(1)*k1 + self%a5(2)*k2 + self%a5(3)*k3 + self%a5(4)*k4) * self%baseRK%massasInvertidas)
    k6 = k1+(self%a6(1)*k1 + self%a6(2)*k2 + self%a6(3)*k3 + self%a6(4)*k4 + self%a6(5)*k5)*self%baseRK%massasInvertidas
    k6 = self % h * k6

    ! calcula o erro tolerado
    IF (controle) THEN
      CALL self % tolerancia (TE, k1, k2, k3, k4, k5, k6)
    ELSE
      TE = 0.0
    ENDIF
    
    IF (TE > self % epsilon) THEN
      ! WRITE (*,*) 'caiu aqui: ', TE, ' / ', self % h
      self % h = self % h * (self%epsilon/(2*TE))**(0.25)
      self % baseRK % h = self % h
    ELSE
      fator = self%b1(1) * k1 + self%b1(2) * k2 + self%b1(3) * k3 + self%b1(4) + k4 + self%b1(5) * k5 + self%b1(6) * k6
      
      ! integra as posicoes
      R1 = R + fator
      ! integra os momentos
      P1 = P + self % h * FSomas

      metodo(1,:,:) = R1
      metodo(2,:,:) = P1

      IF (TE > 0) THEN
        self % h = 0.001_pf
        self % baseRK % h = self % h
      ENDIF
      
      EXIT
    ENDIF 
  
  END DO

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
SUBROUTINE aplicarNVezes (self, R, P, passos, E0, J0)

  IMPLICIT NONE
  CLASS (integracao_rkf45), INTENT(INOUT)                    :: self
  REAL(pf), DIMENSION(self % N, self % dim), INTENT(INOUT) :: R, P
  INTEGER, INTENT(IN)                               :: passos
  REAL(pf), INTENT(IN)                                  :: E0
  REAL(pf), DIMENSION(3), INTENT(IN)                    :: J0
  ! para cada passo
  INTEGER :: i
  ! para as forcas e passos pos-integracao
  REAL(pf), DIMENSION (self % N, self % dim) :: F, R1, P1
  REAL(pf), DIMENSION (2, self % N, self % dim) :: resultado 
  ! para verificar se corrigiu
  LOGICAL :: corrigiu = .FALSE.

  R1 = R
  P1 = P

  DO i = 1, passos
    ! calcula as forcas
    F = self % baseRK % forcas (R1)
    ! aplicada o metodo
    resultado = self % metodo (R1, P1, F, .FALSE.)

    R1 = resultado(1,:,:)
    P1 = resultado(2,:,:)

    ! aplica a correcao geral
    CALL corrigir(self % G, self % m, R1, P1,self%grads,self%gradsT,self%vetorCorrecao, corrigiu, E0, J0)

    IF (.NOT. corrigiu) THEN
      CALL verificar_e_colidir (self % m, R1, P1)
    ENDIF

  END DO

  R = R1
  P = P1

END SUBROUTINE aplicarNVezes

! ************************************************************
!! Aplicacao iterada do metodo com controle de passos
!
! Objetivos:
!   Aplica o metodo iterativamente N vezes com o controle de
!   passo.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
SUBROUTINE aplicarNVezesControleAutomatico (self, R, P, passos, E0, J0, h0)
  
  IMPLICIT NONE
  CLASS(integracao_rkf45), INTENT(INOUT)                    :: self
  REAL(pf), DIMENSION(self % N, self % dim), INTENT(INOUT) :: R, P
  INTEGER, INTENT(IN)                               :: passos
  REAL(pf), INTENT(IN)                                  :: E0, h0
  REAL(pf), DIMENSION(3), INTENT(IN)                    :: J0
  ! para cada passo
  INTEGER :: i
  ! para as forças e passos pós-integração
  REAL(pf), DIMENSION (self % N, self % dim) :: F, R1, P1
  REAL(pf), DIMENSION (2, self % N, self % dim) :: resultado
  REAL(pf) :: h_soma
  
  R1 = R
  P1 = P

  DO i = 1, passos
    
    h_soma = 0.0_pf

    loop_while: DO WHILE (h_soma < h0)

      ! calcula as forcas
      F = self % baseRK % forcas (R1)
      ! aplicada o metodo
      resultado = self % metodo (R1, P1, F, .TRUE.)

      R1 = resultado(1,:,:)
      P1 = resultado(2,:,:)
      
      h_soma = h_soma + self % h

    END DO loop_while

  END DO

  R = R1
  P = P1

END SUBROUTINE aplicarNVezesControleAutomatico

! ************************************************************
!! Tolerancia
!
! Objetivos:
!   Calculo da tolerancia do metodo.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE tolerancia (self, TE, k1, k2, k3, k4, k5, k6)

  IMPLICIT NONE
  CLASS(integracao_rkf45), INTENT(IN) :: self
  REAL(pf), DIMENSION(self % N, self % dim) :: TE_soma, k1, k2, k3, k4, k5, k6
  REAL(pf), DIMENSION(3) :: TE_vet
  REAL(pf) :: TE

  TE_soma = k1*self%CT(1) + k2*self%CT(2) + k3*self%CT(3) + k4*self%CT(4) + k5*self%CT(5) + k6*self%CT(6)
  TE_vet = PRODUCT(TE_soma, DIM = 1)**(0.5)
  TE = MAXVAL(TE_vet)

END SUBROUTINE tolerancia

END MODULE rkf45