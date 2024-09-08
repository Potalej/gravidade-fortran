! ************************************************************
!! METODO NUMERICO: Velocity-Verlet
!
! Objetivos:
!   Aplicacao do metodo simpletico de Velocity-Verlet.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
MODULE verlet
  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  USE OMP_LIB
  USE mecanica
  USE correcao
  USE colisao
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_verlet

  TYPE, EXTENDS(integracao) :: integracao_verlet

    CONTAINS
      PROCEDURE :: Iniciar, metodo, aplicarNVezes, Forcas

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
SUBROUTINE Iniciar (self, massas, G, h, corrigir, corme, cormnt, colidir)
  IMPLICIT NONE
  class(integracao_verlet), INTENT(INOUT) :: self
  REAL(pf), allocatable :: massas(:)
  REAL(pf)              :: G, h
  LOGICAL,INTENT(IN) :: corrigir, colidir
  REAL(pf) :: corme
  INTEGER :: cormnt
  INTEGER :: a, i

  ! Quantidade de particulas
  self % N = SIZE(massas)
  ! Massas
  ALLOCATE(self % m (self % N))
  self % m = massas

  ! gravidade
  self % G = G
  ! Passo
  self % h = h

  ! Se vai ou nao corrigir
  self % corrigir = corrigir
  self % corme = corme
  self % cormnt = cormnt

  ! Se vai ou nao colidir
  self % colidir = colidir  

  ! Alocando variaveis de correcao
  ALLOCATE(self%grads(4, 6*self%N))
  ALLOCATE(self%gradsT(6*self%N,4))
  ALLOCATE(self%vetorCorrecao(1:6*self%N))

END SUBROUTINE Iniciar

! ************************************************************
!! Forcas
!
! Objetivos:
!   Calculo das forcas.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
FUNCTION forcas (self, R)
  IMPLICIT NONE
  class(integracao_verlet), INTENT(IN) :: self
  REAL(pf), DIMENSION(self % N, self % dim), INTENT(IN) :: R
  REAL(pf), DIMENSION(self % dim) :: Fab, dif
  INTEGER :: a, b, thread, threads, qntdPorThread
  REAL(pf) :: distancia
  REAL(pf), DIMENSION(self % N, self % dim) :: forcas, forcas_local
  
  forcas(:,:) = 0

  !$OMP PARALLEL SHARED(forcas) PRIVATE(forcas_local, Fab)
    forcas(:,:) = 0
    forcas_local(:,:) = 0
    !$OMP DO
    DO a = 2, self%N
      DO b = 1, a-1
        ! distancia entre os corpos
        distancia = norm2(R(b,:) - R(a,:))**3
        ! forca entre os corpos a e b
        Fab = self % G * self % m(a) * self % m(b) * (R(b,:) - R(a,:))/distancia
        ! Adiciona na matriz
        forcas_local(a,:) = forcas_local(a,:) + Fab
        forcas_local(b,:) = forcas_local(b,:) - Fab
      END DO
    END DO
    !$OMP END DO

    !$OMP CRITICAL
      forcas = forcas + forcas_local
    !$OMP END CRITICAL
  !$OMP END PARALLEL

END FUNCTION forcas

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
FUNCTION metodo (self, R, P, FSomas_ant)

  IMPLICIT NONE
  class(integracao_verlet), INTENT(IN) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, FSomas_prox
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo
  INTEGER :: a

  ! Integrando as posicoes
  DO a = 1, self % N
    R1(a,:) = R(a,:) + self % h * P(a,:) / self % m(a) + 0.5*(self%h**2)*FSomas_ant(a,:)/self%m(a)
  END DO

  ! Calcula as novas forcas
  FSomas_prox = self%forcas(R1)

  ! Integrando as velocidades
  P1 = P + 0.5*self%h*(FSomas_ant + FSomas_prox)

  metodo(1,:,:) = R1
  metodo(2,:,:) = P1
  metodo(3,:,:) = FSomas_prox

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
  class (integracao_verlet), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P
  INTEGER, INTENT(IN) :: passos_antes_salvar
  REAL(pf), INTENT(IN) :: E0
  REAL(pf)             :: E
  REAL(pf), DIMENSION(3), INTENT(IN) :: J0
  ! Para cada passo
  INTEGER :: i
  ! para verificar se corrigiu
  LOGICAL :: corrigiu = .FALSE.
  ! Para as forcas e passos pos-integracao
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, FSomas_ant
  REAL(pf), DIMENSION(3, self%N, self%dim) :: resultado
  ! Energia total aproximada
  REAL(pf) :: Et_aprox

  ! Salvando as primeiras posicoes e momentos
  R1 = R
  P1 = P

  ! Calcula as forcas
  FSomas_ant = self%forcas(R)

  ! Integrando
  DO i = 1, passos_antes_salvar
    ! Aplica o metodo
    resultado = self % metodo(R1, P1, FSomas_ant)

    ! se tiver colisoes, aplica
    IF (self % colidir .AND. .NOT. corrigiu) THEN
      CALL verificar_e_colidir(self % m, R1, P1)
    ENDIF

    R1 = resultado(1,:,:)
    P1 = resultado(2,:,:)
    FSomas_ant = resultado(3,:,:)

  END DO

  ! Se estiver disposto a corrigir, calcula a energia total para ver se precisa
  IF (self%corrigir) THEN
    E = energia_total(self % G, self % m, R1, P1)
    IF (ABS(E - E0) > self%corme) THEN
      CALL corrigir(self%corme,self%cormnt,self % G, self % m, R1, P1,self%grads,self%gradsT,self%vetorCorrecao, corrigiu, E0, J0)
    END IF
  ENDIF

  R = R1
  P = P1

END SUBROUTINE aplicarNVezes

END MODULE verlet