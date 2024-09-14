! ************************************************************
!! METODO NUMERICO: Ruth 3a Ordem
!
! Objetivos:
!   Aplicacao do metodo simpletico de Ruth com 3a Ordem
!
! Modificado:
!   14 de setembro de 2024
!
! Autoria:
!   oap
!
MODULE ruth3
  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  USE OMP_LIB
  USE mecanica
  USE correcao
  USE colisao
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_ruth3

  TYPE, EXTENDS(integracao) :: integracao_ruth3

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
!   14 de setembro de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE Iniciar (self, massas, G, h, potsoft, corrigir, corme, cormnt, colidir, colmd)
  IMPLICIT NONE
  class(integracao_ruth3), INTENT(INOUT) :: self
  REAL(pf), allocatable :: massas(:)
  REAL(pf)              :: G, h
  LOGICAL,INTENT(IN) :: corrigir, colidir
  REAL(pf) :: corme, potsoft, colmd
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
  ! Softening do potencial
  self % potsoft = potsoft

  ! Se vai ou nao corrigir
  self % corrigir = corrigir
  self % corme = corme
  self % cormnt = cormnt

  ! Se vai ou nao colidir
  self % colidir = colidir
  self % colmd = colmd

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
!   14 de setembro de 2024
!
! Autoria:
!   oap
!
FUNCTION forcas (self, R)
  IMPLICIT NONE
  class(integracao_ruth3), INTENT(IN) :: self
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
        IF (self % potsoft .NE. 0) THEN
          distancia = (norm2(R(b,:) - R(a,:))**2 + self%potsoft**2)**(3/2)
        ELSE
          distancia = norm2(R(b,:) - R(a,:))**3
        ENDIF
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
!   14 de setembro de 2024
!
! Autoria:
!   oap
!
FUNCTION metodo (self, R, P, FSomas_ant)

  IMPLICIT NONE
  class(integracao_ruth3), INTENT(IN) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, FSomas_prox
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo
  INTEGER :: a

  ! i = 3

  ! Momentos: p_ = p + 7/24 F(q) h
  P1 = P + 0.29166666666666666_pf*self%h*FSomas_ant
  
  ! Posicoes: q_ = q + 2/3 p_/m h
  DO a = 1, self % N
    R1(a,:) = R(a,:) + (2/3) * self % h * P1(a,:) / self % m(a)
  END DO


  ! i = 2
  FSomas_prox = self%forcas(R1)

  ! Momentos: p_ = p_ + 3/4 F(q_) h
  P1 = P1 + 0.75_pf * FSomas_prox * self % h

  ! Posicoes: q_ = q_ - 2/3 p_/m h
  DO a = 1, self % N
    R1(a,:) = R1(a,:) - (2/3) * self % h * P1(a,:) / self % m(a)
  END DO


  ! i = 1
  FSomas_prox = self%forcas(R1)

  ! Momentos: p_ = p_ - 1/24 F(q_) h
  P1 = P1 - 0.041666666666666666_pf * FSomas_prox * self % h

  ! Posicoes: q_ = q_ + 1 p_/m h
  DO a = 1, self % N
    R1(a,:) = R1(a,:) + self % h * P1(a,:) / self % m(a)
  END DO


  ! Calcula as novas forcas
  FSomas_prox = self%forcas(R1)


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
!   14 de setembro de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE aplicarNVezes (self, R, P, passos_antes_salvar, E0, J0)

  IMPLICIT NONE
  class (integracao_ruth3), INTENT(INOUT) :: self
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

    R1 = resultado(1,:,:)
    P1 = resultado(2,:,:)
    FSomas_ant = resultado(3,:,:)

    ! se tiver colisoes, aplica
    IF (self % colidir) THEN
      CALL verificar_e_colidir(self % m, R1, P1, self % colmd)
    ENDIF
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

END MODULE ruth3