! ************************************************************
!! SIMULACAO: GERAL
!
! Objetivos:
!   Arquivo base para fazer simulacoes.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
MODULE simulacao
  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  USE OMP_LIB

  ! Auxiliares
  USE mecanica
  USE auxiliares
  USE arquivos
  ! Metodos de integracao numerica
  USE rungekutta4
  USE verlet

  IMPLICIT NONE
  PRIVATE
  PUBLIC simular

  ! Classe de simulacao
  TYPE :: simular
  
    !> N: Quantidade de corpos
    !> dim: Dimensao do problema
    INTEGER :: N, dim = 3, passos_antes_salvar

    !> h: Tamanho do passo de integracao
    !> G: Constante de gravitacao universal
    !> E0: Energia total inicial
    !> mtot: Massa total do sistema
    REAL(pf) :: h, G, E0, mtot
    
    !> M: Massas do sistema
    !> R: Posicoes das particulas
    !> P: Momento linear das particulas
    !> Jtot: Momento angular total do sistema
    !> Ptot: Momento linear total do sistema
    !> Rcm: Centro de massas do sistema
    REAL(pf), allocatable :: M(:), R(:,:), P(:,:), Jtot(:), Ptot(:), Rcm(:)
    
    REAL(pf), DIMENSION(3) :: J0

    ! Configuracoes
    LOGICAL :: corrigir=.FALSE., colidir=.FALSE.

    !> Arquivo
    TYPE(arquivo) :: Arq

    CONTAINS
      PROCEDURE :: Iniciar, rodar_verlet, rodar_rk4
  END TYPE

CONTAINS
  
! ************************************************************
!! Metodo principal
!
! Objetivos:
!   Faz a simulacao.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE Iniciar (self, G, M, R0, P0, h, passos_antes_salvar)

  CLASS(simular), INTENT(INOUT) :: self

  REAL(pf), allocatable :: M(:), R0(:,:), P0(:,:)
  REAL(pf) :: G, h
  INTEGER :: a, i, passos_antes_salvar

  self % passos_antes_salvar = passos_antes_salvar
  
  ! Salva o tamanho dos passos
  self % h = h

  ! Salva a gravidade
  self % G = G

  ! Salva as massas
  self % M = M  
  
  ! Quantidade corpos no sistema
  self % N = SIZE(M)

  ! Massa total do sistema
  self % mtot = SUM(self % M)

  ! Salva as posicoes e os momentos
  self % R = R0
  self % P = P0

  ! Salva a dimensao
  self % dim = SIZE(R0,2)

  ! Salva momento linear total inicial
  self % Ptot = momentoLinear_total(self % P)

  ! Salva a energia inicial
  self % E0 = energia_total(G, self % M, self % R, self % P)

  ! Salva o momento angular inicial
  self % J0 = momento_angular_total(self % R, self % P)

  ! Salva o centro de massas inicial
  self % Rcm = centro_massas(self % M, self % R)

END SUBROUTINE Iniciar

! ************************************************************
!! Roda simulacao com Verlet
!
! Objetivos:
!   Roda uma simulacao com o metodo Velocity-Verlet.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE rodar_verlet (self, qntdPassos)
  IMPLICIT NONE
  class(simular), INTENT(INOUT) :: self
  INTEGER, INTENT(IN) :: qntdPassos
  ! iterador e variavel de tempo que sera o nome do arquivo
  INTEGER :: i, t
  ! Integrador
  TYPE(integracao_verlet) :: integrador
  ! Escritor de arquivos
  TYPE(arquivo) :: Arq
  
  REAL(pf), DIMENSION(2, self % N, self % dim) :: resultado
  REAL(pf), DIMENSION(self % N, self % dim) :: R1, P1

  REAL(pf) :: t0, tf, tempo_total = 0.0_pf

  ! Cria o arquivo onde sera salvo
  CALL Arq % criar(2, self % N, self % dim)

  ! Salva as infos de cabecalho
  CALL Arq % escrever_cabecalho(self % h, self % G, self % M)

  ! Instanciamento do integrador
  WRITE (*,'(a)') 'INTEGRACAO NUMERICA'
  WRITE (*,*) ' > rodando com ', qntdPassos, ' passos'
  WRITE (*, '(a)') '  > instanciando o metodo velocity-verlet'
  CALL integrador % Iniciar(self % M, self % G, self % h, self%corrigir, self%colidir)

  ! Condicoes iniciais
  R1 = self % R
  P1 = self % P

  CALL Arq % escrever((/R1, P1/))

  ! Roda
  WRITE (*, '(a)') '  > iniciando simulacao...'
  DO i = 1, qntdPassos

    ! timer
    t0 = omp_get_wtime()

    ! Integracao
    CALL integrador % aplicarNVezes(R1, P1, self % passos_antes_salvar, self % E0, self % J0)

    ! timer
    tf = omp_get_wtime()
    tempo_total = tempo_total + tf - t0

    CALL Arq % escrever((/R1, P1/))

    IF (mod(i, 500) == 0) THEN
      WRITE (*,*) '     -> Passo:', i, ' / Energia:', energia_total(self % G, self % M, R1, P1)
    ENDIF

  END do
  
  WRITE (*, '(a)') '  > simulacao encerrada!'
  WRITE (*,*)

  CALL Arq % fechar()

END SUBROUTINE rodar_verlet

! ************************************************************
!! Roda simulacao com RK4
!
! Objetivos:
!   Roda uma simulacao com o metodo de Runge-Kutta de ordem 4.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE rodar_rk4 (self, qntdPassos)

  IMPLICIT NONE
  class(simular), INTENT(INOUT) :: self
  INTEGER, INTENT(IN) :: qntdPassos
  ! Iterador e variavel de tempo que sera o nome do arquivo
  INTEGER :: i, t
  ! Integrador
  TYPE(integracao_rk4) :: integrador
  ! Escritor de arquivos
  TYPE(arquivo) :: Arq
  
  REAL(pf), DIMENSION(2, self % N, self % dim) :: resultado
  REAL(pf), DIMENSION(self % N, self % dim) :: R1, P1

  ! Instanciamento do integrador    
  CALL integrador % Iniciar(self % M, self % G, self % h, self%corrigir, self%colidir)

  ! Cria o arquivo onde ficara salvo
  CALL Arq % criar(1, self % N, self % dim)

  ! Salva as infos de cabecalho
  CALL Arq % escrever_cabecalho(self % h, self % G, self % M)

  ! Condicoes iniciais
  R1 = self % R
  P1 = self % P

  CALL Arq % escrever((/R1, P1/))

  ! Roda
  DO i = 1, qntdPassos

    CALL integrador % aplicarNVezes(R1, P1, self % passos_antes_salvar, self % E0, self % J0)

    CALL Arq % escrever((/R1, P1/))

    IF (mod(i, 500) == 0) THEN
      WRITE (*,*) 'energia:', energia_total(self % G, self % M, R1, P1)
      WRITE (*,*) 'passo: ', i
    ENDIF

  END DO 

  CALL Arq % fechar()

END SUBROUTINE rodar_rk4

END module simulacao