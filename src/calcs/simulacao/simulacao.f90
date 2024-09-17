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
  USE integrador
  USE rungekutta4
  USE verlet
  USE eulersimp
  USE ruth3
  USE ruth4
  USE rkn551
  USE rkn671
  USE svcp8s15
  USE svcp10s35
  ! Para configuracoes
  USE leitura

  IMPLICIT NONE
  PRIVATE
  PUBLIC simular

  ! Tipos de metodo de integracao
  TYPE(integracao_verlet),    TARGET, SAVE :: INT_VERLET
  TYPE(integracao_eulersimp), TARGET, SAVE :: INT_EULERSIMP
  TYPE(integracao_rk4),       TARGET, SAVE :: INT_RK4
  TYPE(integracao_ruth3),     TARGET, SAVE :: INT_RUTH3
  TYPE(integracao_ruth4),     TARGET, SAVE :: INT_RUTH4
  TYPE(integracao_rkn551),    TARGET, SAVE :: INT_RKN551
  TYPE(integracao_rkn671),    TARGET, SAVE :: INT_RKN671
  TYPE(integracao_svcp8s15),  TARGET, SAVE :: INT_SVCP8S15
  TYPE(integracao_svcp10s35),  TARGET, SAVE :: INT_SVCP10S35

  ! Classe de simulacao
  TYPE :: simular
  
    !> N: Quantidade de corpos
    !> dim: Dimensao do problema
    INTEGER :: N, dim = 3, passos_antes_salvar, t0, tf

    !> h: Tamanho do passo de integracao
    !> G: Constante de gravitacao universal
    !> E0: Energia total inicial
    !> mtot: Massa total do sistema
    !> potsoft: Softening do potencial
    REAL(pf) :: h, G, E0, mtot, potsoft
    
    !> M: Massas do sistema
    !> R: Posicoes das particulas
    !> P: Momento linear das particulas
    !> Jtot: Momento angular total do sistema
    !> Ptot: Momento linear total do sistema
    !> Rcm: Centro de massas do sistema
    REAL(pf), allocatable :: M(:), R(:,:), P(:,:), Jtot(:), Ptot(:), Rcm(:)
    
    REAL(pf), DIMENSION(3) :: J0

    ! Metodo
    CHARACTER(len=100) :: metodo

    ! Diretorio onde ficara salvo
    CHARACTER(len=256) :: dir

    ! Configuracoes
    LOGICAL :: corrigir=.FALSE., colidir=.FALSE.
    REAL(pf) :: corrigir_margem_erro = 0.1_pf
    INTEGER :: corrigir_max_num_tentativas = 5
    REAL(pf) :: colisoes_max_distancia = 0.1

    !> Arquivo
    TYPE(arquivo) :: Arq

    CONTAINS
      PROCEDURE :: Iniciar, inicializar_metodo, rodar
      
  END TYPE

CONTAINS
  
! ************************************************************
!! Metodo principal
!
! Objetivos:
!   Faz a simulacao.
!
! Modificado:
!   26 de maio de 2024
!
! Autoria:
!   oap
! 
! SUBROUTINE Iniciar (self, G, M, R0, P0, h, potsoft, passos_antes_salvar, metodo, t0, tf)
SUBROUTINE Iniciar (self, configs, M, R0, P0, h)

  CLASS(simular), INTENT(INOUT) :: self
  CLASS(preset_config) :: configs
  REAL(pf), allocatable :: M(:), R0(:,:), P0(:,:)
  REAL(pf) :: h
  INTEGER :: a, i

  self % passos_antes_salvar = configs % passos_antes_salvar
  
  ! Salva o tamanho dos passos
  self % h = h

  ! Salva o softening do potencial
  self % potsoft = configs % potsoft

  ! Salva a gravidade
  self % G = configs % G

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
  self % E0 = energia_total(self % G, self % M, self % R, self % P)

  ! Salva o momento angular inicial
  self % J0 = momento_angular_total(self % R, self % P)

  ! Salva o centro de massas inicial
  self % Rcm = centro_massas(self % M, self % R)

  ! Salva o metodo
  self % metodo = configs % integrador

  self % t0 = configs % t0
  self % tf = configs % tf

  ! Salva o corretor
  self % corrigir = configs%corretor
  self % corrigir_margem_erro = configs%corretor_margem_erro
  self % corrigir_max_num_tentativas = configs%corretor_max_num_tentativas
  
  ! Salva as colisoes
  self % colidir  = configs%colisoes
  self % colisoes_max_distancia = configs%colisoes_max_distancia

  ! Inicializa o metodo
  CALL self % inicializar_metodo()

  ! Copia o arquivo de valores iniciais
  CALL salvar_sorteio('out/data/', self % Arq % dirarq//"/", 'valint.txt', &
  "Sorteio_"//self % Arq % dirarq, &
  self % G,        &
  self % M,        &
  self % R,        &
  self % P,        &
  configs % t0,    &
  configs % tf,    &
  ABS(self % h),   &
  self % potsoft,  &
  self % metodo,   &
  self % corrigir, &
  self % corrigir_margem_erro, &
  self % corrigir_max_num_tentativas, &
  self % colidir,  &
  self % colisoes_max_distancia,  &
  self % passos_antes_salvar)

END SUBROUTINE Iniciar

SUBROUTINE inicializar_metodo (self)
  IMPLICIT NONE
  class(simular), INTENT(INOUT) :: self

  ! Cria o arquivo onde sera salvo
  CALL self % Arq % criar(self % N, self % dim)
  self % dir = self % Arq % dirarq

  ! Salva as infos de cabecalho
  CALL self % Arq % escrever_cabecalho(self % h, self % G, self % M)

  ! Salva as informacoes no info.txt
  CALL self % Arq % inicializar_arquivo_info(self%N, self%metodo, self%G, self%h, self%potsoft, &
    self%t0, self%tf, self % passos_antes_salvar, &
    self%corrigir, self%corrigir_margem_erro, self%corrigir_max_num_tentativas, &
    self%colidir, self%colisoes_max_distancia)

  ! Condicoes iniciais
  CALL self % Arq % escrever((/self % R, self % P/))

  WRITE(*,*) 'E0 = ', self % E0
END SUBROUTINE inicializar_metodo


! ************************************************************
!! Roda simulacao com algum metodo definido dinamicamente
!
! Objetivos:
!   Roda uma simulacao com algum metodo numerico desejado.
!
! Modificado:
!   14 de setembro de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE rodar (self, qntdPassos)
  IMPLICIT NONE
  class(simular), INTENT(INOUT) :: self
  INTEGER, INTENT(IN) :: qntdPassos  
  ! Iterador e variavel de tempo
  INTEGER :: i = 0, t
  
  ! Integrador, definido dinamicamente
  CLASS(integracao), POINTER :: integrador
  
  ! Variaveis locais
  REAL(pf), DIMENSION(self % N, self % dim) :: R1, P1
  REAL(pf) :: t0, tf, tempo_total = 0.0_pf
  INTEGER  :: timestep_inv

  ! Definicao dinamica do metodo
  WRITE (*,'(a)') 'METODO: ', self % metodo
  
  SELECT CASE (TRIM(self % metodo))
    CASE ('verlet');    integrador => INT_VERLET
    CASE ('rk4');       integrador => INT_RK4
    CASE ('eulersimp'); integrador => INT_EULERSIMP
    CASE ('ruth3');     integrador => INT_RUTH3
    CASE ('ruth4');     integrador => INT_RUTH4
    CASE ('rkn551');    integrador => INT_RKN551
    CASE ('rkn671');    integrador => INT_RKN671
    CASE ('svcp8s15');  integrador => INT_SVCP8S15
    CASE ('svcp10s35');  integrador => INT_SVCP10S35
    
    ! Por padrao, sera o VERLET
    CASE DEFAULT;       WRITE(*,*) 'Metodo nao identificado!'
  END SELECT

  ! inicializa o integrador 
  CALL integrador % Iniciar(self % M, self % G, self % h, self % potsoft, &
    self%corrigir, self%corrigir_margem_erro, self%corrigir_max_num_tentativas, &
    self%colidir, self%colisoes_max_distancia)

  ! Condicoes iniciais
  R1 = self % R
  P1 = self % P
  CALL self % Arq % escrever((/R1, P1/))
  timestep_inv = NINT(1/self % h)

  ! Roda
  WRITE (*, '(a)') '  > iniciando simulacao...'
  WRITE(*,*) qntdPassos * timestep_inv
  DO WHILE (i .le. qntdPassos * timestep_inv)
    ! timer
    t0 = omp_get_wtime()
    ! Integracao
    CALL integrador % aplicarNVezes(R1, P1, self % passos_antes_salvar, self % E0, self % J0)
    
    ! timer
    tf = omp_get_wtime()
    tempo_total = tempo_total + tf - t0

    CALL self % Arq % escrever((/R1, P1/))
    IF (mod(i, 10*self%passos_antes_salvar) == 0) THEN
      WRITE (*,*) '     -> Passo:', i, ' / Energia:', energia_total(self % G, self % M, R1, P1), ' / Tempo: ', tempo_total
      CALL self % Arq % arquivo_bkp(i*self%passos_antes_salvar*timestep_inv, tempo_total)
    ENDIF

    i = i + self % passos_antes_salvar
  END DO

  CALL self % Arq % atualizar_arquivo_info(qntdPassos*self%passos_antes_salvar*timestep_inv, tempo_total)
  CALL self % Arq % excluir_bkp()

  WRITE (*, '(a)') '  > simulacao encerrada!'
  WRITE (*,*)

  CALL self % Arq % fechar()

END SUBROUTINE rodar

END module simulacao