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

    ! Metodo
    CHARACTER(len=100) :: metodo

    ! Diretorio onde ficara salvo
    CHARACTER(len=256) :: dir

    ! Configuracoes
    LOGICAL :: corrigir=.FALSE., colidir=.FALSE.
    REAL(pf) :: corrigir_margem_erro = 0.1_pf
    INTEGER :: corrigir_max_num_tentativas = 5

    !> Arquivo
    TYPE(arquivo) :: Arq

    CONTAINS
      PROCEDURE :: Iniciar, inicializar_metodo, rodar_verlet, rodar_rk4
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
SUBROUTINE Iniciar (self, G, M, R0, P0, h, passos_antes_salvar, metodo, t0, tf)

  CLASS(simular), INTENT(INOUT) :: self

  REAL(pf), allocatable :: M(:), R0(:,:), P0(:,:)
  REAL(pf) :: G, h
  INTEGER :: a, i, passos_antes_salvar
  CHARACTER(len=*) :: metodo
  INTEGER, OPTIONAL :: t0, tf

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

  ! Salva o metodo
  self % metodo = metodo

  ! Inicializa o metodo
  CALL self % inicializar_metodo()

  ! Copia o arquivo de valores iniciais
  CALL salvar_sorteio('out/data/', self % Arq % dirarq//"/", 'valint.txt', &
  "Sorteio_"//self % Arq % dirarq, &
  self % G,        &
  self % M,        &
  self % R,        &
  self % P,        &
  t0,              &
  tf,              &
  ABS(self % h),   &
  self % metodo,   &
  self % corrigir, &
  self % corrigir_margem_erro, &
  self % corrigir_max_num_tentativas, &
  self % colidir,  &
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
  CALL self % Arq % inicializar_arquivo_info(self%N, self%metodo, self%G, self%h, self%corrigir, self%colidir)

  ! Condicoes iniciais
  CALL self % Arq % escrever((/self % R, self % P/))
END SUBROUTINE inicializar_metodo

! ************************************************************
!! Roda simulacao com Verlet
!
! Objetivos:
!   Roda uma simulacao com o metodo Velocity-Verlet.
!
! Modificado:
!   25 de maio de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE rodar_verlet (self, qntdPassos)
  IMPLICIT NONE
  class(simular), INTENT(INOUT) :: self
  INTEGER, INTENT(IN) :: qntdPassos
  ! Iterador e variavel de tempo
  INTEGER :: i = 0, t
  ! Integrador
  TYPE(integracao_verlet) :: integrador

  REAL(pf), DIMENSION(self % N, self % dim) :: R1, P1
  REAL(pf) :: t0, tf, tempo_total = 0.0_pf
  INTEGER  :: timestep_inv
  
  ! inicializa o integrador 
  CALL integrador % Iniciar(self % M, self % G, self % h, &
    self%corrigir, self%corrigir_margem_erro, self%corrigir_max_num_tentativas, &
    self%colidir)

  ! Condicoes iniciais
  R1 = self % R
  P1 = self % P
  timestep_inv = NINT(1/self % h)

  ! Roda
  WRITE (*, '(a)') '  > iniciando simulacao...'
  DO WHILE (i .le. qntdPassos)
    ! timer
    t0 = omp_get_wtime()
    ! Integracao
    CALL integrador % aplicarNVezes(R1, P1, self % passos_antes_salvar * timestep_inv, self % E0, self % J0)
    
    ! timer
    tf = omp_get_wtime()
    tempo_total = tempo_total + tf - t0

    CALL self % Arq % escrever((/R1, P1/))
    IF (mod(i, 10*self%passos_antes_salvar) == 0) THEN
      WRITE (*,*) '     -> Passo:', i, ' / Energia:', energia_total(self % G, self % M, R1, P1), ' / Tempo: ', tempo_total
      CALL self % Arq % arquivo_bkp(i*self%passos_antes_salvar, tempo_total)
    ENDIF

    i = i + self % passos_antes_salvar
  END DO

  CALL self % Arq % atualizar_arquivo_info(qntdPassos*self%passos_antes_salvar, tempo_total)
  CALL self % Arq % excluir_bkp()

  WRITE (*, '(a)') '  > simulacao encerrada!'
  WRITE (*,*)

  CALL self % Arq % fechar()
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
  
  REAL(pf), DIMENSION(self % N, self % dim) :: R1, P1
  REAL(pf) :: t0, tf, tempo_total = 0.0_pf

  ! inicializa o integrador 
  CALL integrador % Iniciar(self % M, self % G, self % h, &
    self%corrigir, self%corrigir_margem_erro, self%corrigir_max_num_tentativas, &
    self%colidir)

  ! Condicoes iniciais
  R1 = self % R
  P1 = self % P
  WRITE (*,*) self % passos_antes_salvar
  CALL self % Arq % escrever((/R1, P1/))

  ! Roda
  WRITE (*, '(a)') '  > iniciando simulacao...'
  DO WHILE (i .le. qntdPassos)
    ! timer
    t0 = omp_get_wtime()
    ! Integracao
    CALL integrador % aplicarNVezes(R1, P1, self % passos_antes_salvar, self % E0, self % J0)
    
    ! timer
    tf = omp_get_wtime()
    tempo_total = tempo_total + tf - t0

    CALL self % Arq % escrever((/R1, P1/))
    IF (mod(i, 10) == 0) THEN
      WRITE (*,*) '     -> Passo:', i, ' / Energia:', energia_total(self % G, self % M, R1, P1), ' / Tempo: ', tempo_total
      CALL self % Arq % arquivo_bkp(i*self%passos_antes_salvar, tempo_total)
    ENDIF

    i = i + self % passos_antes_salvar
  END DO

  CALL self % Arq % atualizar_arquivo_info(qntdPassos*self%passos_antes_salvar, tempo_total)
  CALL self % Arq % excluir_bkp()

  WRITE (*, '(a)') '  > simulacao encerrada!'
  WRITE (*,*)

END SUBROUTINE rodar_rk4

END module simulacao