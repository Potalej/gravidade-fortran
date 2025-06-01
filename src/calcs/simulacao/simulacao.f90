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
  USE tipos
  USE OMP_LIB

  ! Auxiliares
  USE mecanica
  USE auxiliares
  USE arquivos
  USE arquivos_json
  ! Metodos de integracao numerica
  USE integrador
  USE rungekutta2
  USE rungekutta3
  USE rungekutta4
  USE euler_simp
  USE verlet
  USE ruth3
  USE ruth4
  USE rkn551
  USE rkn671
  USE svcp8s15
  USE svcp10s35
  USE euler_exp
  USE euler_imp
  ! Para configuracoes
  USE json_utils_mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC simular, rodar_simulacao

  ! Tipos de metodo de integracao
  TYPE(integracao_verlet),    TARGET, SAVE :: INT_VERLET
  TYPE(integracao_euler_simp), TARGET, SAVE :: INT_EULER_SIMP
  TYPE(integracao_rk2),       TARGET, SAVE :: INT_RK2
  TYPE(integracao_rk3),       TARGET, SAVE :: INT_RK3
  TYPE(integracao_rk4),       TARGET, SAVE :: INT_RK4
  TYPE(integracao_ruth3),     TARGET, SAVE :: INT_RUTH3
  TYPE(integracao_ruth4),     TARGET, SAVE :: INT_RUTH4
  TYPE(integracao_rkn551),    TARGET, SAVE :: INT_RKN551
  TYPE(integracao_rkn671),    TARGET, SAVE :: INT_RKN671
  TYPE(integracao_svcp8s15),  TARGET, SAVE :: INT_SVCP8S15
  TYPE(integracao_svcp10s35), TARGET, SAVE :: INT_SVCP10S35
  TYPE(integracao_euler_exp), TARGET, SAVE :: INT_EULER_EXP
  TYPE(integracao_euler_imp), TARGET, SAVE :: INT_EULER_IMP

  ! Classe de simulacao
  TYPE :: simular
  
    !> N: Quantidade de corpos
    !> dim: Dimensao do problema
    INTEGER :: N, dim = 3, qntd_checkpoints, t0, tf

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
    REAL(pf) :: m_esc, m_inv, m2
    LOGICAL :: mi ! Massas iguais
    
    REAL(pf), DIMENSION(3) :: J0

    ! Metodo
    CHARACTER(LEN=:), ALLOCATABLE :: metodo

    ! Diretorio onde ficara salvo
    CHARACTER(LEN=:), ALLOCATABLE :: dir

    ! Configuracoes
    LOGICAL :: corrigir=.FALSE., colidir=.FALSE., paralelo=.FALSE.
    REAL(pf) :: corrigir_margem_erro = 0.1_pf
    INTEGER :: corrigir_max_num_tentativas = 5
    REAL(pf) :: colisoes_max_distancia = 0.1
    CHARACTER(LEN=:), ALLOCATABLE :: colisoes_modo

    !> Arquivo
    TYPE(arquivo) :: Arq

    CONTAINS
      PROCEDURE :: Iniciar, inicializar_metodo, rodar
      
  END TYPE

  ! Instanciamento da classe
  TYPE(simular) :: Simulador

CONTAINS


SUBROUTINE rodar_simulacao (infos, massas, posicoes, momentos)
  TYPE(json_value), POINTER :: infos
  REAL(pf), INTENT(IN) :: massas(:), posicoes(:,:), momentos(:,:)
  INTEGER  :: t0, tf
  REAL(pf) :: timestep

  CALL json % get(infos, 'integracao.t0', t0)
  CALL json % get(infos, 'integracao.tf', tf)
  timestep = json_get_float(infos, 'integracao.timestep')

  IF (t0 < 0) THEN
    IF (tf == 0) THEN
      ! Roda apenas o passado
      WRITE (*,*) " > Intervalo [", t0, ",", tf, "]"
      CALL rodar_simulacao_intervalo(infos, massas, posicoes, momentos, -timestep, t0, 0)

    ELSE IF (tf > 0) THEN     
      ! Roda o passado e o futuro
      WRITE (*,*) " > Intervalo [", t0, ",", 0, "]"
      CALL rodar_simulacao_intervalo(infos, massas, posicoes, momentos, -timestep, t0, 0)
      
      WRITE (*,*) " > Intervalo [", 0, ",", tf, "]"
      CALL rodar_simulacao_intervalo(infos, massas, posicoes, momentos, timestep, 0, tf)

    ENDIF
  
  ! Se for positivo, apenas roda normal
  ELSE
    ! Roda apenas o futuro
    WRITE (*,*) " > Intervalo [", 0, ",", tf, "]"
    CALL rodar_simulacao_intervalo(infos, massas, posicoes, momentos, timestep, 0, tf)
  ENDIF
END SUBROUTINE rodar_simulacao

SUBROUTINE rodar_simulacao_intervalo (infos, massas, posicoes, momentos, timestep, t0, tf)
  TYPE(json_value), POINTER :: infos
  REAL(pf), INTENT(IN) :: massas(:), posicoes(:,:), momentos(:,:)
  REAL(pf), INTENT(IN) :: timestep
  INTEGER, INTENT(IN)  :: t0, tf
  REAL(pf) :: timer_0, timer_1
  INTEGER :: qntd_total_passos
  CHARACTER(LEN=:), ALLOCATABLE :: metodo

  qntd_total_passos = (tf - t0)/timestep

  ! Timer
  timer_0 = omp_get_wtime()
  
  ! Inicializando a classe de simulacao
  CALL Simulador%Iniciar(infos, massas, posicoes, momentos, timestep)

  ! Agora roda
  CALL Simulador%rodar(tf - t0)

  metodo = json_get_string(infos, 'integracao.metodo')

  timer_1 = omp_get_wtime()
  WRITE (*,*) ' * tempo ', metodo, ': ', timer_1 - timer_0
  WRITE (*,*) ' * tempo por passo: ', (timer_1 - timer_0) / qntd_total_passos
END SUBROUTINE rodar_simulacao_intervalo



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
SUBROUTINE Iniciar (self, infos, M, R0, P0, h)

  CLASS(simular), INTENT(INOUT) :: self
  TYPE(json_value), POINTER :: infos
  REAL(pf) :: M(:), R0(:,:), P0(:,:)
  REAL(pf) :: h, colmd, densidade
  INTEGER  :: a, i
  REAL(pf) :: PI = 4.D0*DATAN(1.D0)
  LOGICAL :: encontrado

  CALL json % get(infos, 'integracao.passos_por_rodada', self % qntd_checkpoints)
  
  ! Salva o tamanho dos passos
  self % h = json_get_float(infos, 'integracao.timestep')

  ! Salva o softening do potencial
  self % potsoft = json_get_float(infos, 'integracao.amortecedor')

  ! Salva a gravidade
  self % G = json_get_float(infos, 'G')

  ! Salva as massas
  self % M = M
  
  ! Quantidade corpos no sistema
  self % N = SIZE(M)

  ! Massa total do sistema
  self % mtot = SUM(self % M)

  ! Salva se as massas sao iguais
  CALL json % get(infos, 'massas_iguais', self % mi, encontrado)
  IF (.NOT. encontrado) self % mi = .FALSE.
  IF (self % mi) THEN
    self % m_esc = self % M(1)
    self % m2 = self % m_esc * self % m_esc
    self % m_inv = 1 / self % m_esc
  ENDIF

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
  self % metodo = json_get_string(infos, 'integracao.metodo')
  CALL json % get(infos, 'integracao.t0', self % t0)
  CALL json % get(infos, 'integracao.tf', self % tf)
  
  ! Salva o corretor
  CALL json % get(infos, 'correcao.corrigir', self % corrigir)
  self % corrigir_margem_erro = json_get_float(infos, 'correcao.margem_erro')
  CALL json % get(infos, 'correcao.max_num_tentativas', self % corrigir_max_num_tentativas)

  ! Salva as colisoes
  CALL json % get(infos, 'colisoes.colidir', self % colidir)
  self % colisoes_modo = json_get_string(infos, 'colisoes.metodo')
  densidade = json_get_float(infos, 'colisoes.densidade')
  colmd = (0.75_pf / (PI * densidade))**(1.0_pf/3.0_pf)
  self % colisoes_max_distancia = colmd
  WRITE(*,*) 'COLMD: ', colmd

  ! Salva o uso de paralelizacao
  CALL json % get(infos, 'paralelo', self % paralelo)

  ! Inicializa o metodo
  CALL self % inicializar_metodo()

  ! Copia o arquivo de valores iniciais
  CALL diretorio_data()
  CALL salvar_vi_json('out/data/'//self % Arq % dirarq//'/vi', infos, M, R0, P0, .FALSE.)
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
    self%t0, self%tf, self % qntd_checkpoints, &
    self%corrigir, self%corrigir_margem_erro, self%corrigir_max_num_tentativas, &
    self%colisoes_modo, self%colisoes_max_distancia, self % paralelo)

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
  REAL(pf64) :: t0, tf, tempo_total = 0.0_pf64
  REAL(pf) :: E, Icm, Ep, Ec, ati, virial = 0.0_pf, difE, av
  REAL(pf) :: Ec_media = 0.0_pf, V_media = 0.0_pf, inst_t = 0.0_pf
  ! INTEGER  :: timestep_inv, inst = 0
  INTEGER  :: inst = 0
  INTEGER  :: qntd_por_rodada
  REAL(pf) :: timestep_inv

  ! Definicao dinamica do metodo
  WRITE (*,'(a)') 'METODO: ', self % metodo
  
  SELECT CASE (TRIM(self % metodo))
    CASE ('verlet');     integrador => INT_VERLET
    CASE ('rk2');        integrador => INT_RK2
    CASE ('rk3');        integrador => INT_RK3
    CASE ('rk4');        integrador => INT_RK4
    CASE ('euler_simp'); integrador => INT_EULER_SIMP
    CASE ('ruth3');      integrador => INT_RUTH3
    CASE ('ruth4');      integrador => INT_RUTH4
    CASE ('rkn551');     integrador => INT_RKN551
    CASE ('rkn671');     integrador => INT_RKN671
    CASE ('svcp8s15');   integrador => INT_SVCP8S15
    CASE ('svcp10s35');  integrador => INT_SVCP10S35
    CASE ('euler_exp');  integrador => INT_EULER_EXP
    CASE ('euler_imp');  integrador => INT_EULER_IMP
    
    ! Por padrao, sera o VERLET
    CASE DEFAULT;       WRITE(*,*) 'Metodo nao identificado!'
  END SELECT

  ! inicializa o integrador 
  CALL integrador % Iniciar(self % M, self % G, self % h, self % potsoft, &
    self % E0, self % J0, &
    self%corrigir, self%corrigir_margem_erro, self%corrigir_max_num_tentativas, &
    self%colidir, self%colisoes_modo, self%colisoes_max_distancia, self%paralelo, &
    self%mi)

  ! Condicoes iniciais
  R1 = self % R
  P1 = self % P
  CALL self % Arq % escrever((/R1, P1/))
  timestep_inv = 1/self % h

  ! Roda
  WRITE (*, '(a)') '  > iniciando simulacao...'
  WRITE(*,*) self % qntd_checkpoints

  qntd_por_rodada  = NINT(ABS(qntdPassos * timestep_inv) / self % qntd_checkpoints)

  DO WHILE (i < self % qntd_checkpoints)
    ! Instante
    inst_t = inst_t + qntd_por_rodada * self % h

    ! timer
    t0 = omp_get_wtime()
    ! Integracao
    CALL integrador % aplicarNVezes(R1, P1, qntd_por_rodada)
    
    ! timer
    tf = omp_get_wtime()
    tempo_total = tempo_total + tf - t0

    CALL self % Arq % escrever((/R1, P1/))

    ! Energia
    IF (self % mi) THEN
      Ep = energia_potencial_esc(self % G, self % m, R1)
      Ec = energia_cinetica_esc(self % m, P1)
    ELSE
      Ep = energia_potencial_vec(self % G, self % M, R1)
      Ec = energia_cinetica_vec(self % M, P1)
    ENDIF
    difE = Ec + Ep - self%E0

    ! virial
    Ec_media = (inst * Ec_media + Ec) / (inst + 1)
    V_media = (inst * V_media + Ep) / (inst + 1)
    virial = (Ec_media + Ec_media)/(V_media) + 1

    ! anisotropia
    ! ati = anisotropia_tensor_inercia(self % M, R1)
    ! av = anisotropia_velocidades(self % M, R1, P1)

    ! tamanho do sistema
    Icm = momento_inercia(self % M, R1)

    WRITE (*,*) '     -> Passo:', i, ' / Tempo: ', tempo_total, ' / t = ', inst_t
    WRITE (*,*) '        E-E0:', difE, ' / Virial:', virial, ' / I: ', Icm
    WRITE (*,*)
    inst = inst + 1
    CALL self % Arq % arquivo_bkp(i*self%qntd_checkpoints*NINT(timestep_inv), tempo_total)

    i = i + 1
  END DO

  CALL self % Arq % atualizar_arquivo_info(qntdPassos*self%qntd_checkpoints*NINT(timestep_inv), tempo_total)
  CALL self % Arq % excluir_bkp()

  WRITE (*, '(a)') '  > simulacao encerrada!'
  WRITE (*,*)

  CALL self % Arq % fechar()

END SUBROUTINE rodar

END module simulacao