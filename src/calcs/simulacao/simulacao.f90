! ************************************************************
!! SIMULACAO: GERAL
!
! Objetivos:
!   Arquivo base para fazer simulacoes.
!
! Modificado:
!   11 de novembro de 2025
!
! Autoria:
!   oap
! 
MODULE simulacao

  USE OMP_LIB
!> Versionamento
  USE version
!> Tipos
  USE tipos
!> Json-Fortran e arquivos
  USE arquivos_mod
!> Biblioteca de mecanica
  USE utilidades
!> Integracao numerica
  USE integrador
  USE integradores
!> Para plot em tempo real
  USE conexao
!> Colisoes
  USE colisao
!> Correcao numerica
  USE correcao

  IMPLICIT NONE
  PRIVATE
  PUBLIC simular

!> Classe de simulacao
  TYPE :: simular
    
    !> Ponteiro para o dicionario de valores iniciais
    TYPE(json_value), POINTER :: infos => NULL()

    !> N: Quantidade de corpos
    !> dim: Dimensao do problema
    INTEGER :: N, dim = 3, qntd_checkpoints, t0, tf

    !> h: Tamanho do passo de integracao
    !> G: Constante de gravitacao universal
    !> mtot: Massa total do sistema
    !> potsoft: Softening do potencial
    !> virial: 2*ec + <F,q>
    REAL(pf) :: h, G, mtot, potsoft, virial
    
    !> M: Massas do sistema
    !> R0: Posicoes das particulas (inicial)
    !> P0: Momento linear das particulas (inicial)
    REAL(pf), ALLOCATABLE :: M(:), R0(:,:), P0(:,:)
    REAL(pf) :: m_esc, m_inv, m2
    LOGICAL :: mi ! Massas iguais

    !> Integrais primeiras
    ! E0: Energia total inicial
    ! J0: Momento angular total inicial
    ! Ptot0: Momento linear total inicial
    ! Rcm0: Centro de massas inicial
    REAL(pf), ALLOCATABLE :: E0, J0(:), Ptot0(:), Rcm0(:)
    
    ! Metodo
    CHARACTER(LEN=:), ALLOCATABLE :: metodo
    CLASS(integracao), POINTER :: integrador

    ! Diretorio onde ficara salvo
    CHARACTER(LEN=:), ALLOCATABLE :: dir

    ! Correcao
    LOGICAL :: corrigir=.FALSE., paralelo=.FALSE., gpu=.FALSE.
    REAL(pf) :: corme = 0.1_pf
    INTEGER :: cormnt = 5
    ! Colisoes
    LOGICAL :: colidir = .FALSE.
    REAL(pf) :: colmd
    CHARACTER(LEN=:), ALLOCATABLE :: colisoes_modo
    REAL(pf), ALLOCATABLE :: raios(:)

    ! Exibir
    LOGICAL :: exibir = .FALSE.
    TYPE(conexao_socket) :: conexao

    !> Arquivo
    TYPE(arquivo) :: Arq

    CONTAINS
      PROCEDURE :: iniciar, &  ! Inicializacao principal
                   inicializar_data, &  ! Inicializacao do arquivo "data"
                   inicializar_metodo, & ! Inicializacao do integrador e do socket
                   rodar, &   
                   output_passo

  END TYPE simular

  ! Instanciamento
  TYPE(simular) :: simulador

CONTAINS

SUBROUTINE iniciar (self, infos, m, R0, P0, h, out_dir, out_ext)
  CLASS(simular), INTENT(INOUT) :: self
  CHARACTER(LEN=*), INTENT(IN) :: out_dir, out_ext
  TYPE(json_value), POINTER :: infos
  REAL(pf) :: m(:), R0(:,:), P0(:,:)
  REAL(pf) :: h, colmd, densidade, ec, f_prod_q
  REAL(pf) :: PI = 4.D0*DATAN(1.D0)
  LOGICAL :: encontrado
  INTEGER :: a

  !> Faz uma copia do dicionario de informacoes
  CALL json_clone(infos, self % infos)

  !## Variaveis e constantes do sistema ##!
  !> Constante de gravitacao universal
  self % G = json_get_float(infos, "G")
  !> Salva o softening do potencial
  self % potsoft = json_get_float(infos, 'integracao.amortecedor')
  !> Salva a dimensao
  self % dim = SIZE(R0,2)
  !> Variaveis de estado
  self % m = m  ! massas
  self % R0 = R0 ! posicoes
  self % P0 = P0 ! momentos
  !> Salva se as massas sao iguais
  CALL json % get(infos, 'massas_iguais', self % mi, encontrado)
  IF (.NOT. encontrado) self % mi = .FALSE.
  IF (self % mi) THEN
    self % m_esc = self % M(1)
    self % m2 = self % m_esc * self % m_esc
    self % m_inv = 1 / self % m_esc
  ENDIF
  !> Quantidade de corpos
  self % N = SIZE(m)
  !> Integrais primeiras
  self % E0 = energia_total(self%G, self%m, R0, P0, self%potsoft)
  self % J0 = momento_angular_total(R0, P0)
  self % Ptot0 = momento_linear_total(P0)
  self % rcm0 = centro_massas(self%m, R0)
  !> Equilibrio
  ec = energia_cinetica(self%m, P0)
  IF (self%potsoft == 0) THEN
    self % virial = ec + self % E0
  ELSE
    f_prod_q = virial_potencial_amortecido(self%G, self%m, R0, self%potsoft)
    self % virial = ec + ec + f_prod_q
  ENDIF

  !## Uso da paralelizacao ##!
  CALL json % get(infos, 'paralelo', self%paralelo)
  CALL json % get(infos, 'gpu', self%gpu, encontrado)
  IF (.NOT. encontrado) self%gpu = .FALSE.

  !## Sobre a integracao numerica ##!
  CALL self % inicializar_metodo(h)

  !## Sobre a correcao numerica ##!
  CALL json % get(infos, 'correcao.corrigir', self % corrigir)
  self % corme = json_get_float(infos, 'correcao.margem_erro')
  CALL json % get(infos, 'correcao.max_num_tentativas', self % cormnt)

  !## Sobre as colisoes ##!
  CALL json % get(infos, 'colisoes.colidir', self % colidir)
  self % colisoes_modo = json_get_string(infos, 'colisoes.metodo')
  densidade = json_get_float(infos, 'colisoes.densidade')
  colmd = (0.75_pf / (PI * densidade))**(1.0_pf/3.0_pf)
  self % colmd = colmd
  ! Deixando os raios calculados de antemao
  ALLOCATE(self % raios(self % N))
  DO a = 1, self % N
    self % raios(a) = self % colmd * m(a)**(1.0_pf / 3.0_pf)
  END DO

  !> Copia o arquivo de valores iniciais
  CALL diretorio_data(out_dir)
  !> Inicializa o arquivo data.csv
  CALL self % inicializar_data(out_dir, out_ext)
  !> Salva os valores iniciais
  CALL salvar_vi_json(out_dir//'/data/'//self % Arq % dir_arq//'/vi', infos, M, R0, P0, .FALSE.)

END SUBROUTINE iniciar

! ************************************************************
!! Inicializar o arquivo "data"
!
! Modificado:
!   11 de novembro de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE inicializar_data (self, out_dir, out_ext)
  CLASS(simular), INTENT(INOUT) :: self
  CHARACTER(LEN=*), INTENT(IN) :: out_dir
  CHARACTER(LEN=*), INTENT(IN) :: out_ext

  ! Define o diretorio de saida
  CALL self % Arq % definir_diretorio_saida(out_dir, out_ext)

  ! Cria o arquivo onde sera salvo
  CALL self % Arq % criar_data(self % N, self % dim)
  self % dir = self % Arq % dir_arq

  ! Salva as infos de cabecalho
  CALL self % Arq % escrever_cabecalho_data(self % h, self % G, self % M)

  ! Salva as informacoes no info.txt
  CALL self % Arq % inicializar_arquivo_info(self % infos, version_string, precisao)

  ! Condicoes iniciais
  CALL self % Arq % escrever_data((/self % R0, self % P0/))
END SUBROUTINE inicializar_data

! ************************************************************
!! Inicializa o integrador e o socket
!
! Modificado:
!   08 de agosto de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE inicializar_metodo (self, h)
  CLASS(simular), INTENT(INOUT) :: self
  REAL(pf), INTENT(IN) :: h
  LOGICAL :: encontrado

  self % h = h  ! Timestep
  !> Metodo
  self % metodo = json_get_string(self%infos, 'integracao.metodo')
  CALL json % get(self%infos, 'integracao.t0', self % t0)
  CALL json % get(self%infos, 'integracao.tf', self % tf)
  !> Quantidade de checkpoints
  CALL json % get(self%infos, 'integracao.checkpoints', self % qntd_checkpoints)
  !> Se quer ou nao plotar durante a simulacao 
  CALL json % get(self%infos, 'exibir', self % exibir, encontrado)
  IF (.NOT. encontrado) self % exibir = .FALSE.
  !> Definindo o metodo
  CALL definir_metodo(self%integrador, self%metodo)

  ! Inicializando o integrador
  CALL self % integrador % iniciar(self%infos, self%h, self%M, self%E0, self%J0)

  ! Atualizando as constantes se for necessario
  CALL self % integrador % atualizar_constantes()

  ! Se for plotar em tempo real, precisa inicializar tambem
  IF (self % exibir) CALL self % conexao % inicializar_plot_tempo_real(self % N)
END SUBROUTINE


SUBROUTINE rodar (self, qntdPassos)
  CLASS(simular), INTENT(INOUT) :: self
  INTEGER, INTENT(IN) :: qntdPassos
  INTEGER :: qntd_por_rodada, passo, sub_passo
  REAL(pf) :: inst_t, subinst_t, t0, tf
  REAL(pf64)  :: tempo_total
  REAL(pf) :: R1(self%N,self%dim), P1(self%N,self%dim)
  REAL(pf) :: E
  LOGICAL :: corrigiu

  WRITE (*, '(a)') '  > iniciando simulacao...'
  R1 = self % R0
  P1 = self % P0

  qntd_por_rodada = NINT(ABS(qntdPassos / self%h) / self%qntd_checkpoints)

  passo = 0
  inst_t = 0
  subinst_t = 0
  tempo_total = 0

  DO WHILE (passo < self % qntd_checkpoints)
    !> Atualizando o instante
    inst_t = inst_t + qntd_por_rodada * self%h
    !> Timer
    t0 = OMP_GET_WTIME()

    !> Simulacao em si
    DO sub_passo = 1, qntd_por_rodada
      !> Integracao
      subinst_t = subinst_t + self % h
      CALL self % integrador % aplicar(R1, P1)

      !> Colide, se for o caso
      IF (self % colidir) THEN
        CALL verificar_e_colidir(self%m, R1, P1, self%paralelo, &
                                self%raios, self%colisoes_modo, &
                                self%integrador%distancias)
      END IF

      !> Se for exibir, envia os dados
      IF (self % exibir) THEN
        !> Envia nos subpassos pares
        IF (MOD(sub_passo,2) == 0) CALL self%conexao%enviar(subinst_t, R1, P1)
      END IF
    END DO

    !> Se for corrigir, calcula a energia total para ver se precisa
    IF (self % corrigir) THEN
      E = energia_total(self%G, self%m, R1, P1, self%potsoft, self%integrador%distancias)

      !> Se o erro for maior que o permitido, corrige
      IF (ABS(E - self%E0) >= self%corme) THEN
        !> (Corrige energia total e momento angular total)[CUSTOSO]
        ! CALL corrigir(self%corme, self%cormnt, self%G, self%m, R1, P1, corrigiu, self%E0, self%J0)

        !> (Corrige somente a energia total)
        CALL corrigir_apenas_energia(self % corme, self % cormnt, self % G, &
                                  self % m, R1, P1, corrigiu, self % E0, self % J0, E, &
                                  self % potsoft)
      ENDIF
    END IF

    !> Timer
    tf = OMP_GET_WTIME()
    tempo_total = tempo_total + (tf - t0)

    !> Salvando o output
    CALL self % arq % escrever_data((/R1, P1/))
    CALL self % arq % atualizar_arquivo_bkp(passo * self%qntd_checkpoints*NINT(1.0_pf/self%h), tempo_total)

    !> Imprimindo informacoes
    CALL self % output_passo(passo+1, tempo_total, inst_t, R1, P1)

    !> Atualiza o passo
    passo = passo + 1
  END DO

  CALL self % arq % atualizar_arquivo_info(qntdPassos*self%qntd_checkpoints*NINT(1.0_pf/self%h), tempo_total)
  CALL self % arq % excluir_bkp()

  IF (self % exibir) CALL self % conexao % encerrar_conexao()

  WRITE (*, '(a)') '  > simulacao encerrada!'
  WRITE (*,*)

  CALL self % arq % fechar()

END SUBROUTINE

! ************************************************************
!! Mensagens na tela que aparecem durante a simulacao
!
! Modificado:
!   25 de julho de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE output_passo (self, i, tempo_total, inst_t, R, P)
  CLASS(simular), INTENT(INOUT) :: self
  INTEGER :: i
  REAL(pf64) :: tempo_total
  REAL(pf) :: inst_t
  REAL(pf), DIMENSION(self % N, 3) :: R, P

  REAL(pf) :: EP, EC, errene
  REAL(pf) :: virial, f_prod_q
  REAL(pf) :: momine, momdil
  CHARACTER(100) :: char_real
  CHARACTER(300) :: saida

  ! Se tiver amortecimento, calcula o potencial de um jeito diferente
  IF (self % potsoft .NE. 0) THEN
    f_prod_q = virial_potencial_amortecido(self%G, self%M, R, self%potsoft, Ep)
    Ec = energia_cinetica_vec(self % M, P)
    virial = Ec + Ec + f_prod_q
  ELSE
    ! Energia
    IF (self % mi) THEN
      Ep = energia_potencial_esc(self % G, self % m_esc, R, self % potsoft)
      Ec = energia_cinetica_esc(self % m_esc, P)
    ELSE
      Ep = energia_potencial_vec(self % G, self % M, R, self % potsoft)
      Ec = energia_cinetica_vec(self % M, P)
    ENDIF
    virial = Ec + Ec + Ep
  ENDIF

  ! Erro na energia
  errene = EC + EP - self % E0

  ! Tamanho do sistema
  momine = momento_inercia(self % M, R)
  momdil = momento_dilatacao(R, P)

  ! Virial
  self % virial = (i * self % virial + virial)/(i+1)
  
  ! Montando a string de output
  WRITE(saida, '(5X,A,I8,4X,A,F12.3,2X,A,F10.2)') '-> Passo:', i, ' / Tempo:', tempo_total, ' / t:', inst_t

  WRITE(char_real, '(9X,A,E12.4,A,E12.4)') 'E-E0: ', errene, ' / Virial: ', self % virial
  saida = TRIM(saida)//char(10)//TRIM(char_real)

  WRITE(char_real, '(A,E12.4,A,E12.4)') 'I: ', momine, ' / D: ', momdil
  saida = TRIM(saida)//" / "//TRIM(char_real)
  
  WRITE (*,*) TRIM(saida)
  WRITE (*,*)

END SUBROUTINE

END MODULE simulacao
