! ************************************************************
!! INTEGRADOR NUMERICO
!
! Objetivos:
!   Estrutura basica de um integrador numerico. As configuracoes
!   gerais de um metodo devem vir aqui, como tamanho de passo,
!   dimensao, massas, se corrige ou nao, se colide ou nao, etc.
!
! Modificado:
!   24 de julho de 2025
!
! Autoria:
!   oap
! 
MODULE integrador

  USE tipos
  USE OMP_LIB
  USE funcoes_forca
  USE funcoes_forca_mi
  USE mecanica
  USE correcao
  USE colisao
  USE octree
  USE json_utils_mod
  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao

  TYPE :: integracao
  
    ! m: Massas
    REAL(pf), ALLOCATABLE :: m(:)
    REAL(pf), ALLOCATABLE :: massasInvertidas(:,:)
    
    ! Massas iguais
    LOGICAL  :: mi
    REAL(pf) :: m_esc, m_inv, m2
    REAL(pf128) :: m_esc_128, m_inv_128, m2_128

    ! Distancias entre os corpos
    REAL(pf), ALLOCATABLE :: distancias(:)

    ! h: Passo de integracao
    ! G: Constante de gravitacao
    ! potsoft: Softening do potencial
    ! E0: Energia total inicial
    ! J0: Momento angular total inicial
    REAL(pf) :: h, G, potsoft, potsoft2, E0, J0(3)

    ! dim: Dimensao do problema
    ! N: Quantidade de partículas
    INTEGER :: dim = 3, N

    ! Se vai ou nao corrigir
    LOGICAL  :: corrigir = .FALSE.
    REAL(pf) :: corme ! margem de erro
    INTEGER  :: cormnt ! max num tentativas

    ! Se vai ou nao colidir
    LOGICAL       :: colidir = .FALSE.
    CHARACTER(:), ALLOCATABLE :: colisoes_modo
    REAL(pf)      :: colmd ! max dist colisoes
    REAL(pf), ALLOCATABLE :: raios(:)

    ! Se vai ou nao usar paralelizacao
    LOGICAL :: paralelo = .FALSE., gpu = .FALSE.

    ! Funcao de forcas (aceleracao)
    PROCEDURE(forcas_funcbase), POINTER, NOPASS    :: forcas_funcao    => NULL()
    PROCEDURE(forcas_mi_funcbase), POINTER, NOPASS :: forcas_mi_funcao => NULL()

    ! Arvore octree
    TYPE(arvore_octo), ALLOCATABLE :: arvore

    CONTAINS
      PROCEDURE :: Iniciar, aplicarNVezes, metodo, metodo_mi, forcas, atualizar_constantes
  
  END TYPE integracao

  ABSTRACT INTERFACE
    FUNCTION forcas_funcbase (m, R, G, N, dim, potsoft2, distancias)
        IMPORT :: pf
        IMPLICIT NONE
        INTEGER,                     INTENT(IN) :: N, dim
        REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
        REAL(pf), DIMENSION(N),      INTENT(IN) :: m
        REAL(pf),                    INTENT(IN) :: G, potsoft2
        REAL(pf), DIMENSION(INT(N*(N-1)/2)), INTENT(INOUT) :: distancias
        REAL(pf), DIMENSION(N, dim) :: forcas_funcbase
    END FUNCTION forcas_funcbase

    FUNCTION forcas_mi_funcbase (R, G, N, dim, potsoft2, distancias)
        IMPORT :: pf
        IMPLICIT NONE
        INTEGER,                     INTENT(IN) :: N, dim
        REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
        REAL(pf),                    INTENT(IN) :: G, potsoft2
        REAL(pf), DIMENSION(INT(N*(N-1)/2)), INTENT(INOUT) :: distancias
        REAL(pf), DIMENSION(N, dim) :: forcas_mi_funcbase
    END FUNCTION forcas_mi_funcbase
  END INTERFACE
CONTAINS

! ************************************************************
!! Construtor da classe
!
! Objetivos:
!   Define o principal, salvando os valores e inicializando o
!   metodo.
!
! Modificado:
!   24 de julho de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE Iniciar (self, infos, timestep, massas, E0, J0)
  IMPLICIT NONE
  class(integracao), INTENT(INOUT) :: self
  TYPE(json_value), POINTER :: infos
  REAL(pf) :: timestep
  REAL(pf), allocatable :: massas(:)
  REAL(pf)              :: E0, J0(3)
  INTEGER :: a, i

  LOGICAL :: encontrado
  REAL(pf) :: colmd, densidade
  REAL(pf) :: PI = 4.D0*DATAN(1.D0)

  ! Quantidade de particulas
  self % N = SIZE(massas)
  
  ! Massas
  CALL json % get(infos, 'massas_iguais', self % mi, encontrado)
  IF (.NOT. encontrado) self % mi = .FALSE.
  ALLOCATE(self % m(self % N))
  self % m = massas

  IF (self % mi) THEN
    self % m_esc = massas(1)
    self % m_inv = 1.0_pf/self % m_esc
    self % m2 = self % m_esc * self % m_esc

    self % m_esc_128 = REAL(massas(1), KIND=pf128)
    self % m_inv_128 = 1.0_pf128 / self % m_esc_128
    self % m2_128 = self % m_esc_128 * self % m_esc_128
    
    ! Se 1/m - N < 1e-10, assume que m = 1/N
    ! Nesse caso, podemos melhorar a precisao
    IF (ABS(self % m_inv - self % N) < 1E-10) THEN
      self % m_esc = 1.0_pf / (self % N)
      self % m2 = 1.0_pf / (self % N * self % N)
      self % m_inv = self % N

      self % m_esc_128 = 1.0_pf128 / (self % N)
      self % m2_128 = 1.0_pf128 / (self % N * self % N)
      self % m_inv_128 = self % N
    ! Se m = 1
    ELSE IF (ABS(self % m_esc - 1) < 1E-10) THEN
      self % m_esc = 1.0_pf
      self % m2 = 1.0_pf
      self % m_inv = 1.0_pf

      self % m_esc_128 = 1.0_pf128
      self % m2_128 = 1.0_pf128
      self % m_inv_128 = 1.0_pf128
    ENDIF
  ELSE  
    ! vetor de massas invertidas
    ALLOCATE(self % massasInvertidas (self % N, self % dim))
    DO a = 1, self % N
      DO i = 1, self % dim
        self % massasInvertidas(a,i) = 1.0_pf/(massas(a))
      END DO
    END DO
  ENDIF  

  ! gravidade
  self % G = json_get_float(infos, 'G')
  ! Passo
  self % h = timestep
  ! Softening do potencial
  self % potsoft = json_get_float(infos, 'integracao.amortecedor')
  self % potsoft2 = self%potsoft * self%potsoft

  ! Valores iniciais
  self % E0 = E0
  self % J0 = J0

  ! Distancias
  ALLOCATE(self % distancias(INT(self%N * (self%N-1)/2)))

  ! Se vai ou nao corrigir
  CALL json % get(infos, 'correcao.corrigir', self % corrigir)
  self % corme = json_get_float(infos, 'correcao.margem_erro')
  CALL json % get(infos, 'correcao.max_num_tentativas', self % cormnt)

  ! Colisoes
  CALL json % get(infos, 'colisoes.colidir', self % colidir)
  self % colisoes_modo = json_get_string(infos, 'colisoes.metodo')
  densidade = json_get_float(infos, 'colisoes.densidade')
  colmd = (0.75_pf / (PI * densidade))**(1.0_pf/3.0_pf)
  self % colmd = colmd

  ! Deixando os raios calculados de antemao
  ALLOCATE(self % raios(self % N))
  DO a = 1, self % N
    self % raios(a) = self % colmd * massas(a)**(1.0_pf / 3.0_pf)
  END DO

  ! Codigo paralelo
  CALL json % get(infos, 'paralelo', self % paralelo)
  CALL json % get(infos, 'gpu', self % gpu, encontrado)
  IF (.NOT. encontrado) self % gpu = .FALSE.

  IF (self % GPU) THEN
#ifdef USAR_GPU
    IF (self % mi) THEN
      self % forcas_mi_funcao => forcas_mi_par_gpu
    ELSE
      self % forcas_funcao => forcas_par_gpu
    ENDIF
#else
    WRITE (*,*) "GPU nao compilada. Recompile o programa com -DUSAR_GPU=ON"
    WRITE (*,*) "ou desative a opcao de gpu."
    STOP 0
#endif
  ELSE
    IF (self % paralelo) THEN
      IF (self % mi) THEN
        self % forcas_mi_funcao => forcas_mi_par
      ELSE
        self % forcas_funcao => forcas_par
      ENDIF
    ELSE
      IF (self % mi) THEN
        self % forcas_mi_funcao => forcas_mi_seq
      ELSE
        self % forcas_funcao => forcas_seq
      ENDIF
    ENDIF
  ENDIF

END SUBROUTINE Iniciar

FUNCTION forcas (self, R)
  IMPLICIT NONE
  class(integracao), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self % N, self % dim), INTENT(IN) :: R
  REAL(pf), DIMENSION(self % N, self % dim) :: forcas
  
  IF (self % mi) THEN
    forcas = self % forcas_mi_funcao(R, self%G, self%N, self%dim, &
                    self%potsoft2, self%distancias)
  ELSE
    forcas = self % forcas_funcao(self % m, R, self%G, self%N, self%dim, &
                    self%potsoft2, self%distancias)
  ENDIF

END FUNCTION forcas

! ************************************************************
!! Aplicacao iterada do metodo
!
! Objetivos:
!   Aplica o metodo iterativamente N vezes.
!
! Modificado:
!   24 de julho de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE aplicarNVezes (self, R, P, qntd_passos)

  IMPLICIT NONE
  class (integracao), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P
  INTEGER, INTENT(IN) :: qntd_passos
  REAL(pf)             :: E
  ! Para cada passo
  INTEGER :: i
  ! para verificar se corrigiu
  LOGICAL :: corrigiu = .FALSE.
  ! Para as forcas e passos pos-integracao
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, FSomas_ant
  REAL(pf), DIMENSION(3, self%N, self%dim) :: resultado
  ! Consumo de tempo com correcao
  REAL(pf) :: t0_cor, tempo_correcao
  INTEGER :: contagem_correcao

  contagem_correcao = 0

  ! Salvando as primeiras posicoes e momentos
  R1 = R
  P1 = P

  ! Calcula as forcas
  FSomas_ant = self%forcas(R)

  ! Integrando (massas iguais)
  IF (self % mi) THEN  
    DO i = 1, qntd_passos
      ! Aplica o metodo
      resultado = self % metodo_mi(R1, P1, FSomas_ant)

      R1 = resultado(1,:,:)
      P1 = resultado(2,:,:)
      FSomas_ant = resultado(3,:,:)

      ! se tiver colisoes, aplica
      IF (self % colidir) THEN
        CALL verificar_e_colidir(self%m, R1, P1, self%paralelo, &
                                self%raios, self%arvore, self%colisoes_modo, &
                                self%distancias)
      ENDIF
    END DO
  ! Integrando (massas diferentes)
  ELSE
    DO i = 1, qntd_passos
      ! Aplica o metodo
      resultado = self % metodo(R1, P1, FSomas_ant)

      R1 = resultado(1,:,:)
      P1 = resultado(2,:,:)
      FSomas_ant = resultado(3,:,:)

      ! se tiver colisoes, aplica
      IF (self % colidir) THEN
        CALL verificar_e_colidir(self%m, R1, P1, self%paralelo, &
                                self%raios, self%arvore, self%colisoes_modo, &
                                self%distancias)
      ENDIF
    END DO
  ENDIF

  ! Se estiver disposto a corrigir, calcula a energia total para ver se precisa
  IF (self%corrigir) THEN
    E = energia_total(self % G, self % m, R1, P1, self % potsoft, self % distancias)
        
    IF (ABS(E - self%E0) >= self%corme) THEN
      
      t0_cor = omp_get_wtime()

      ! Correcao com energia total e momento angular total (desativado)
      ! CALL corrigir(self%corme,self%cormnt,self % G,self % m,R1,P1,corrigiu,E0,J0)

      ! Corrige somente a energia total
      CALL corrigir_apenas_energia(self % corme, self % cormnt, self % G, &
                                  self % m, R1, P1, corrigiu, self % E0, self % J0, E, &
                                  self % potsoft)
      
      IF (corrigiu) THEN
        contagem_correcao = contagem_correcao + 1
      END IF
      
      tempo_correcao = tempo_correcao + (omp_get_wtime() - t0_cor)
    END IF
  ENDIF

  R = R1
  P = P1

  IF (contagem_correcao > 0) THEN
    WRITE (*,*) contagem_correcao, ' correcoes, tempo: ', tempo_correcao
  ENDIF

END SUBROUTINE aplicarNVezes

! ************************************************************
!! Metodo numerico
!
! Objetivos:
!   Aplicacao do metodo em si.
!
! Modificado:
!   12 de julho de 2025
!
! Autoria:
!   oap
!
FUNCTION metodo (self, R, P, FSomas_ant)
  IMPLICIT NONE
  class(integracao), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo
  
  ! Cada integrador precisa ter um metodo definido, que substitui esta funcao vazia
  WRITE (*,*) 'OPS'

END FUNCTION metodo

! ************************************************************
!! Metodo numerico (massas iguais)
!
! Objetivos:
!   Aplicacao do metodo em si.
!
! Modificado:
!   12 de julho de 2025
!
! Autoria:
!   oap
!
FUNCTION metodo_mi (self, R, P, FSomas_ant)
  IMPLICIT NONE
  class(integracao), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo_mi
  
  ! Cada integrador precisa ter um metodo definido, que substitui esta funcao vazia
  WRITE (*,*) 'OPS'

END FUNCTION metodo_mi

! ************************************************************
!! Atualiza as constantes se necessario
!
! Modificado:
!   03 de junho de 2025
!
! Autoria:
!   oap
!
SUBROUTINE atualizar_constantes (self)
  IMPLICIT NONE
  class(integracao), INTENT(IN) :: self
END SUBROUTINE atualizar_constantes

END MODULE integrador
