! ************************************************************
!! INTEGRADOR NUMERICO
!
! Objetivos:
!   Estrutura basica de um integrador numerico. As configuracoes
!   gerais de um metodo devem vir aqui, como tamanho de passo,
!   dimensao, massas, se corrige ou nao, se colide ou nao, etc.
!
! Modificado:
!   03 de junho de 2025
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
    LOGICAL :: paralelo = .FALSE.

    ! Funcao de forcas (aceleracao)
    PROCEDURE(forcas_funcbase), POINTER, NOPASS    :: forcas_funcao    => NULL()
    PROCEDURE(forcas_mi_funcbase), POINTER, NOPASS :: forcas_mi_funcao => NULL()

    ! Arvore octree
    TYPE(arvore_octo), ALLOCATABLE :: arvore

    CONTAINS
      PROCEDURE :: Iniciar, aplicarNVezes, metodo, metodo_mi, forcas, atualizar_constantes
  
  END TYPE integracao

  ABSTRACT INTERFACE
    FUNCTION forcas_funcbase (m, R, G, N, dim, potsoft, potsoft2)
        IMPORT :: pf
        IMPLICIT NONE
        INTEGER,                     INTENT(IN) :: N, dim
        REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
        REAL(pf), DIMENSION(N),      INTENT(IN) :: m
        REAL(pf),                    INTENT(IN) :: G, potsoft, potsoft2
        REAL(pf), DIMENSION(N, dim) :: forcas_funcbase
    END FUNCTION forcas_funcbase

    FUNCTION forcas_mi_funcbase (R, G, N, dim, potsoft, potsoft2)
        IMPORT :: pf
        IMPLICIT NONE
        INTEGER,                     INTENT(IN) :: N, dim
        REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
        REAL(pf),                    INTENT(IN) :: G, potsoft, potsoft2
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
!   03 de junho de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE Iniciar (self, massas, G, h, potsoft, E0, J0, corrigir, corme, cormnt, colidir, colmodo, colmd, paralelo, mi)
  IMPLICIT NONE
  class(integracao), INTENT(INOUT) :: self
  REAL(pf), allocatable :: massas(:)
  REAL(pf)              :: G, h, E0, J0(3)
  LOGICAL,INTENT(IN) :: corrigir, paralelo
  LOGICAL, INTENT(IN) :: colidir
  CHARACTER(LEN=*), INTENT(IN) :: colmodo
  REAL(pf) :: corme, potsoft, colmd
  INTEGER :: cormnt
  INTEGER :: a, i
  LOGICAL :: mi

  ! Quantidade de particulas
  self % N = SIZE(massas)
  
  ! Massas
  self % mi = mi
  ALLOCATE(self % m(self % N))
  self % m = massas

  IF (mi) THEN
    self % m_esc = massas(1)
    self % m_inv = 1.0_pf/self % m_esc
    self % m2 = self % m_esc * self % m_esc

    self % m_esc_128 = REAL(massas(1), KIND=pf128)
    self % m_inv_128 = 1.0_pf128 / self % m_esc_128
    self % m2_128 = self % m_esc_128 * self % m_esc_128
    
    ! Se 1/m - N < 1e-10, assume que m = 1/N
    ! Nesse caso, podemos melhorar a precisao
    IF (self % m_inv - self % N < 1E-10) THEN
        self % m_esc = 1.0_pf / (self % N)
        self % m2 = 1.0_pf / (self % N * self % N)
        self % m_inv = self % N

        self % m_esc_128 = 1.0_pf128 / (self % N)
        self % m2_128 = 1.0_pf128 / (self % N * self % N)
        self % m_inv_128 = self % N
    ! Se m = 1
    ELSE IF (self % m_esc - 1 < 1E-10) THEN
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
        self % massasInvertidas(a,i) = 1/(massas(a))
      END DO
    END DO
  ENDIF  

  ! gravidade
  self % G = G
  ! Passo
  self % h = h
  ! Softening do potencial
  self % potsoft = potsoft
  self % potsoft2 = potsoft*potsoft

  ! Valores iniciais
  self % E0 = E0
  self % J0 = J0

  ! Se vai ou nao corrigir
  self % corrigir = corrigir
  self % corme = corme
  self % cormnt = cormnt

  ! Se vai ou nao colidir
  self % colidir = colidir
  self % colisoes_modo = colmodo
  self % colmd = colmd
  ALLOCATE(self % raios(self % N))
  DO a = 1, self % N
    self % raios(a) = self % colmd * massas(a)**(1.0_pf / 3.0_pf)
  END DO

  ! Codigo paralelo
  self % paralelo = paralelo
  IF (paralelo) THEN
    IF (mi) THEN
      self % forcas_mi_funcao => forcas_mi_par
    ELSE
      self % forcas_funcao => forcas_par
    ENDIF
  ELSE
    IF (mi) THEN
      self % forcas_mi_funcao => forcas_mi_seq
    ELSE
      self % forcas_funcao => forcas_seq
    ENDIF
  ENDIF
END SUBROUTINE Iniciar

FUNCTION forcas (self, R)
  IMPLICIT NONE
  class(integracao), INTENT(IN) :: self
  REAL(pf), DIMENSION(self % N, self % dim), INTENT(IN) :: R
  REAL(pf), DIMENSION(self % N, self % dim) :: forcas
  
  IF (self % mi) THEN
    forcas = self % forcas_mi_funcao(R, self%G, self%N, self%dim, self%potsoft, self%potsoft2)
  ELSE
    forcas = self % forcas_funcao(self % m, R, self%G, self%N, self%dim, self%potsoft, self%potsoft2)
  ENDIF

END FUNCTION forcas

! ************************************************************
!! Aplicacao iterada do metodo
!
! Objetivos:
!   Aplica o metodo iterativamente N vezes.
!
! Modificado:
!   03 de junho de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE aplicarNVezes (self, R, P, qntd_checkpoints)

  IMPLICIT NONE
  class (integracao), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P
  INTEGER, INTENT(IN) :: qntd_checkpoints
  REAL(pf)             :: E
  REAL(pf), DIMENSION(3) :: J
  ! Para cada passo
  INTEGER :: i, a, b
  ! para verificar se corrigiu
  LOGICAL :: corrigiu = .FALSE.
  ! Para as forcas e passos pos-integracao
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, FSomas_ant
  REAL(pf), DIMENSION(3, self%N, self%dim) :: resultado
  ! Energia total aproximada
  REAL(pf) :: Et_aprox
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
    DO i = 1, qntd_checkpoints
      ! Aplica o metodo
      resultado = self % metodo_mi(R1, P1, FSomas_ant)

      R1 = resultado(1,:,:)
      P1 = resultado(2,:,:)
      FSomas_ant = resultado(3,:,:)

      ! se tiver colisoes, aplica
      IF (self % colidir) THEN
        CALL verificar_e_colidir(self%m, R1, P1, self%colmd, self%paralelo, &
                                self%raios, self%arvore, self%colisoes_modo)
      ENDIF
    END DO
  ! Integrando (massas diferentes)
  ELSE
    DO i = 1, qntd_checkpoints
      ! Aplica o metodo
      resultado = self % metodo(R1, P1, FSomas_ant)

      R1 = resultado(1,:,:)
      P1 = resultado(2,:,:)
      FSomas_ant = resultado(3,:,:)

      ! se tiver colisoes, aplica
      IF (self % colidir) THEN
        CALL verificar_e_colidir(self%m, R1, P1, self%colmd, self%paralelo, &
                                self%raios, self%arvore, self%colisoes_modo)
      ENDIF
    END DO
  ENDIF

  ! Se estiver disposto a corrigir, calcula a energia total para ver se precisa
  IF (self%corrigir) THEN
    E = energia_total(self % G, self % m, R1, P1)
        
    IF (ABS(E - self%E0) >= self%corme) THEN
      
      t0_cor = omp_get_wtime()

      ! Correcao com energia total e momento angular total (desativado)
      ! CALL corrigir(self%corme,self%cormnt,self % G,self % m,R1,P1,corrigiu,E0,J0)

      ! Corrige somente a energia total
      CALL corrigir_apenas_energia(self % corme, self % cormnt, self % G, &
                                  self % m, R1, P1, corrigiu, self % E0, self % J0, E)
      
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
!   14 de setembro de 2024
!
! Autoria:
!   oap
!
FUNCTION metodo (self, R, P, FSomas_ant)
  IMPLICIT NONE
  class(integracao), INTENT(IN) :: self
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
!   01 de junho de 2025
!
! Autoria:
!   oap
!
FUNCTION metodo_mi (self, R, P, FSomas_ant)
  IMPLICIT NONE
  class(integracao), INTENT(IN) :: self
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