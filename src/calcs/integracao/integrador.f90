! ************************************************************
!! INTEGRADOR NUMERICO
!
! Objetivos:
!   Estrutura basica de um integrador numerico. As configuracoes
!   gerais de um metodo devem vir aqui, como tamanho de passo,
!   dimensao, massas, se corrige ou nao, se colide ou nao, etc.
!
! Modificado:
!   10 de novembro de 2024
!
! Autoria:
!   oap
! 
MODULE integrador

  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  USE OMP_LIB
  USE funcoes_forca
  USE mecanica
  USE correcao
  USE colisao
  USE octree
  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao

  TYPE :: integracao
  
    ! m: Massas
    REAL(pf), ALLOCATABLE :: m(:), massasInvertidas(:,:)

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

    ! Vetor de forcas
    PROCEDURE(forcas_funcbase), POINTER, NOPASS :: forcas_funcao => NULL()

    ! Arvore octree
    TYPE(arvore_octo), ALLOCATABLE :: arvore

    CONTAINS
      PROCEDURE :: Iniciar, aplicarNVezes, metodo, forcas
  
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
!   14 de setembro de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE Iniciar (self, massas, G, h, potsoft, E0, J0, corrigir, corme, cormnt, colidir, colmodo, colmd, paralelo)
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

  ! Quantidade de particulas
  self % N = SIZE(massas)
  ! Massas
  ALLOCATE(self % m (self % N))
  self % m = massas
  
  ! vetor de massas invertidas
  ALLOCATE(self % massasInvertidas (self % N, self % dim))
  DO a = 1, self % N
    DO i = 1, self % dim
      self % massasInvertidas(a,i) = 1/(massas(a))
    END DO
  END DO

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
    self % forcas_funcao => forcas_par
  ELSE
    self % forcas_funcao => forcas_seq
  ENDIF
END SUBROUTINE Iniciar


FUNCTION forcas (self, R)
  IMPLICIT NONE
  class(integracao), INTENT(IN) :: self
  REAL(pf), DIMENSION(self % N, self % dim), INTENT(IN) :: R
  REAL(pf), DIMENSION(self % N, self % dim) :: forcas
  
  forcas = self % forcas_funcao(self%m, R, self%G, self%N, self%dim, self%potsoft, self%potsoft2)

END FUNCTION forcas

! ************************************************************
!! Aplicacao iterada do metodo
!
! Objetivos:
!   Aplica o metodo iterativamente N vezes.
!
! Modificado:
!   10 de novembro de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE aplicarNVezes (self, R, P, passos_antes_salvar)

  IMPLICIT NONE
  class (integracao), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P
  INTEGER, INTENT(IN) :: passos_antes_salvar
  REAL(pf)             :: E
  REAL(pf), DIMENSION(3) :: J
  ! Para cada passo
  INTEGER :: i
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

  ! Integrando
  DO i = 1, passos_antes_salvar
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

  ! Se estiver disposto a corrigir, calcula a energia total para ver se precisa
  IF (self%corrigir) THEN
    E = energia_total(self % G, self % m, R1, P1)
        
    IF (ABS(E - self%E0) >= self%corme) THEN
      
      t0_cor = omp_get_wtime()

      ! Correcao com energia total e momento angular total (desativado)
      ! CALL corrigir(self%corme,self%cormnt,self % G,self % m,R1,P1,corrigiu,E0,J0)

      ! Corrige somente a energia total
      CALL corrigir_apenas_energia(self % corme, self % cormnt, self % G, &
                                  self % m, R1, P1, corrigiu, self % E0, self % J0)
      
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

END MODULE integrador