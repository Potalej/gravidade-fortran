! ************************************************************
!! INTEGRADOR NUMERICO
!
! Objetivos:
!   Estrutura basica de um integrador numerico. As configuracoes
!   gerais de um metodo devem vir aqui, como tamanho de passo,
!   dimensao, massas, etc.
!
! Modificado:
!   29 de janeiro de 2026
!
! Autoria:
!   oap
! 
MODULE integrador
	
  USE tipos
  USE funcoes_forca
  USE json_utils_mod
  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao, iniciar_base

!> Esta eh a classe de integracao, da qual todos os metodos de 
!  integracao devem ser filhos.
  TYPE, ABSTRACT :: integracao
    ! m: Massas
    REAL(pf), ALLOCATABLE :: m(:)
    REAL(pf), ALLOCATABLE :: massasInvertidas(:,:)
    
    ! Massas iguais
    LOGICAL  :: mi
    REAL(pf) :: m_esc, m_inv, m2
    REAL(pf128) :: m_esc_128, m_inv_128, m2_128

    ! Vetor de forcas entre os corpos
    REAL(pf), ALLOCATABLE :: fs(:,:)

    ! h: Passo de integracao
    ! G: Constante de gravitacao
    ! potsoft: Softening do potencial
    REAL(pf) :: h, G, potsoft, potsoft2

    ! dim: Dimensao do problema
    ! N: Quantidade de partículas
    INTEGER :: dim = 3, N

    ! Se vai ou nao usar paralelizacao
    LOGICAL :: paralelo = .FALSE., gpu = .FALSE.

    ! Funcao de forcas (aceleracao)
    PROCEDURE(forcas_funcbase), POINTER, NOPASS    :: forcas_funcao    => NULL()
    PROCEDURE(forcas_mi_funcbase), POINTER, NOPASS :: forcas_mi_funcao => NULL()

    ! Se eh um metodo multipasso. 0 se nao eh, > 0 se for
    INTEGER :: multipasso = 0
    REAL(pf), ALLOCATABLE :: Ps_ant(:,:,:), Fs_ant(:,:,:)
    INTEGER :: mp_rb_idx ! indice do ring buffer para metodos multipasso

    CONTAINS
      PROCEDURE :: iniciar => iniciar_base
      PROCEDURE :: inicializar_massas, atualizar_constantes, &
                   forcas, metodo, metodo_mi, aplicar, &
                   mp_rb ! multipasso ring-buffer
                  
  END TYPE integracao

CONTAINS

! ************************************************************
!! Construtor da classe
!
! Objetivos:
!   Define o principal, salvando os valores e inicializando o
!   metodo.
!
! Modificado:
!   29 de janeiro de 2026
!
! Autoria:
!   oap
! 
SUBROUTINE iniciar_base (self, infos, timestep, massas)
  IMPLICIT NONE
  class(integracao), INTENT(INOUT) 			:: self
  TYPE(json_value), POINTER, INTENT(IN) :: infos
  REAL(pf), INTENT(IN), allocatable :: massas(:)
  REAL(pf), INTENT(IN) :: timestep
  LOGICAL :: encontrado

  !> Quantidade de particulas
  self % N = SIZE(massas)

  !> Inicializando as variaveis ligadas as massas dos corpos
  CALL self % inicializar_massas(infos, massas)

  !> Constantes
  self % G = json_get_float(infos, 'G') ! Const. de Gravitacao
  self % h = timestep                   ! Tamanho de passo
  !> Amortecedor
  self % potsoft = json_get_float(infos, 'integracao.amortecedor')
  self % potsoft2 = self % potsoft * self % potsoft

  !> Uso (ou nao) do paralelismo
  CALL json % get(infos, 'paralelo', self % paralelo)
  CALL json % get(infos, 'gpu', self % gpu, encontrado)
  IF (.NOT. encontrado) self % gpu = .FALSE.

  !> Determinando as funcoes de forca a se utilizar
  !> Modulo: funcoes_forca
  CALL inicializar_forcas(self%mi, self%paralelo, self%gpu, &
	  				  self%forcas_funcao, self%forcas_mi_funcao)

  ! Mesmo que nao seja um metodo multipasso, inicia o indice do ring buffer
  self % mp_rb_idx = 1
END SUBROUTINE iniciar_base

! ************************************************************
!! Inicializa as massas
!
! Modificado:
!   29 de janeiro de 2026
!
! Autoria:
!   oap
! 
SUBROUTINE inicializar_massas (self, infos, massas)
  CLASS(integracao), INTENT(INOUT) :: self
  TYPE(json_value), POINTER, INTENT(IN) :: infos
  REAL(pf), ALLOCATABLE :: massas(:)
  INTEGER :: a, i
  LOGICAL :: encontrado

  !> Massas iguais
  CALL json % get(infos, 'massas_iguais', self % mi, encontrado)
  IF (.NOT. encontrado) self % mi = .FALSE.

  !> Alocando vetor de massas
  IF (ALLOCATED(self % m)) DEALLOCATE(self % m)
  ALLOCATE(self % m (self % N))
  self % m = massas

  !> Define com maior precisao se tiver massas iguais.
  !> Se nao tiver, define o vetor de massas invertidas
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
    IF (ALLOCATED(self % massasInvertidas)) DEALLOCATE(self % massasInvertidas)
    ALLOCATE(self % massasInvertidas (self % N, self % dim))
    DO a = 1, self % N
      DO i = 1, self % dim
        self % massasInvertidas(a,i) = 1.0_pf/(massas(a))
      END DO
    END DO
  ENDIF
END SUBROUTINE inicializar_massas

! ************************************************************
!! Calculo das forcas conforme as massas
!
! Modificado:
!   27 de janeiro de 2026
!
! Autoria:
!   oap
! 
FUNCTION forcas (self, R)
  IMPLICIT NONE
  class(integracao), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self % N, self % dim), INTENT(IN) :: R
  REAL(pf), DIMENSION(self % N, self % dim) :: forcas
  
  IF (self % mi) THEN
    forcas = self % forcas_mi_funcao(R, self%G, self%N, self%dim, &
                    self%potsoft2)
  ELSE
    forcas = self % forcas_funcao(self % m, R, self%G, self%N, self%dim, &
                    self%potsoft2)
  ENDIF
END FUNCTION forcas

! ************************************************************
!! Metodo numerico
!
! Objetivos:
!   Aplicacao do metodo em si.
!
! Modificado:
!   18 de janeiro de 2026
!
! Autoria:
!   oap
!
SUBROUTINE metodo (self, R, P, FSomas)
  IMPLICIT NONE
  class(integracao), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  
  ! Cada integrador precisa ter um metodo definido, que substitui esta funcao vazia
  WRITE (*,*) 'OPS'

END SUBROUTINE metodo

! ************************************************************
!! Metodo numerico (massas iguais)
!
! Objetivos:
!   Aplicacao do metodo em si.
!
! Modificado:
!   18 de janeiro de 2026
!
! Autoria:
!   oap
!
SUBROUTINE metodo_mi (self, R, P, FSomas)
  IMPLICIT NONE
  class(integracao), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  
  ! Cada integrador precisa ter um metodo definido, que substitui esta funcao vazia
  WRITE (*,*) 'OPS'

END SUBROUTINE metodo_mi

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

! ************************************************************
!! Aplicacao do metodo
!
! Modificado:
!   18 de janeiro de 2026
!
! Autoria:
!   oap
!
SUBROUTINE aplicar (self, R, P)
  class (integracao), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P

  IF (.NOT. ALLOCATED(self % fs)) THEN
    ALLOCATE(self % fs(self % N, 3))
    self % fs = self % forcas(R)
  ENDIF

  IF (self % mi) THEN
    CALL self % metodo_mi(R, P, self % fs)
  ELSE
    CALL self % metodo(R, P, self % fs)
  ENDIF

  ! Atualiza o indice do ring buffer se for multipasso
  IF (self % multipasso > 0) self % mp_rb_idx = self % mp_rb()
END SUBROUTINE

! ************************************************************
!! (MULTIPASSOS) Indice do ring buffer
!
! Modificado:
!   29 de janeiro de 2026
!
! Autoria:
!   oap
!
PURE INTEGER FUNCTION mp_rb (self, k)
  CLASS (integracao), INTENT(IN) :: self
  INTEGER, OPTIONAL, INTENT(IN)  :: k

  IF (.NOT. PRESENT(k)) THEN
    mp_rb = modulo(self%mp_rb_idx-2, self%multipasso-1) + 1
  ELSE  
    mp_rb = modulo(self%mp_rb_idx+k-2, self%multipasso-1) + 1
  ENDIF
END FUNCTION mp_rb


END MODULE integrador
