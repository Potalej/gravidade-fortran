! Integração
! 
! Aqui consta a classe de integração numérica via método de Runge-Kutta de 
! ordem 4.
! 

module rungekutta4
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use rungekutta
  use hamiltoniano
  use angular
  use correcao

  implicit none
  private
  public integracao


  type :: integracao

    ! base do runge-kutta
    type(RK) :: baseRK
  
    ! m: Massas
    real(pf), allocatable :: m(:)

    ! h: Passo de integração
    ! G: Constante de gravitação
    real(pf) :: h, G

    ! dim: Dimensão do problema
    ! N: Quantidade de partículas
    integer :: dim = 3, N

    ! vetores para aplicar a correcao
    real(pf), allocatable :: grads(:,:), gradsT(:,:), vetorCorrecao(:)

    contains
      procedure :: Iniciar, metodo, aplicarNVezes

  end type

contains

  ! Construtor da classe, para definir o principal
  subroutine Iniciar (self, massas, G, h)
    implicit none
    class(integracao), intent(inout) :: self
    real(pf), allocatable :: massas(:)
    real(pf)              :: G, h
    integer :: a, i

    ! quantidade de partículas
    self % N = size(massas)
    ! massas
    allocate(self % m (self % N))
    self % m = massas

    ! gravidade
    self % G = G
    ! passo
    self % h = h

    ! inicia o método
    call self % baseRK % Iniciar(self % n, self % m, self % G, self % h)

    ! alocando variaveis de correcao
    allocate(self%grads(10, 6*self%N))
    allocate(self%gradsT(6*self%N,10))
    allocate(self%vetorCorrecao(1:6*self%N))

  end subroutine Iniciar

  ! Método em si
  function metodo (self, R, P, FSomas)

    implicit none
    class(integracao), intent(in) :: self
    real(pf), dimension(self % N, self % dim), intent(in) :: R, P, FSomas
    real(pf), dimension(self % N, self % dim) :: R1, P1
    real(pf), dimension(2, self % N, self % dim) :: metodo

    ! componentes da integração (kappas)
    real(pf), dimension(self % N, self % dim) :: k1, k2, k3, k4, fator

    ! faz a integração sobre as equações x'
    k1 = P * self % baseRK % massasInvertidas
    k2 = k1 * self % baseRK % massasInvertidas
    k3 = k2 * self % baseRK % massasInvertidas
    k4 = k3 * self % baseRK % massasInvertidas

    ! fator para integração
    fator = (self % h / 6) * (6*k1 + 3*self % h*k2 + self % h**2 * k3 + 0.25 * self % h**3 * k4)

    ! integra as posições
    R1 = R + fator

    ! integra os momentos
    P1 = P + self % h * FSomas

    metodo(1,:,:) = R1
    metodo(2,:,:) = P1

  end function metodo


  ! Aplicador do método com correção (para aplicar várias vezes)
  subroutine aplicarNVezes (self, R, P, passos, E0, J0)

    implicit none
    class (integracao), intent(inout)                    :: self
    real(pf), dimension(self % N, self % dim), intent(inout) :: R, P
    integer, intent(in)                               :: passos
    real(pf), intent(in)                                  :: E0
    real(pf), dimension(3), intent(in)                    :: J0
    ! para cada passo
    integer :: i
    ! para as forças e passos pós-integração
    real(pf), dimension (self % N, self % dim) :: F, R1, P1
    real(pf), dimension (2, self % N, self % dim) :: resultado  
    R1 = R
    P1 = P

    do i = 1, passos
      ! calcula as forças
      F = self % baseRK % forcas (R1)
      ! aplicada o método
      resultado = self % metodo (R1, P1, F)

      R1 = resultado(1,:,:)
      P1 = resultado(2,:,:)
      
      ! aplica a correção de energia (ainda nao funciona direito)
      ! call energia_correcao(self % m, R1, P1, E0, self % G)

      ! aplica a correção de momento angular
      ! call angular_correcao(self % m, R1, P1, J0)

      ! aplica a correcao geral
      call corrigir(self % G, self % m, R1, P1,self%grads,self%gradsT,self%vetorCorrecao)
    end do

    R = R1
    P = P1

    ! aplicarNVezes(1,:,:) = R1
    ! aplicarNVezes(2,:,:) = P1

  end subroutine aplicarNVezes


end module rungekutta4