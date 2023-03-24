! Integração
! 
! Aqui consta a classe de integração numérica via método de Runge-Kutta de 
! ordem 4.
! 

module rungekutta4

  use hamiltoniano
  use angular

  implicit none
  private
  public integracao

  type :: integracao
  
    ! m: Massas
    ! massasInvertidas : Matriz com o inverso das massas para facilitar a integração dos momentos
    real, allocatable :: m(:), massasInvertidas(:,:)

    ! h: Passo de integração
    ! G: Constante de gravitação
    real :: h = 0.05, G = 3.0

    ! dim: Dimensão do problema
    ! N: Quantidade de partículas
    integer :: dim = 3, N

    contains
      procedure :: Iniciar, metodo, forcas, aplicarNVezes

  end type

contains

  ! Construtor da classe, para definir o principal
  subroutine Iniciar (self, massas)
    implicit none
    class(integracao), intent(inout) :: self
    real, allocatable :: massas(:)
    integer :: a, i

    ! quantidade de partículas
    self % N = size(massas)
    ! massas
    allocate(self % m (self % N))
    self % m = massas
    ! vetor de massas invertidas
    allocate(self % massasInvertidas (self % N, self % dim))
    do a = 1, self % N
      do i = 1, self % dim
        self % massasInvertidas(a,i) = 1/(massas(a))
      end do 
    end do
  end subroutine Iniciar


  ! Matriz de forças
  function forcas (self, R)
    implicit none
    class(integracao), intent(in) :: self
    real, dimension(self % N, self % dim), intent(in) :: R
    real, dimension(self % dim) :: Fab
    integer :: a, b
    real :: distancia
    real, dimension(self % N, self % dim) :: forcas
    
    forcas(:,:) = 0

    do a = 2, self % N
      do b = 1, a - 1
        ! distância entre os corpos
        distancia = norm2(R(b,:) - R(a,:))**3
        ! força entre os corpos a e b
        Fab = - self % G * self % m(a) * self % m(b) * (R(b,:) - R(a,:))/distancia
        ! Adiciona na matriz
        forcas(a,:) = forcas(a,:) - Fab
        forcas(b,:) = forcas(b,:) + Fab
      end do
    end do

  end function forcas

  
  ! Método em si
  function metodo (self, R, P, FSomas)

    implicit none
    class(integracao), intent(in) :: self
    real, dimension(self % N, self % dim), intent(in) :: R, P, FSomas
    real, dimension(self % N, self % dim) :: R1, P1
    real, dimension(2, self % N, self % dim) :: metodo

    ! componentes da integração (kappas)
    real, dimension(self % N, self % dim) :: k1, k2, k3, k4, fator

    ! faz a integração sobre as equações x'
    k1 = P * self % massasInvertidas
    k2 = k1 * self % massasInvertidas
    k3 = k2 * self % massasInvertidas
    k4 = k3 * self % massasInvertidas

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
  function aplicarNVezes (self, R, P, passos, E0, J0)

    implicit none
    class (integracao), intent(in)                    :: self
    real, dimension(self % N, self % dim), intent(in) :: R, P
    integer, intent(in)                               :: passos
    real, intent(in)                                  :: E0
    real, dimension(3), intent(in)                    :: J0
    ! para cada passo
    integer :: i
    ! para as forças e passos pós-integração
    real, dimension (self % N, self % dim) :: F, R1, P1
    real, dimension (2, self % N, self % dim) :: resultado , aplicarNVezes
    R1 = R
    P1 = P

    do i = 1, passos
      ! calcula as forças
      F = self % forcas (R1)
      ! aplicada o método
      resultado = self % metodo (R1, P1, F)

      R1 = resultado(1,:,:)
      P1 = resultado(2,:,:)
      
      ! aplica a correção de energia (ainda nao funciona direito)
      ! call energia_correcao(self % m, R1, P1, E0, self % G)

      ! aplica a correção de momento angular
      call angular_correcao(self % m, R1, P1, J0)

    end do

    aplicarNVezes(1,:,:) = R1
    aplicarNVezes(2,:,:) = P1

  end function aplicarNVezes


end module rungekutta4