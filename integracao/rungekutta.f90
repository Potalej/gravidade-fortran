! Integração
! 
! Aqui consta o básico da integração numérica com o método de Runge-Kutta.
! Os métodos básicos aqui servem apenas para construir outras classes que
! usem de fato o método, como RK4 e RKF45.

module rungekutta
  use, intrinsic :: iso_fortran_env, only: pf=>real64

  implicit none
  private
  public RK

  type :: RK

    ! m: Massas
    ! massasInvertidas : Matriz como inverso das massas para facilitar a integração dos momentos
    real(pf), allocatable :: m(:), massasInvertidas(:,:)

    ! h: Passo de integração
    ! G: Constante de gravitação
    real(pf) :: h, G

    ! dim: dimensão do problema
    ! N: quantidade de partículas
    integer :: dim = 3, N

    contains
      procedure :: Iniciar, forcas

  end type

contains

  ! Construtor da classe, para definir o principal
  subroutine Iniciar (self, N, massas, G, h)
    implicit none
    class(RK), intent(inout) :: self
    real(pf), allocatable :: massas(:)
    integer, intent(inout) :: N
    real(pf)               :: G, h
    integer :: a, i

    ! quantidade de particulas
    self % N = N
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

    ! gravidade
    self % G = G
    ! passo
    self % h = h
  end subroutine Iniciar  

  ! Matriz de forças
  function forcas (self, R)
    implicit none
    class(RK), intent(in) :: self
    real(pf), dimension(self % N, self % dim), intent(in) :: R
    real(pf), dimension(self % dim) :: Fab
    integer :: a, b
    real(pf) :: distancia
    real(pf), dimension(self % N, self % dim) :: forcas
    
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


end module rungekutta