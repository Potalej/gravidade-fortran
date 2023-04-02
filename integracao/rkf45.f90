! Integração
! 
! Aqui consta a classe de integração numérica via método de Runge-Kutta de 
! ordem 4.
! 

module rkf45
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use rungekutta
  use hamiltoniano
  use angular

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

    ! constantes do método
    real, dimension(6) :: c = (/0, 1/4, 3/8, 12/13, 1, 1/2/)
    real, dimension(6) :: a2 = (/0, 0, 0, 0, 0, 0/)    
    real, dimension(6) :: a3 = (/3/32, 9/32, 0, 0, 0, 0/)
    real, dimension(6) :: a4 = (/1932/2197, -7200/2197, 7296/2197, 0, 0, 0/)
    real, dimension(6) :: a5 = (/ 439/216, -8, 3680/513, -845/4104, 0, 0 /)
    real, dimension(6) :: a6 = (/ -8/27, 2, -3544/2565, 1859/4104, -11/40, 0 /)
    real, dimension(6) :: b1 = (/16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55 /) 
    real, dimension(6) :: b2 = (/25/216, 0, 1408/2565, 2197/4104, -1/5, 0 /)

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
  end subroutine Iniciar

  ! Método em si
  function metodo (self, R, P, FSomas)

    implicit none
    class(integracao), intent(in) :: self
    real(pf), dimension(self % N, self % dim), intent(in) :: R, P, FSomas
    real(pf), dimension(self % N, self % dim) :: R1, P1
    real(pf), dimension(2, self % N, self % dim) :: metodo

    ! componentes da integração (kappas)
    real(pf), dimension(self % N, self % dim) :: k1, k2, k3, k4, k5, k6, fator

    k1 = P * self % baseRK % massasInvertidas
    k2 = k1+self%h * (self%a2(1)*k1) * self%baseRK%massasInvertidas
    k3 = k1+self%h * (self%a3(1)*k1 + self%a3(2)*k2) * self%baseRK%massasInvertidas
    k4 = k1+self%h * (self%a4(1)*k1 + self%a4(2)*k2 + self%a4(3)*k3) * self%baseRK%massasInvertidas
    k5 = k1+self%h * (self%a5(1)*k1 + self%a5(2)*k2 + self%a5(3)*k3 + self%a5(4)*k4) * self%baseRK%massasInvertidas
    k6 = k1+self%h * (self%a6(1)*k1 + self%a6(2)*k2 + self%a6(3)*k3 + self%a6(4)*k4 + self%a6(5)*k5)*self%baseRK%massasInvertidas

    fator = self%b1(1) * k1 + self%b1(2) * k2 + self%b1(3) * k3 + self%b1(4) + k4 + self%b1(5) * k5 + self%b1(6) * k6

    ! integra as posições
    R1 = R + self % h * fator

    ! integra os momentos
    P1 = P + self % h * FSomas

    metodo(1,:,:) = R1
    metodo(2,:,:) = P1

  end function metodo


  ! Aplicador do método com correção (para aplicar várias vezes)
  function aplicarNVezes (self, R, P, passos, E0, J0)

    implicit none
    class (integracao), intent(in)                    :: self
    real(pf), dimension(self % N, self % dim), intent(in) :: R, P
    integer, intent(in)                               :: passos
    real(pf), intent(in)                                  :: E0
    real(pf), dimension(3), intent(in)                    :: J0
    ! para cada passo
    integer :: i
    ! para as forças e passos pós-integração
    real(pf), dimension (self % N, self % dim) :: F, R1, P1
    real(pf), dimension (2, self % N, self % dim) :: resultado , aplicarNVezes
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

    end do

    aplicarNVezes(1,:,:) = R1
    aplicarNVezes(2,:,:) = P1

  end function aplicarNVezes


end module rkf45