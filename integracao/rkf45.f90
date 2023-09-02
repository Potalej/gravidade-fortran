! Integracao via Runge-Kutta-Fehlberg (RKF45)
! 
! Aqui consta a classe de integracao numerica via metodo de Runge-Kutta-Fehlberg de 
! ordem 4/5.
! 

module rkf45
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use rungekutta
  use hamiltoniano
  use angular
  use correcao
  use colisao
  use integrador

  implicit none
  private
  public integracao_rkf45

  type, extends(integracao) :: integracao_rkf45

    ! base do runge-kutta
    type(RK) :: baseRK
  
    ! constantes do metodo
    real, dimension(6) :: c = (/0, 1/4, 3/8, 12/13, 1, 1/2/)
    real, dimension(6) :: a2 = (/1/4, 0, 0, 0, 0, 0/)    
    real, dimension(6) :: a3 = (/3/32, 9/32, 0, 0, 0, 0/)
    real, dimension(6) :: a4 = (/1932/2197, -7200/2197, 7296/2197, 0, 0, 0/)
    real, dimension(6) :: a5 = (/ 439/216, -8, 3680/513, -845/4104, 0, 0 /)
    real, dimension(6) :: a6 = (/ -8/27, 2, -3544/2565, 1859/4104, -11/40, 0 /)
    real, dimension(6) :: b1 = (/16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55 /) 
    real, dimension(6) :: b2 = (/25/216, 0, 1408/2565, 2197/4104, -1/5, 0 /)
    ! real, dimension(6) :: CT = (/ 1.0/360.0, 0.0, -128.0/4275.0, -2197.0/75240.0, 1.0/50.0, 2.0/55.0 /)
    real(pf),dimension(6)::CT=(/1.0_pf/360.0_pf,0.0_pf,-128.0_pf/4275.0_pf,-2197.0_pf/75240.0_pf,1.0_pf/50.0_pf,2.0_pf/55.0_pf/)
    real(pf) :: epsilon = 0.00001_pf

    contains
      procedure :: Iniciar, metodo, aplicarNVezes, tolerancia, aplicarNVezesControleAutomatico

  end type

contains

  ! Construtor da classe, para definir o principal
  subroutine Iniciar (self, massas, G, h, corrigir, colidir)
    implicit none
    class(integracao_rkf45), intent(inout) :: self
    real(pf), allocatable :: massas(:)
    real(pf)              :: G, h
    logical,intent(in) :: corrigir, colidir
    integer :: a, i

    ! quantidade de particulas
    self % N = size(massas)
    ! massas
    allocate(self % m (self % N))
    self % m = massas

    ! gravidade
    self % G = G
    ! passo
    self % h = h

    ! Se vai ou nao corrigir
    self % corrigir = corrigir

    ! Se vai ou nao colidir
    self % colidir = colidir

    ! inicia o metodo
    call self % baseRK % Iniciar(self % n, self % m, self % G, self % h)

    ! alocando variaveis de correcao
    allocate(self%grads(10, 6*self%N))
    allocate(self%gradsT(6*self%N,10))
    allocate(self%vetorCorrecao(1:6*self%N))
  end subroutine Iniciar

  ! Metodo em si
  function metodo (self, R, P, FSomas, controle)

    implicit none
    class(integracao_rkf45), intent(inout)                      :: self
    real(pf), dimension(self % N, self % dim), intent(in) :: R, P, FSomas
    real(pf), dimension(self % N, self % dim)             :: R1, P1
    real(pf), dimension(2, self % N, self % dim)          :: metodo
    logical                                               :: controle

    ! componentes da integracao (kappas)
    real(pf), dimension(self % N, self % dim) :: k1, k2, k3, k4, k5, k6, fator
    real(pf) :: TE

    do while (.true.)

      ! metodo RKF45
      k1 = self % h * (P * self % baseRK % massasInvertidas)
      k2 = self % h * (k1 + (self%a2(1)*k1) * self%baseRK%massasInvertidas)
      k3 = self % h * (k1 + (self%a3(1)*k1 + self%a3(2)*k2) * self%baseRK%massasInvertidas)
      k4 = self % h * (k1 + (self%a4(1)*k1 + self%a4(2)*k2 + self%a4(3)*k3) * self%baseRK%massasInvertidas)
      k5 = self % h * (k1 + (self%a5(1)*k1 + self%a5(2)*k2 + self%a5(3)*k3 + self%a5(4)*k4) * self%baseRK%massasInvertidas)
      k6 = k1+(self%a6(1)*k1 + self%a6(2)*k2 + self%a6(3)*k3 + self%a6(4)*k4 + self%a6(5)*k5)*self%baseRK%massasInvertidas
      k6 = self % h * k6

      ! calcula o erro tolerado
      if (controle) then
        call self % tolerancia (TE, k1, k2, k3, k4, k5, k6)
      else
        TE = 0.0
      end if
      
      if (TE > self % epsilon) then
        ! print *, 'caiu aqui: ', TE, ' / ', self % h
        self % h = self % h * (self%epsilon/(2*TE))**(0.25)
        self % baseRK % h = self % h
      else
        fator = self%b1(1) * k1 + self%b1(2) * k2 + self%b1(3) * k3 + self%b1(4) + k4 + self%b1(5) * k5 + self%b1(6) * k6
        
        ! integra as posicoes
        R1 = R + fator
        ! integra os momentos
        P1 = P + self % h * FSomas

        metodo(1,:,:) = R1
        metodo(2,:,:) = P1

        if (TE > 0) then
          self % h = 0.001_pf
          self % baseRK % h = self % h
        end if
        
        exit 
      end if 
    
    end do 

  end function metodo


  ! Aplicador do metodo com correcao (para aplicar varias vezes)
  subroutine aplicarNVezes (self, R, P, passos, E0, J0)

    implicit none
    class (integracao_rkf45), intent(inout)                    :: self
    real(pf), dimension(self % N, self % dim), intent(inout) :: R, P
    integer, intent(in)                               :: passos
    real(pf), intent(in)                                  :: E0
    real(pf), dimension(3), intent(in)                    :: J0
    ! para cada passo
    integer :: i
    ! para as forcas e passos pos-integracao
    real(pf), dimension (self % N, self % dim) :: F, R1, P1
    real(pf), dimension (2, self % N, self % dim) :: resultado 
    ! para verificar se corrigiu
    logical :: corrigiu = .FALSE.

    R1 = R
    P1 = P

    do i = 1, passos
      ! calcula as forcas
      F = self % baseRK % forcas (R1)
      ! aplicada o metodo
      resultado = self % metodo (R1, P1, F, .false.)

      R1 = resultado(1,:,:)
      P1 = resultado(2,:,:)

      ! aplica a correcao geral
      call corrigir(self % G, self % m, R1, P1,self%grads,self%gradsT,self%vetorCorrecao, corrigiu)

      if (.NOT. corrigiu) then
        call verificar_e_colidir (self % m, R1, P1)
      endif

    end do

    R = R1
    P = P1

  end subroutine aplicarNVezes

  ! Aplicador do metodo com correcao (para aplicar varias vezes) e 
  ! com controle de passo atraves do parametro h0
  subroutine aplicarNVezesControleAutomatico (self, R, P, passos, E0, J0, h0)
    
    implicit none
    class (integracao_rkf45), intent(inout)                    :: self
    real(pf), dimension(self % N, self % dim), intent(inout) :: R, P
    integer, intent(in)                               :: passos
    real(pf), intent(in)                                  :: E0, h0
    real(pf), dimension(3), intent(in)                    :: J0
    ! para cada passo
    integer :: i
    ! para as forças e passos pós-integração
    real(pf), dimension (self % N, self % dim) :: F, R1, P1
    real(pf), dimension (2, self % N, self % dim) :: resultado
    real(pf) :: h_soma
    
    R1 = R
    P1 = P

    do i = 1, passos
      
      h_soma = 0.0_pf

      loop_while: do while (h_soma < h0)

        ! calcula as forcas
        F = self % baseRK % forcas (R1)
        ! aplicada o metodo
        resultado = self % metodo (R1, P1, F, .true.)

        R1 = resultado(1,:,:)
        P1 = resultado(2,:,:)
        
        h_soma = h_soma + self % h

      end do loop_while

    end do

    R = R1
    P = P1

  end subroutine aplicarNVezesControleAutomatico

  subroutine tolerancia (self, TE, k1, k2, k3, k4, k5, k6)

    implicit none
    class(integracao_rkf45), intent(in) :: self
    real(pf), dimension(self % N, self % dim) :: TE_soma, k1, k2, k3, k4, k5, k6
    real(pf), dimension(3) :: TE_vet
    real(pf) :: TE

    TE_soma = k1*self%CT(1) + k2*self%CT(2) + k3*self%CT(3) + k4*self%CT(4) + k5*self%CT(5) + k6*self%CT(6)
    TE_vet = PRODUCT(TE_soma, DIM = 1)**(0.5)
    TE = MAXVAL(TE_vet)

  end subroutine tolerancia

end module rkf45