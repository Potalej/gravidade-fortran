! Integracao
! 
! Aqui consta a classe de integracao numerica via metodo de Runge-Kutta de 
! ordem 4.
! 

module rungekutta4
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use rungekutta
  use correcao
  use colisao
  use integrador

  implicit none
  private
  public integracao_rk4

  type, extends(integracao) :: integracao_rk4
    
    ! Base do Runge-Kutta
    type(RK) :: baseRK
    ! Modulos adicionais
    contains
      procedure :: Iniciar, metodo, aplicarNVezes

  end type integracao_rk4

contains

  ! Construtor da classe, para definir o principal
  subroutine Iniciar (self, massas, G, h, corrigir, colidir)
    implicit none
    class(integracao_rk4), intent(inout) :: self
    logical,intent(in) :: corrigir, colidir
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

    ! Se vai ou nao corrigir
    self % corrigir = corrigir

    ! Se vai ou nao colidir
    self % colidir = colidir

    ! alocando variaveis de correcao
    allocate(self%grads(10, 6*self%N))
    allocate(self%gradsT(6*self%N,10))
    allocate(self%vetorCorrecao(1:6*self%N))
    
    ! Inicia o base do RK 
    call self % baseRK % Iniciar(self % n, self % m, self % G, self % h)

  end subroutine Iniciar


  ! Método em si
  function metodo (self, R, P, FSomas)

    implicit none
    class(integracao_rk4), intent(in) :: self
    real(pf), dimension(self % N, self % dim), intent(in) :: R, P, FSomas
    real(pf), dimension(self % N, self % dim) :: R1, P1
    real(pf), dimension(2, self % N, self % dim) :: metodo

    ! componentes da integracao (kappas)
    real(pf), dimension(self % N, self % dim) :: k1, k2, k3, k4, fator

    ! faz a integracao sobre as equacoes x'
    k1 = P * self % baseRK % massasInvertidas
    k2 = k1 * self % baseRK % massasInvertidas
    k3 = k2 * self % baseRK % massasInvertidas
    k4 = k3 * self % baseRK % massasInvertidas

    ! fator para integracao
    fator = (self % h / 6) * (6*k1 + 3*self % h*k2 + self % h**2 * k3 + 0.25 * self % h**3 * k4)

    ! integra as posicoes
    R1 = R + fator

    ! integra os momentos
    P1 = P + self % h * FSomas

    metodo(1,:,:) = R1
    metodo(2,:,:) = P1

  end function metodo


  ! Aplicador do metodo com correcao (para aplicar varias vezes)
  subroutine aplicarNVezes (self, R, P, passos, E0, J0)

    implicit none
    class (integracao_rk4), intent(inout)                    :: self
    real(pf), dimension(self % N, self % dim), intent(inout) :: R, P
    integer, intent(in)                               :: passos
    real(pf), intent(in)                                  :: E0
    real(pf), dimension(3), intent(in)                    :: J0
    ! para cada passo
    integer :: i, qntd = 10
    ! para verificar se corrigiu
    logical :: corrigiu = .FALSE.
    ! para as forcas e passos pos-integracao
    real(pf), dimension (self % N, self % dim) :: F, R1, P1
    real(pf), dimension (2, self % N, self % dim) :: resultado  
    R1 = R
    P1 = P

    do i = 1, passos
      ! calcula as forcas
      F = self % baseRK % forcas (R1)
      ! aplicada o metodo
      resultado = self % metodo (R1, P1, F)
      
      R1 = resultado(1,:,:)
      P1 = resultado(2,:,:)

      ! aplica a correcao geral, se solicitado
      if (self % corrigir) then
        call corrigir(self % G, self % m, R1, P1,self%grads,self%gradsT,self%vetorCorrecao, corrigiu)
      end if

      ! se tiver colisoes, aplica
      if (self % colidir .AND. .NOT. corrigiu) then
        call verificar_e_colidir(self % m, R1, P1)
      end if

    end do

    R = R1
    P = P1

  end subroutine aplicarNVezes


end module rungekutta4