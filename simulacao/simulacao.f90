module simulacao
  use hamiltoniano
  implicit none
  private
  public simular

  ! classe de simulação
  type :: simular
  
  !> N: Quantidade de corpos
  !> passos: Quantidade de passos por integração
  !> dim: Dimensão do problema
  integer :: N, passos = 100, dim = 3

  !> h: Tamanho do passo de integração
  !> G: Constante de gravitação universal
  !> E0: Energia total inicial
  !> mtot: Massa total do sistema
  real :: h = 0.05, G = 1.0, E0, mtot
  
  !> M: Massas do sistema
  !> R: Posições das partículas
  !> P: Momento linear das partículas
  !> Jtot: Momento angular total do sistema
  !> Ptot: Momento linear total do sistema
  real, allocatable :: M(:), R(:,:), P(:,:), Jtot(:), Ptot(:)

  contains
    procedure :: Iniciar

  end type

contains
  
  ! Construtor da classe, para definir o principal
  subroutine Iniciar (self, M, R0, P0)

    class(simular), intent(inout) :: self

    real, allocatable :: M(:), R0(:,:), P0(:,:)
    
    integer :: a, i
    
    ! salva as massas
    self % M = M  
    
    ! quantidade corpos no sistema
    self % N = size(M)

    ! massa total do sistema
    self % mtot = sum(self % M)

    ! salvas as posições e os momentos
    self % R = R0
    self % P = P0

    ! salva a dimensão
    self % dim = size(R0,1)

    ! salva momento linear total inicial
    self % Ptot = [(0, i = 1, self % dim)]
    ! percorre os corpos
    do a = 1, self % N
      self % Ptot = self % Ptot + self % P(a,:)
    end do

    ! salva a energia inicial
    self % E0 = energia_total(self % M, self % R, self % P)
  
  end subroutine

end module simulacao