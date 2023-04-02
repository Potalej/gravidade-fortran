module simulacao
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use hamiltoniano
  use angular
  ! use rungekutta4
  use rkf45
  use arquivos

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
  real(pf) :: h = 0.01_pf, G = 30.0_pf, E0, mtot
  
  !> M: Massas do sistema
  !> R: Posições das partículas
  !> P: Momento linear das partículas
  !> Jtot: Momento angular total do sistema
  !> Ptot: Momento linear total do sistema
  real(pf), allocatable :: M(:), R(:,:), P(:,:), Jtot(:), Ptot(:)
  
  real(pf), dimension(3) :: J0

  !> Arquivo
  type(arquivo) :: Arq

  contains
    procedure :: Iniciar, rodar

  end type

contains
  
  ! Construtor da classe, para definir o principal
  subroutine Iniciar (self, M, R0, P0)

    class(simular), intent(inout) :: self

    real(pf), allocatable :: M(:), R0(:,:), P0(:,:)
    
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
    self % dim = size(R0,2)

    ! salva momento linear total inicial
    self % Ptot = [(0, i = 1, self % dim)]
    ! percorre os corpos
    do a = 1, self % N
      self % Ptot = self % Ptot + self % P(a,:)
    end do

    ! salva a energia inicial
    self % E0 = energia_total(self % M, self % R, self % P)

    ! salva o momento angular inicial
    self % J0 = angular_geral(self % R, self % P)
  
  end subroutine


  ! Para simular de fato uma determinada quantidade de passos
  subroutine rodar (self, qntdPassos)

    implicit none
    class(simular), intent(inout) :: self
    integer, intent(in) :: qntdPassos
    ! iterador e variável de tempo que será o nome do arquivo
    integer :: i, t
    
    real(pf), dimension(2, self % N, self % dim) :: resultado
    real(pf), dimension(self % N, self % dim) :: R1, P1

    ! instanciamento do método de integração
    type(integracao) :: RK4
    call RK4 % Iniciar(self % M, self % G, self % h)

    ! cria o arquivo onde ficará salvo
    call self % Arq % criar(1, self % N, self % dim)

    ! salva as massas
    call self % Arq % escrever_massas(self % M)

    ! condições iniciais
    R1 = self % R
    P1 = self % P

    call self % Arq % escrever((/R1, P1/))

    ! roda
    do i = 1, qntdPassos

      resultado = RK4 % aplicarNVezes(R1, P1, self % passos, self % E0, self % J0)
      R1 = resultado(1,:,:)
      P1 = resultado(2,:,:)

      call self % Arq % escrever((/R1, P1/))

    end do 

    call self % Arq % fechar()

  end subroutine rodar

end module simulacao