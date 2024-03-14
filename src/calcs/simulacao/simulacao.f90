module simulacao
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use OMP_LIB

  ! Auxiliares
  use mecanica
  use auxiliares
  use arquivos
  ! Metodos de integracao numerica
  use rungekutta4
  use verlet

  implicit none
  private
  public simular

  ! Classe de simulacao
  type :: simular
  
    !> N: Quantidade de corpos
    !> dim: Dimensao do problema
    integer :: N, dim = 3, passos_antes_salvar

    !> h: Tamanho do passo de integracao
    !> G: Constante de gravitacao universal
    !> E0: Energia total inicial
    !> mtot: Massa total do sistema
    real(pf) :: h, G, E0, mtot
    
    !> M: Massas do sistema
    !> R: Posicoes das particulas
    !> P: Momento linear das particulas
    !> Jtot: Momento angular total do sistema
    !> Ptot: Momento linear total do sistema
    !> Rcm: Centro de massas do sistema
    real(pf), allocatable :: M(:), R(:,:), P(:,:), Jtot(:), Ptot(:), Rcm(:)
    
    real(pf), dimension(3) :: J0

    ! Configuracoes
    logical :: corrigir=.FALSE., colidir=.FALSE.

    !> Arquivo
    type(arquivo) :: Arq

    contains
      procedure :: Iniciar, rodar_verlet, rodar_rk4
  end type

contains
  
  ! Construtor da classe, para definir o principal
  subroutine Iniciar (self, G, M, R0, P0, h, passos_antes_salvar)

    class(simular), intent(inout) :: self

    real(pf), allocatable :: M(:), R0(:,:), P0(:,:)
    real(pf) :: G, h
    integer :: a, i, passos_antes_salvar

    self % passos_antes_salvar = passos_antes_salvar
    
    ! Salva o tamanho dos passos
    self % h = h

    ! Salva a gravidade
    self % G = G

    ! Salva as massas
    self % M = M  
    
    ! Quantidade corpos no sistema
    self % N = size(M)

    ! Massa total do sistema
    self % mtot = sum(self % M)

    ! Salva as posicoes e os momentos
    self % R = R0
    self % P = P0

    ! Salva a dimensao
    self % dim = size(R0,2)

    ! Salva momento linear total inicial
    self % Ptot = momentoLinear_total(self % P)

    ! Salva a energia inicial
    self % E0 = energia_total(G, self % M, self % R, self % P)

    ! Salva o momento angular inicial
    self % J0 = momento_angular_total(self % R, self % P)

    ! Salva o centro de massas inicial
    self % Rcm = centro_massas(self % M, self % R)

  end subroutine

  ! Simulacao com o metodo Velocity Verlet (Simpletico)
  subroutine rodar_verlet (self, qntdPassos)
    implicit none
    class(simular), intent(inout) :: self
    integer, intent(in) :: qntdPassos
    ! iterador e variavel de tempo que sera o nome do arquivo
    integer :: i, t
    ! Integrador
    type(integracao_verlet) :: integrador
    ! Escritor de arquivos
    type(arquivo) :: Arq
    
    real(pf), dimension(2, self % N, self % dim) :: resultado
    real(pf), dimension(self % N, self % dim) :: R1, P1

    real(pf) :: t0, tf, tempo_total = 0.0_pf

    ! Cria o arquivo onde sera salvo
    call Arq % criar(2, self % N, self % dim)

    ! Salva as infos de cabecalho
    call Arq % escrever_cabecalho(self % h, self % G, self % M)

    ! Instanciamento do integrador
    WRITE (*,'(a)') 'INTEGRACAO NUMERICA'
    WRITE (*,*) ' > rodando com ', qntdPassos, ' passos'
    WRITE (*, '(a)') '  > instanciando o metodo velocity-verlet'
    call integrador % Iniciar(self % M, self % G, self % h, self%corrigir, self%colidir)

    ! Condicoes iniciais
    R1 = self % R
    P1 = self % P

    call Arq % escrever((/R1, P1/))

    ! Roda
    WRITE (*, '(a)') '  > iniciando simulacao...'
    do i = 1, qntdPassos

      ! timer
      t0 = omp_get_wtime()

      ! Integracao
      call integrador % aplicarNVezes(R1, P1, self % passos_antes_salvar, self % E0, self % J0)

      ! timer
      tf = omp_get_wtime()
      tempo_total = tempo_total + tf - t0

      call Arq % escrever((/R1, P1/))

      if (mod(i, 500) == 0) then
        WRITE (*,*) '     -> Passo:', i, ' / Energia:', energia_total(self % G, self % M, R1, P1)
      end if

    end do
    
    WRITE (*, '(a)') '  > simulacao encerrada!'
    WRITE (*,*)

    call Arq % fechar()

  end subroutine rodar_verlet

  ! Simulacao com o metodo de Runge-Kutta de ordem 4 (RK4)
  subroutine rodar_rk4 (self, qntdPassos)

    implicit none
    class(simular), intent(inout) :: self
    integer, intent(in) :: qntdPassos
    ! Iterador e variavel de tempo que sera o nome do arquivo
    integer :: i, t
    ! Integrador
    type(integracao_rk4) :: integrador
    ! Escritor de arquivos
    type(arquivo) :: Arq
    
    real(pf), dimension(2, self % N, self % dim) :: resultado
    real(pf), dimension(self % N, self % dim) :: R1, P1

    ! Instanciamento do integrador    
    call integrador % Iniciar(self % M, self % G, self % h, self%corrigir, self%colidir)

    ! Cria o arquivo onde ficara salvo
    call Arq % criar(1, self % N, self % dim)

    ! Salva as infos de cabecalho
    call Arq % escrever_cabecalho(self % h, self % G, self % M)

    ! Condicoes iniciais
    R1 = self % R
    P1 = self % P

    call Arq % escrever((/R1, P1/))

    ! Roda
    do i = 1, qntdPassos

      call integrador % aplicarNVezes(R1, P1, self % passos_antes_salvar, self % E0, self % J0)

      call Arq % escrever((/R1, P1/))

      if (mod(i, 500) == 0) then
        WRITE (*,*) 'energia:', energia_total(self % G, self % M, R1, P1)
        WRITE (*,*) 'passo: ', i
      end if

    end do 

    call Arq % fechar()

  end subroutine rodar_rk4

end module simulacao