! Integracao via metodo de Verlet
! 
! Aqui consta a lasse de integracao numerica via metodo de Verlet (simpletico)
!

module verlet
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use OMP_LIB
  use correcao
  use colisao
  use integrador

  implicit none
  private
  public integracao_verlet

  type, extends(integracao) :: integracao_verlet

    contains
      procedure :: Iniciar, metodo, aplicarNVezes, Forcas

  end type
  
contains

  ! Construtor da classe
  subroutine Iniciar (self, massas, G, h, corrigir, colidir)
    implicit none
    class(integracao_verlet), intent(inout) :: self
    real(pf), allocatable :: massas(:)
    real(pf)              :: G, h
    logical,intent(in) :: corrigir, colidir
    integer :: a, i

    ! Quantidade de particulas
    self % N = size(massas)
    ! Massas
    allocate(self % m (self % N))
    self % m = massas

    ! gravidade
    self % G = G
    ! Passo
    self % h = h

    ! Se vai ou nao corrigir
    self % corrigir = corrigir

    ! Se vai ou nao colidir
    self % colidir = colidir  

    ! Alocando variaveis de correcao
    allocate(self%grads(4, 6*self%N))
    allocate(self%gradsT(6*self%N,4))
    allocate(self%vetorCorrecao(1:6*self%N))

  end subroutine Iniciar

  ! Calculo das forcas
  function forcas (self, R)
    implicit none
    class(integracao_verlet), intent(in) :: self
    real(pf), dimension(self % N, self % dim), intent(in) :: R
    real(pf), dimension(self % dim) :: Fab, dif
    integer :: a, b, thread, threads, qntdPorThread
    real(pf) :: distancia
    real(pf), dimension(self % N, self % dim) :: forcas
    real(pf), dimension(self%N,self%N,self%dim) :: MF
    
    forcas(:,:) = 0
    MF(:,:,:) = 0
    
    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(a, b, distancia, Fab) SHARED(MF, self, R) REDUCTION(+:forcas)
      !$OMP DO
      do a = 2, self%N
        do b = 1, a-1
          ! distancia entre os corpos
          distancia = norm2(R(b,:) - R(a,:))**3 + self % h
          ! forca entre os corpos a e b
          Fab = - self % G * self % m(a) * self % m(b) * (R(b,:) - R(a,:))/distancia
          
          ! Adiciona na matriz
          ! MF(a,b,:) = -Fab
          ! MF(b,a,:) =  Fab
          forcas(a,:) = forcas(a,:) - Fab
          forcas(b,:) = forcas(b,:) + Fab
        end do
      end do
      !$OMP END DO
    !$OMP END PARALLEL

    ! do a = 1, self%N
    !   forcas(a,1) = SUM(MF(a,:,1))
    !   forcas(a,2) = SUM(MF(a,:,2))
    !   forcas(a,3) = SUM(MF(a,:,3))
    ! enddo
          

  end function forcas

  ! Metodo em si
  function metodo (self, R, P, FSomas_ant)

    implicit none
    class(integracao_verlet), intent(in) :: self
    real(pf), dimension(self%N, self%dim), intent(in) :: R, P, FSomas_ant
    real(pf), dimension(self%N, self%dim) :: R1, P1, FSomas_prox
    real(pf), dimension(3, self%N, self%dim) :: metodo
    integer :: a

    ! Integrando as posicoes
    do a = 1, self % N
      R1(a,:) = R(a,:) + self % h * P(a,:) / self % m(a) + 0.5*(self%h**2)*FSomas_ant(a,:)/self%m(a)
    end do

    ! Calcula as novas forcas
    FSomas_prox = self%forcas(R1)

    ! Integrando as velocidades
    P1 = P + 0.5*self%h*(FSomas_ant + FSomas_prox)

    metodo(1,:,:) = R1
    metodo(2,:,:) = P1
    metodo(3,:,:) = FSomas_prox

  end function metodo

  ! Aplicacao do metodo (para aplicar varias vezes)
  subroutine aplicarNVezes (self, R, P, passos, E0, J0)

    implicit none
    class (integracao_verlet), intent(inout) :: self
    real(pf), dimension(self%N, self%dim), intent(inout) :: R, P
    integer, intent(in) :: passos
    real(pf), intent(in) :: E0
    real(pf), dimension(3), intent(in) :: J0
    ! Para cada passo
    integer :: i
    ! para verificar se corrigiu
    logical :: corrigiu = .FALSE.
    ! Para as forcas e passos pos-integracao
    real(pf), dimension(self%N, self%dim) :: R1, P1, FSomas_ant
    real(pf), dimension(3, self%N, self%dim) :: resultado

    ! Salvando as primeiras posicoes e momentos
    R1 = R
    P1 = P

    ! Calcula as forcas
    FSomas_ant = self%forcas(R)

    ! Integrando
    do i = 1, passos
      ! Aplica o metodo
      resultado = self % metodo(R1, P1, FSomas_ant)

      ! se tiver colisoes, aplica
      if (self % colidir .AND. .NOT. corrigiu) then
        call verificar_e_colidir(self % m, R1, P1)
      end if

      ! Separando as variaveis
      R1 = resultado(1,:,:)
      P1 = resultado(2,:,:)
      FSomas_ant = resultado(3,:,:)

      if (self%corrigir) then
        call corrigir(self % G, self % m, R1, P1,self%grads,self%gradsT,self%vetorCorrecao, corrigiu)
      endif

    end do

    R = R1
    P = P1

  end subroutine aplicarNVezes

end module verlet