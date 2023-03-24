program main
  use condicoesArtigo
  use arquivos
  use analise
  use simulacao
  use rungekutta4
  implicit none

  ! instanciamento da classe
  type(simular) :: Sim
  type(integracao) :: RK4
  
  ! arrays padrão
  ! real, allocatable :: M(:), R(:,:,:), P(:,:,:), FSomas(:,:), R1(:,:,:), P1(:,:)
  real, allocatable :: M(:), R(:,:), P(:,:), FSomas(:,:), R1(:,:,:), P1(:,:)
  ! quantidade de corpos
  integer :: N = 0

  ! tempo de execução
  real :: t0, tf

  ! condições iniciais
  N = 3 

  allocate(M(N))
  
  allocate(R(N,3))
  R(:,:) = 0
  
  allocate(P(N,3))
  P(:,:) = 0

  ! EXEMPLO DE GERAÇÃO CONDICIONADA
  call gerar_condicionado(N, M, R, P, -100, 100, -10, 10, 100, 100)
  call exibir_informacoes(M, R, P)

  ! salva as informações
  call Sim % Iniciar (M, R, P)

  ! integra
  call Sim % rodar(100)

end program