program main
  use condicoesArtigo
  use analise
  use simulacao
  use rungekutta4
  implicit none

  ! instanciamento da classe
  type(simular) :: Sim
  type(integracao) :: RK4
  
  ! arrays padrão
  real, allocatable :: M(:), R(:,:), P(:,:), FSomas(:,:), R1(:,:,:), P1(:,:)
  ! quantidade de corpos
  integer :: N

  ! tempo de execução
  real :: t0, tf

  ! condições iniciais
  M = [100.0,100.0,80.0]
  N = size(M)
  
  allocate(R(N,3))
  R(:,:) = 0
  
  allocate(P(N,3))
  P(:,:) = 0

  call gerar_condicionado(3, M, R, P, -100, 100, -10, 10, 100, 100)

  call exibir_informacoes(M, R, P)

end program