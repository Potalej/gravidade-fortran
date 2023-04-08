program main
  use, intrinsic :: iso_fortran_env, only: pf=>real64
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
  real(pf), allocatable :: M(:), R(:,:), P(:,:), FSomas(:,:), R1(:,:,:), P1(:,:)
  ! quantidade de corpos
  integer :: N = 0

  ! tempo de execução
  real :: t0, tf

  ! condições iniciais
  ! N = 250 
  ! N = 1000
  N = 200  

  allocate(M(N))
  
  allocate(R(N,3))
  R(:,:) = 0.0_pf
  
  allocate(P(N,3))
  P(:,:) = 0.0_pf

  ! EXEMPLO DE GERAÇÃO CONDICIONADA
  ! call gerar_condicionado(N, M, R, P, -10000, 10000, -1, 1, 100, 300)
  call gerar_condicionado(N, M, R, P, -5000, 5000, -10, 10, 1000, 3000)
  call exibir_informacoes(M, R, P)  

  call cpu_time(t0)

  ! salva as informações
  call Sim % Iniciar (M, R, P)

  ! integra
  call Sim % rodar(3000) 

  call cpu_time(tf)

  print *
  print *, 'Tempo: ', tf - t0

end program