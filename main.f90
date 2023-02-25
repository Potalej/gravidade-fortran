program main
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
  R = transpose(reshape((/   &
       100.0,  100.0, 100.0, &
      -100.0, -100.0, 0.0,   &
      -200.0,  200.0, 0.0    /), (/3,N/)))

  P = transpose(reshape((/ &
      70.0,  0.0,  0.0,    &
      0.0,  30.0,  10.0,   &
      0.0,  -5.0, -10.0   /), (/3,N/)))

  call cpu_time(t0)

  ! salva as informações
  call Sim % Iniciar (M, R, P)

  call Sim % rodar(1000)

  call cpu_time(tf)

  print *
  print *, 'Tempo: ', tf - t0

end program