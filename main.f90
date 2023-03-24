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
      157.14285714285714, 42.85714285714285, 64.28571428571428, &
      -42.85714285714285, -157.14285714285714, -35.71428571428571,   &
      -142.85714285714286, 142.85714285714286, -35.71428571428571    /), (/3,N/)))

  P = transpose(reshape((/ &
      83.07091521588325, -11.168958975925863, 29.755473333965476, & 
      -31.814968947314256, 45.88086855383708, -6.1955047860132035, &
      -51.25594626856899, -34.71190957791123, -23.559968547952273 /), (/3,N/)))

  call cpu_time(t0)

  ! salva as informações
  call Sim % Iniciar (M, R, P)

  call Sim % rodar(500)

  call cpu_time(tf)

  print *
  print *, 'Tempo: ', tf - t0

end program