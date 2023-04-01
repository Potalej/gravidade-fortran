program main
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use simulacao
  use rungekutta4
  implicit none

  ! instanciamento da classe
  type(simular) :: Sim
  type(integracao) :: RK4
  
  ! arrays padrão
  real(pf), allocatable :: M(:), R(:,:), P(:,:), FSomas(:,:), R1(:,:,:), P1(:,:)
  ! quantidade de corpos
  integer :: N

  ! tempo de execução
  real :: t0, tf

  ! condições iniciais
  M = [100.0_pf,100.0_pf,80.0_pf, 10.0_pf]
  N = size(M)
  R = transpose(reshape((/   &
      157.14285714285714_pf, 42.85714285714285_pf, 64.28571428571428_pf, &
      -42.85714285714285_pf, -157.14285714285714_pf, -35.71428571428571_pf,   &
      -142.85714285714286_pf, 142.85714285714286_pf, -35.71428571428571_pf, &
      -10000.0_pf, -10000.0_pf, -10000.0_pf    /), (/3,N/)))

  P = transpose(reshape((/ &
      83.07091521588325_pf, -11.168958975925863_pf, 29.755473333965476_pf, & 
      -31.814968947314256_pf, 45.88086855383708_pf, -6.1955047860132035_pf, &
      -51.25594626856899_pf, -34.71190957791123_pf, -23.559968547952273_pf, &
      -10.0_pf, -10.0_pf, -10.0_pf /), (/3,N/)))   

    
  call cpu_time(t0)

  ! salva as informações
  call Sim % Iniciar (M, R, P)

  call Sim % rodar(800)

  call cpu_time(tf)

  print *
  print *, 'Tempo: ', tf - t0

end program