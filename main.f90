program main
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use simulacao
  use condicoesArtigo
  implicit none

  ! instanciamento da classe
  type(simular) :: Sim_rk4, Sim_verlet, Sim_corrigir
  
  ! arrays padrao
  real(pf), allocatable :: M(:), R(:,:), P(:,:)
  ! gravidade
  real(pf) :: G = 5.0_pf, h = 0.0001_pf
  ! quantidade de corpos
  integer :: N = 10, passos = 100

  ! tempo de execucao
  real :: t0, tf

  ! quantidade total de passos
  integer :: qntd_total_passos = 5000

  ! condições iniciais
  ! M = [100.0_pf,100.0_pf,80.0_pf]
  ! N = size(M)
  ! R = transpose(reshape((/   &
  !     157.14285714285714_pf, 42.85714285714285_pf, 64.28571428571428_pf, &
  !     -42.85714285714285_pf, -157.14285714285714_pf, -35.71428571428571_pf,   &
  !     -142.85714285714286_pf, 142.85714285714286_pf, -35.71428571428571_pf /), (/3,N/)))

  ! P = transpose(reshape((/ &
  !     83.07091521588325_pf, -11.168958975925863_pf, 29.755473333965476_pf, & 
  !     -31.814968947314256_pf, 45.88086855383708_pf, -6.1955047860132035_pf, &
  !     -51.25594626856899_pf, -34.71190957791123_pf, -23.559968547952273_pf /), (/3,N/)))   

  ! M = [100.0_pf, 100.0_pf, 100.0_pf, 100.0_pf]
  ! N = size(M)
  ! R = transpose(reshape((/   &
  !     -10.0_pf, -10.0_pf, 0.0_pf, &
  !      10.0_pf,  10.0_pf, 0.0_pf, &
  !      10.0_pf, -10.0_pf, 0.0_pf, &
  !     -10.0_pf,  10.0_pf, 0.0_pf /), (/3,N/)))
  ! P = transpose(reshape((/   &
  !      1062.0551580910137_pf, -515.8553625013496_pf, -60.68886617662932_pf, &
  !     -576.5442286779789_pf,   637.2330948546082_pf, -60.68886617662932_pf, &
  !     -394.47763014809084_pf, -1304.8106227975313_pf, 60.68886617662918_pf, &
  !     -91.03329926494408_pf,   1183.4328904442725_pf, 60.68886617662945_pf /), (/3,N/)))  

  ! M = [10, 11, 12]
  ! N = size(M)
  ! R = transpose(reshape((/   &
  !     -30.0_pf, 0.0_pf, 0.0_pf,&
  !        0.0_pf, 0.0_pf, 10.0_pf, &
  !        0.0_pf, 0.0_pf, 0.0_pf /), (/3,N/)))

  ! P = transpose(reshape((/   &
  !     1000.0_pf, 1000.0_pf, 0.0_pf, &
  !     -1000.0_pf, -10.0_pf, 0.0_pf,  &
  !     0.0_pf, 0.0_pf, 0.0_pf /), (/3,N/)))  

  ! call condicionar (G, M, R, P)

  call gerar_condicionado (G, N, M, R, P, -500, 500, -100, 100, 10, 30)
  print *, 'gerou'

  ! ! CONDICOES INICIAIS DE UMA LEMNISCATA
  ! G = 1.0_pf
  ! M = [1.0_pf, 1.0_pf, 1.0_pf]
  ! N = size(M)
  ! R = transpose(reshape((/   &
  !     -0.97000436_pf, 0.24308753_pf, 0.0_pf,&
  !        0.0_pf, 0.0_pf, 0.0_pf, &
  !      0.97000436_pf, -0.24308753_pf, 0.0_pf /), (/3,N/)))

  ! P = transpose(reshape((/   &
  !     0.4662036850_pf, 0.4323657300_pf, 0.0_pf, &
  !     -0.93240737_pf, -0.86473146_pf, 0.0_pf,  &
  !     0.4662036850_pf, 0.4323657300_pf, 0.0_pf /), (/3,N/)))  

  ! call condicionar (G, M, R, P)

  ! Roda o RK4 corrigindo
  call cpu_time(t0)
  Sim_corrigir % corrigir = .TRUE.
  Sim_corrigir % colidir = .TRUE.
  call Sim_corrigir % Iniciar (G, M, R, P, h, passos)
  call Sim_corrigir % rodar_rk4(qntd_total_passos)
  call cpu_time(tf)
  print *
  print *, 'Tempo RK4 corrigido: ', tf - t0 

  ! Roda o Verlet
  call cpu_time(t0)
  call Sim_verlet % Iniciar (G, M, R, P, h, passos)
  call Sim_verlet % rodar_verlet(qntd_total_passos)
  call cpu_time(tf)
  print *
  print *, 'Tempo Verlet: ', tf - t0 

  ! Roda o RK4
  call cpu_time(t0)
  call Sim_rk4 % Iniciar (G, M, R, P, h, passos)
  call Sim_rk4 % rodar_rk4(qntd_total_passos)
  call cpu_time(tf)
  print *
  print *, 'Tempo RK4: ', tf - t0 


end program