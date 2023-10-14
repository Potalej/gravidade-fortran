program main
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use OMP_LIB
  use simulacao
  use condicoesArtigo
  implicit none

  ! instanciamento da classe
  type(simular) :: Sim_rk4, Sim_verlet, Sim_corrigir
  
  ! arrays padrao
  real(pf), allocatable :: M(:), R(:,:), P(:,:)
  ! gravidade
  real(pf) :: G = 10.0_pf, h = 0.001_pf
  ! quantidade de corpos
  integer :: N = 10, passos = 100

  ! tempo de execucao
  real :: t0, tf

  ! quantidade total de passos
  integer :: qntd_total_passos = 50000


  call gerar_condicionado (G, N, M, R, P, -2000, 2000, -10, 10, 100, 800)
  print *, 'gerou'

  ! Roda ao contrario
  t0 = omp_get_wtime()
  ! call cpu_time(t0)
  Sim_verlet % corrigir = .FALSE.
  Sim_verlet % colidir = .FALSE.
  call Sim_verlet % Iniciar (G, M, R, P, -h, passos)
  call Sim_verlet % rodar_verlet(qntd_total_passos)
  ! call cpu_time(tf)
  tf = omp_get_wtime() 
  print *, 'Tempo Verlet: ', tf - t0
  print *, 'Tempo por passo: ', (tf-t0)/(qntd_total_passos*passos)
  print *
  ! print *, 'Tempo Verlet ao contrario: ', tf - t0 

  ! Roda o Verlet
  t0 = omp_get_wtime()
  ! call cpu_time(t0)
  Sim_verlet % corrigir = .TRUE.
  Sim_verlet % colidir = .FALSE.
  call Sim_verlet % Iniciar (G, M, R, P, h, passos)
  call Sim_verlet % rodar_verlet(qntd_total_passos)
  call cpu_time(tf)
  print *
  tf = omp_get_wtime() 
  print *, 'Tempo Verlet: ', tf - t0
  print *, 'Tempo por passo: ', (tf-t0)/(qntd_total_passos*passos)


end program