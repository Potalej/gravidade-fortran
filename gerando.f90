program main_gerando

  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use simulacao
  use condicoesArtigo
  implicit none

  ! Instanciamento da classe
  type(simular) :: Sim

  ! Arrays padrao
  real(pf), allocatable :: M(:), R(:,:), P(:,:)
  ! Gravidade
  real(pf) :: G = 5.0_pf, h = 0.0001_pf
  ! Quantidade de corpos
  integer :: N = 10, passos = 100

  ! Tempo de execucao
  real :: t0, tf

  ! Quantidade total de passos
  integer :: qntd_total_passos = 5000

  ! Gera
  call gerar_condicionado (G, N, M, R, P, -500, 500, -100, 100, 10, 30)
  print *, 'gerou'

  ! Roda o metodo de Verlet
  call cpu_time(t0)
  call Sim % Iniciar (G, M, R, P, h, passos)
  call Sim % rodar_verlet(qntd_total_passos)
  call cpu_time(tf)
  print *, 'Tempo verlet: ', tf - t0

end program main_gerando