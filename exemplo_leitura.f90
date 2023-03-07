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
  real, allocatable :: M(:), R(:,:,:), P(:,:,:), FSomas(:,:), R1(:,:,:), P1(:,:)
  ! quantidade de corpos
  integer :: N

  ! tempo de execução
  real :: t0, tf

  call ler_csv('20230307_10.csv', M, R, P)

  print *, 'massas: ', M
  print *, ''
  print *, 'posições: ', R
  print *, ''
  print *, 'momentos: ', P

end program