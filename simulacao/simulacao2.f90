module simulacao
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use hamiltoniano
  use angular
  ! use rungekutta4
  use rkf45
  use arquivos

  implicit none
  private
  public simular

  ! classe de simulacao
  type :: simular

  !> N: quantidade de corpos
  !> passos: quantidade de passos por integracao