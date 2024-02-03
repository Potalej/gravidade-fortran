! Estrutura basica de um integrador
! As configuracoes gerais da integracao devem vir aqui

module integrador

  use, intrinsic :: iso_fortran_env, only: pf=>real64
  implicit none
  private
  public integracao

  type :: integracao
    ! m: Massas
    real(pf), allocatable :: m(:)

    ! h: Passo de integracao
    ! G: Constante de gravitacao
    real(pf) :: h, G

    ! dim: Dimensao do problema
    ! N: Quantidade de partículas
    integer :: dim = 3, N

    ! Se vai ou nao corrigir
    logical :: corrigir = .FALSE.

    ! Se vai ou nao colidir
    logical :: colidir = .FALSE.

    ! vetores para aplicar a correcao
    real(pf), allocatable :: grads(:,:), gradsT(:,:), vetorCorrecao(:)
  end type integracao

end module integrador