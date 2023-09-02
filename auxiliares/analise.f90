! Para facilitar na analise
! 

module analise
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use angular
  use hamiltoniano
  use auxiliares

  implicit none
  private
  public exibir_informacoes

contains

  ! exibe as informacoes
  subroutine exibir_informacoes (G, massas, posicoes, momentos)

    implicit none
    real(pf), intent(in) :: G, massas(:), posicoes(:,:), momentos(:,:)

    print *, 'Momento angular total: ', angular_geral(posicoes, momentos)
    print *, 'Energia total: ', energia_total(G, massas, posicoes, momentos)
    print *, 'Centro de massas: ', centro_massas(massas, posicoes)
    print *, 'Momento linear total: ', momentoLinear_total(momentos)

  end subroutine exibir_informacoes

end module analise