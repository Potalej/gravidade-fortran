! Para facilitar na analise
! 

module analise
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use mecanica
  use auxiliares

  implicit none
  private
  public exibir_informacoes

contains

  ! exibe as informacoes
  subroutine exibir_informacoes (G, massas, posicoes, momentos)

    implicit none
    real(pf), intent(in) :: G, massas(:), posicoes(:,:), momentos(:,:)

    WRITE (*,*) 'Momento angular total: ', momento_angular_total(posicoes, momentos)
    WRITE (*,*) 'Energia total: ', energia_total(G, massas, posicoes, momentos)
    WRITE (*,*) 'Centro de massas: ', centro_massas(massas, posicoes)
    WRITE (*,*) 'Momento linear total: ', momentoLinear_total(momentos)

  end subroutine exibir_informacoes

end module analise