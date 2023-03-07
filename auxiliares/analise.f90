! Para facilitar na análise
! 

module analise

  use angular
  use hamiltoniano
  use auxiliares

  implicit none
  private
  public exibir_informacoes

contains

  ! exibe as informações
  subroutine exibir_informacoes (massas, posicoes, momentos)

    implicit none
    real, intent(in) :: massas(:), posicoes(:,:), momentos(:,:)

    print *, 'Momento angular total: ', angular_geral(posicoes, momentos)
    print *, 'Energia total: ', energia_total(massas, posicoes, momentos)
    print *, 'Centro de massas: ', centro_massas(massas, posicoes)
    print *, 'Momento linear total: ', momentoLinear_total(momentos)

  end subroutine exibir_informacoes

end module analise