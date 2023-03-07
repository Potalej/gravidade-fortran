module condicoesArtigo
  use angular
  use hamiltoniano
  use auxiliares

  implicit none
  private
  public condicionar, zerar_momentoAngular, zerar_energiaTotal

contains

  function gerar_vetores3d (N, min, max)

    integer, intent(in)     :: N, min, max
    real, dimension(N,3)    :: gerar_vetores3d
    integer, dimension(N,3) :: ajuste

    ajuste(:,:) = min

    call random_seed()
    call random_number(gerar_vetores3d)

    ! agora condiciona no intervalo
    gerar_vetores3d = gerar_vetores3d * (max - min + 1) + ajuste

  end function gerar_vetores3d

  function gerar_massas (N, min, max)

    integer, intent(in) :: N, min, max
    real, dimension(N)  :: gerar_massas
    integer          :: ajuste(N)

    ajuste(:) = min

    call random_seed()
    call random_number(gerar_massas)

    ! agora condiciona no intervalo
    gerar_massas = gerar_massas * (max - min + 1) + ajuste

  end function gerar_massas

  subroutine gerarValores (N, massas, posicoes, momentos, minMassas, maxMassas, minPos, maxPos, minMom, maxMom)

    implicit none
    integer, intent(in)      :: N, minPos, maxPos, minMom, maxMom, minMassas, maxMassas
    real, intent(inout) :: posicoes(N,3), momentos(N,3), massas(N)

    ! gera massas
    massas = gerar_massas(N, minMassas, maxMassas)

    ! gera as posições
    posicoes = gerar_vetores3d(N, minPos, maxPos)

    ! gera os momentos
    momentos = gerar_vetores3d(N, minMom, maxMom)

  end subroutine gerarValores

  ! zerando o momento angular
  subroutine zerar_momentoAngular (massas, posicoes, momentos)
    
    implicit none
    real, intent(inout)    :: posicoes(:,:), momentos(:,:), massas(:)
    integer                :: a
    real                :: momentoAngular_total(3), vetorRotacao(3), tensorInercia(3,3)

    ! calcula o momento angular total
    momentoAngular_total = angular_geral (posicoes,momentos)

    ! calcular o tensor de inércia
    tensorInercia = tensor_inercia_geral(massas, posicoes)

    ! calcula o vetor de rotação (resolve sistema linear)
    vetorRotacao = sistema_linear3 (tensorInercia, - momentoAngular_total)

    ! percorre os corpos
    do a = 1, size(posicoes,1)
      ! produto vetorial da posição pelo vetor de rotação
      momentos(a,:) = momentos(a,:) + massas(a) * produto_vetorial(posicoes(a,:), vetorRotacao)
    end do

  end subroutine zerar_momentoAngular

  ! zerando a energia total
  subroutine zerar_energiaTotal (massas, posicoes, momentos)

    implicit none
    real, intent(inout) :: posicoes(:,:), momentos(:,:), massas(:)
    real                :: EP, EC, fator

    ! calcula as energias
    EP = energia_potencial(massas, posicoes)
    EC = energia_cinetica(massas, momentos)

    ! calcula o fator
    fator = (-EP/EC)**0.5

    ! aplica sobre os momentos
    momentos = fator * momentos

  end subroutine zerar_energiaTotal

  ! subroutine condicionar (N, minPos, maxPos, minMom, maxMom, minMassas, maxMassas)
  subroutine condicionar (massas, posicoes, momentos)
    implicit none

    ! integer, intent(in)      :: N, minPos, maxPos, minMom, maxMom, minMassas, maxMassas
    real, intent(inout) :: posicoes(:,:), momentos(:,:), massas(:)

    ! gera os valores
    ! gerarValores(N, massas, posicoes, momentos, minMassas, maxMassas, minPos, maxPos, minMom, maxMom)

    ! zera o momento angular

    ! zera o centro de massas

    ! zera o momento linear

    ! zera a energia total


  end subroutine condicionar

end module condicoesArtigo