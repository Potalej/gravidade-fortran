module condicoesArtigo
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use angular
  use hamiltoniano
  use auxiliares

  implicit none
  private
  public condicionar, zerar_momentoAngular, zerar_energiaTotal, zerar_centroMassas, zerar_momentoLinear, gerar_condicionado

contains

  function gerar_vetores3d (N, min, max)

    integer, intent(in)     :: N, min, max
    real(pf), dimension(N,3)    :: gerar_vetores3d
    integer, dimension(N,3) :: ajuste

    ajuste(:,:) = min

    call random_seed()
    call random_number(gerar_vetores3d)

    ! agora condiciona no intervalo
    gerar_vetores3d = gerar_vetores3d * (max - min + 1) + ajuste
    
    ! arruma
    gerar_vetores3d = transpose(reshape(gerar_vetores3d, (/3,N/)))

  end function gerar_vetores3d

  function gerar_massas (N, min, max)

    integer, intent(in) :: N, min, max
    real(pf), dimension(N)  :: gerar_massas
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
    real(pf), intent(inout) :: posicoes(N,3), momentos(N,3), massas(N)

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
    real(pf), intent(inout)    :: posicoes(:,:), momentos(:,:), massas(:)
    integer                :: a
    real(pf)                :: momentoAngular_total(3), vetorRotacao(3), tensorInercia(3,3)

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
  subroutine zerar_energiaTotal (G, massas, posicoes, momentos)

    implicit none
    real(pf), intent(inout) :: G, posicoes(:,:), momentos(:,:), massas(:)
    real(pf)               :: EP, EC, fator

    ! calcula as energias
    EP = energia_potencial(G, massas, posicoes)
    EC = energia_cinetica(massas, momentos)

    ! calcula o fator
    fator = (-EP/EC)**0.5

    ! aplica sobre os momentos
    momentos = fator * momentos

  end subroutine zerar_energiaTotal

  ! zerando o centro de massas
  subroutine zerar_centroMassas (massas, posicoes)

    implicit none
    real(pf), intent(inout) :: massas(:), posicoes(:,:)
    real(pf), dimension(3)  :: rcm
    integer             :: a

    rcm = centro_massas(massas, posicoes)
    do a = 1, size(massas)
      posicoes(a,:) = posicoes(a,:) - rcm
    end do

  end subroutine zerar_centroMassas

  ! zerando o momento linear total
  subroutine zerar_momentoLinear (massas, momentos)

    implicit none
    real(pf), intent(inout) :: massas(:), momentos(:,:)
    real(pf), dimension(3) :: pcm ! análogo ao rcm
    integer            :: a

    ! usa o mesmo método porque a ideia é exatamente igual
    pcm = momentoLinear_total(momentos) / sum(massas)
    ! substitui
    do a = 1, size(massas)
      momentos(a,:) = momentos(a,:) - massas(a)*pcm
    end do

  end subroutine

  ! condiciona vetores já existentes
  subroutine condicionar (G, massas, posicoes, momentos)
    
    implicit none
    real(pf), intent(inout) :: G, posicoes(:,:), momentos(:,:), massas(:)

    ! zera o centro de massas
    call zerar_centroMassas(massas, posicoes)
    
    ! zera o momento linear
    call zerar_momentoLinear(massas, momentos)

    ! zera o momento angular
    call zerar_momentoAngular(massas, posicoes, momentos)

    ! zera a energia total
    call zerar_energiaTotal(G, massas, posicoes, momentos)

    call zerar_momentoAngular(massas, posicoes, momentos)
    
  end subroutine condicionar

  ! gerar valores aleatórios condicionados
  subroutine gerar_condicionado (G, N, massas, posicoes, momentos, minPos, maxPos, minMom, maxMom, minMassas, maxMassas)

    implicit none
    integer, intent(in) :: N, minPos, maxPos, minMom, maxMom, minMassas, maxMassas
    real(pf), intent(inout) :: G, massas(:), posicoes(:,:), momentos(:,:)

    ! gera os valores
    call gerarValores(N, massas, posicoes, momentos, minMassas, maxMassas, minPos, maxPos, minMom, maxMom)

    ! condiciona
    call condicionar(G, massas, posicoes, momentos)

  end subroutine gerar_condicionado

end module condicoesArtigo