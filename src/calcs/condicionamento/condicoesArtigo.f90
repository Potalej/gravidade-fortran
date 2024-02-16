! Funcoes para condicionamento de particulas
! 
! Funcoes
! 
! = function gerar_vetores3d (N, min, max)
! Gera N vetores tridimensionais limitados no retangulo [min, max]3
! 
! = function gerar_massas (N, min, max)
! Gera vetor de massas das particulas.
! 
! = subroutine gerarVaores (N, massas, posicoes, momentos, minMassas, maxMassas, minPos, maxPos, minMom, maxMom)
! Gera condicoes iniciais aleatorias.
! 
! = subroutine zerar_momentoAngular (massas, posicoes, momentos)
! Zera o momento angular dadas as massas, posicoes e momentos.
! 
! = subroutine zerar_energiaTotal (G, massas, posicoes, momentos)
! Zera a energia total dada a gravidade, massas, posicoes e momentos.
! 
! = subroutine zerar_centroMassas (massas, posicoes)
! Faz a translacao dos corpos para zerar o centro de massas.
! 
! = subroutine zerar_momentoLinear (massas, momentos)
! Zera o momento linear dadas as massas e momentos.
! 
! = subroutine condicionar (G, massas, posicoes, momentos)
! Faz o condicionamento dadas as condicoes.
! 
! = gerar_condicionado (G, N, massas, posicoes, momentos, minPos, maxPos, minMom, maxMom, minMassas, maxMassas)
! Gera condicoes iniciais aleatorias e depois condiciona.
! 

module condicoesArtigo
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use mecanica
  use auxiliares

  IMPLICIT NONE

contains

  function gerar_vetores3d (N, min, max)

    integer, intent(in)     :: N
    real(pf), intent(in)    :: min, max
    real(pf), dimension(N,3)    :: gerar_vetores3d
    integer, dimension(N,3) :: ajuste

    ajuste(:,:) = min

    call random_seed()
    call random_number(gerar_vetores3d)

    ! Agora condiciona no intervalo
    gerar_vetores3d = gerar_vetores3d * (max - min + 1) + ajuste
    
    ! Arruma
    gerar_vetores3d = transpose(reshape(gerar_vetores3d, (/3,N/)))

  end function gerar_vetores3d

  function gerar_massas (N, min, max)

    integer, intent(in) :: N
    real(pf), intent(in) :: min, max
    real(pf), dimension(N)  :: gerar_massas
    integer          :: ajuste(N)

    ajuste(:) = min

    call random_seed()
    call random_number(gerar_massas)
    
    ! Agora condiciona no intervalo
    gerar_massas = gerar_massas * (max - min + 1) + ajuste

  end function gerar_massas

  subroutine gerarValores (N, massas, posicoes, momentos, int_posicoes, int_momentos, int_massas)

    implicit none
    integer, intent(in)      :: N
    real(pf), dimension(2), intent(in) :: int_posicoes, int_momentos, int_massas
    real(pf), intent(inout) :: posicoes(N,3), momentos(N,3), massas(N)

    ! Gera massas
    WRITE (*,*) '    * gerando massas'
    massas = gerar_massas(N, int_massas(1), int_massas(2))

    ! Gera as posições
    WRITE (*,*) '    * gerando posicoes'
    posicoes = gerar_vetores3d(N, int_posicoes(1), int_posicoes(2))

    ! Gera os momentos
    WRITE (*,*) '    * gerando momentos'
    momentos = gerar_vetores3d(N, int_momentos(1), int_momentos(2))

  end subroutine gerarValores

  ! Condicionando o momento angular
  subroutine condicionar_momentoAngular (J, massas, posicoes, momentos)
    
    implicit none
    real(pf), intent(inout) :: posicoes(:,:), momentos(:,:), massas(:)
    integer                 :: a
    real(pf)                :: momentoAngular_total(3), vetorRotacao(3), tensorInercia(3,3), J(3)

    ! Calcula o momento angular total
    momentoAngular_total = momento_angular_total (posicoes,momentos)

    ! Calcular o tensor de inercia
    tensorInercia = tensor_inercia_geral(massas, posicoes)

    ! Calcula o vetor de rotacao (resolve sistema linear)
    vetorRotacao = sistema_linear3 (tensorInercia, - momentoAngular_total + J)

    ! Percorre os corpos
    do a = 1, size(posicoes,1)
      ! Produto vetorial da posicao pelo vetor de rotacao
      momentos(a,:) = momentos(a,:) + massas(a) * produto_vetorial(posicoes(a,:), vetorRotacao)
    end do

  end subroutine condicionar_momentoAngular

  ! Condicionando a energia total
  subroutine condicionar_energiaTotal (H, G, massas, posicoes, momentos)
    ! Outra forma de H>0 seria somente aplicar o fator ((H-EP)/EC)**0.5

    implicit none
    real(pf), intent(inout) :: H, G, posicoes(:,:), momentos(:,:), massas(:)
    real(pf)                :: EP, EC, fator

    ! Calcula as energias
    EP = energia_potencial(G, massas, posicoes)
    EC = energia_cinetica(massas, momentos)

    ! Calcula o fator para zerar
    fator = (-EP/EC)**0.5
      
    ! Aplica sobre os momentos
    momentos = fator * momentos

    ! Se nao for zero, aplica homotetia nas posicoes
    IF (H .NE. 0) THEN
      fator = 1/(H/EP + 1)
      posicoes = fator * posicoes
    END IF

  end subroutine condicionar_energiaTotal

  ! Zerando o centro de massas
  subroutine zerar_centroMassas (massas, posicoes)

    implicit none
    real(pf), intent(inout) :: massas(:), posicoes(:,:)
    real(pf), dimension(3)  :: rcm
    integer                 :: a

    rcm = centro_massas(massas, posicoes)
    do a = 1, size(massas)
      posicoes(a,:) = posicoes(a,:) - rcm
    end do

  end subroutine zerar_centroMassas

  ! Condicionando o momento linear total
  subroutine condicionar_momentoLinear (P, massas, momentos)

    implicit none
    real(pf), intent(inout) :: massas(:), momentos(:,:)
    real(pf)                :: P(3)
    real(pf), dimension(3)  :: pcm ! analogo ao rcm
    integer                 :: a

    ! Usa o mesmo metodo porque a ideia eh exatamente igual
    pcm = (momentoLinear_total(momentos) - P)/ sum(massas)
    ! Substitui
    do a = 1, size(massas)
      momentos(a,:) = momentos(a,:) - massas(a)*pcm
    end do

  end subroutine condicionar_momentoLinear

  ! Condiciona vetores ja existentes
  subroutine condicionar (G, massas, posicoes, momentos, H, J, P)
    
    implicit none
    real(pf), intent(inout) :: G, posicoes(:,:), momentos(:,:), massas(:)
    REAL(pf) :: H, J(3), P(3)
    REAL(pf) :: erro_0, erro_1, ENERGIA, LINEAR(3), ANGULAR(3) ! Para medir a taxa de erro
    INTEGER :: i = 0

    ! Zera o centro de massas
    call zerar_centroMassas(massas, posicoes)
        
    ! Condiciona o momento linear
    call condicionar_momentoLinear(P, massas, momentos)

    ! Condiciona o momento angular
    call condicionar_momentoAngular(J, massas, posicoes, momentos)

    ! ! Condiciona a energia total
    call condicionar_energiaTotal(H, G, massas, posicoes, momentos)

    ! Calculo dos erros
    ENERGIA = energia_total(G,massas,posicoes,momentos) - H
    LINEAR = momentoLinear_total(momentos) - P
    ANGULAR = momento_angular_total(posicoes,momentos) - J
    erro_1 = NORM2((/ENERGIA,LINEAR(1),LINEAR(2),LINEAR(3),ANGULAR(1),ANGULAR(2),ANGULAR(3)/))

    IF (erro_1 >= 0.00005) THEN
      erro_0 = ABS(erro_1) + 1.0 ! Para entrar no loop abaixo
      
      DO WHILE (ABS(erro_0 - erro_1) >= 0.00005 .AND. i <= 10)
        i = i + 1
        erro_0 = erro_1

        ! Condiciona o momento linear
        call condicionar_momentoLinear(P, massas, momentos)

        ! Condiciona o momento angular
        call condicionar_momentoAngular(J, massas, posicoes, momentos)

        ! ! Condiciona a energia total
        call condicionar_energiaTotal(H, G, massas, posicoes, momentos)

        ! Calculo dos erros
        ENERGIA = energia_total(G,massas,posicoes,momentos) - H
        LINEAR = momentoLinear_total(momentos) - P
        ANGULAR = momento_angular_total(posicoes,momentos) - J
        erro_1 = NORM2((/ENERGIA,LINEAR(1),LINEAR(2),LINEAR(3),ANGULAR(1),ANGULAR(2),ANGULAR(3)/))
      END DO
    END IF
    WRITE (*,*) ' > condicionamento aplicado ', i, ' vezes para obter o erro ', erro_1
  end subroutine condicionar

  ! Gerar valores aleatorios condicionados
  subroutine gerar_condicionado (G, N, massas, posicoes, momentos, int_posicoes, int_momentos, int_massas, H, J, P)

    implicit none
    integer, intent(in)     :: N
    real(pf), dimension(2), intent(in) :: int_posicoes, int_momentos, int_massas
    real(pf), intent(inout) :: G, H, J(3), P(3)
    real(pf), intent(inout), allocatable :: massas(:), posicoes(:,:), momentos(:,:)

    WRITE (*,'(a)') "GERACAO DAS CONDICOES INICIAIS"
    
    ! Gera os valores
    WRITE (*,'(a)') '  > gerando valores...'
    allocate(massas (N))
    allocate(posicoes (N, 3))
    allocate(momentos (N, 3))
    call gerarValores(N, massas, posicoes, momentos, int_posicoes, int_momentos, int_massas)

    ! Condiciona
    WRITE (*,'(a)') '  > condicionando...'
    call condicionar(G, massas, posicoes, momentos, H, J, P)

    ! Exibe as integrais primeiras do sistema
    WRITE (*,*) '    * H   =', energia_total(G,massas,posicoes,momentos) 
    WRITE (*,*) '    * Rcm =', centro_massas(massas,posicoes) 
    WRITE (*,*) '    * P   =', momentoLinear_total(momentos) 
    WRITE (*,*) '    * J   =', momento_angular_total(posicoes,momentos) 

    WRITE (*,'(a)') '  > condicoes iniciais geradas!'
    WRITE (*,*)

  end subroutine gerar_condicionado

end module condicoesArtigo