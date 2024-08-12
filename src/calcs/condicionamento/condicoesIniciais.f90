! *****************************************************************
!! CONDICIONAMENTO DE VALORES INICIAIS
!
! Objetivos:
!   Funcoes para o condicionamento de particulas.
! 
! Modificado:
!   15 de marco de 2024
! 
! Autoria:
!   oap
! 
MODULE condicoesIniciais
  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  USE mecanica
  USE auxiliares
  IMPLICIT NONE
CONTAINS

! ************************************************************
!! Gera vetores 3d
!
! Objetivos:
!   Gera vetores 3d utilizados para as posicoes e velocidades.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
FUNCTION gerar_vetores3d (N, min, max)

  INTEGER, INTENT(IN)     :: N
  REAL(pf), INTENT(IN)    :: min, max
  REAL(pf), DIMENSION(N,3)    :: gerar_vetores3d
  INTEGER, DIMENSION(N,3) :: ajuste

  ajuste(:,:) = min

  CALL RANDOM_SEED()
  CALL RANDOM_NUMBER(gerar_vetores3d)

  ! Agora condiciona no intervalo
  gerar_vetores3d = gerar_vetores3d * (max - min + 1) + ajuste
  
  ! Arruma
  gerar_vetores3d = TRANSPOSE(RESHAPE(gerar_vetores3d, (/3,N/)))

END FUNCTION gerar_vetores3d

! ************************************************************
!! Gera vetor de massas
!
! Objetivos:
!   Gera vetor de massas conforme um intervalo informado.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
FUNCTION gerar_massas (N, min, max)

  INTEGER, INTENT(IN)    :: N
  REAL(pf), INTENT(IN)   :: min, max
  REAL(pf), DIMENSION(N) :: gerar_massas
  INTEGER                :: ajuste(N)

  ajuste(:) = min

  CALL RANDOM_SEED()
  CALL RANDOM_NUMBER(gerar_massas)
  
  ! Agora condiciona no intervalo
  gerar_massas = gerar_massas * (max - min + 1) + ajuste

END FUNCTION gerar_massas

! ************************************************************
!! Gera valores
!
! Objetivos:
!   Gera valores iniciais aleatorios.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE gerarValores (N, massas, posicoes, momentos, int_posicoes, int_momentos, int_massas)

  IMPLICIT NONE
  INTEGER, INTENT(IN)      :: N
  REAL(pf), DIMENSION(2), INTENT(IN) :: int_posicoes, int_momentos, int_massas
  REAL(pf), INTENT(INOUT) :: posicoes(N,3), momentos(N,3), massas(N)

  ! Gera massas
  WRITE (*,*) '    * gerando massas'
  massas = gerar_massas(N, int_massas(1), int_massas(2))

  ! Gera as posições
  WRITE (*,*) '    * gerando posicoes'
  posicoes = gerar_vetores3d(N, int_posicoes(1), int_posicoes(2))

  ! Gera os momentos
  WRITE (*,*) '    * gerando momentos'
  momentos = gerar_vetores3d(N, int_momentos(1), int_momentos(2))

END SUBROUTINE gerarValores

! ************************************************************
!! Condicionamento: momento angular
!
! Objetivos:
!   Condiciona o momento angular para assumir determinado vetor
!   desejado.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE condicionar_momentoAngular (J, massas, posicoes, momentos)
  
  IMPLICIT NONE
  REAL(pf), INTENT(INOUT) :: posicoes(:,:), momentos(:,:), massas(:)
  INTEGER                 :: a
  REAL(pf)                :: momentoAngular_total(3), vetorRotacao(3), tensorInercia(3,3), J(3)

  ! Calcula o momento angular total
  momentoAngular_total = momento_angular_total (posicoes,momentos)

  ! Calcular o tensor de inercia
  tensorInercia = tensor_inercia_geral(massas, posicoes)

  ! Calcula o vetor de rotacao (resolve sistema linear)
  vetorRotacao = sistema_linear3 (tensorInercia, - momentoAngular_total + J)

  ! Percorre os corpos
  DO a = 1, SIZE(posicoes,1)
    ! Produto vetorial da posicao pelo vetor de rotacao
    momentos(a,:) = momentos(a,:) + massas(a) * produto_vetorial(posicoes(a,:), vetorRotacao)
  END DO

END SUBROUTINE condicionar_momentoAngular

! ************************************************************
!! Condicionamento: energia total
!
! Objetivos:
!   Condiciona a energia total para assumir determinado valor
!   desejado.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE condicionar_energiaTotal (H, G, massas, posicoes, momentos)
  ! Outra forma de H>0 seria somente aplicar o fator ((H-EP)/EC)**0.5

  IMPLICIT NONE
  REAL(pf), INTENT(IN)    :: H, G
  REAL(pf), INTENT(INOUT) :: posicoes(:,:), momentos(:,:), massas(:)
  REAL(pf)                :: EP, EC, fator

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
  ENDIF

END SUBROUTINE condicionar_energiaTotal

! ************************************************************
!! Condicionamento: centro de massas (anula)
!
! Objetivos:
!   Condiciona o centro de massas para assumir determinado vetor
!   nulo (origem).
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE zerar_centroMassas (massas, posicoes)

  IMPLICIT NONE
  REAL(pf), INTENT(INOUT) :: massas(:), posicoes(:,:)
  REAL(pf), DIMENSION(3)  :: rcm
  INTEGER                 :: a

  rcm = centro_massas(massas, posicoes)
  DO a = 1, SIZE(massas)
    posicoes(a,:) = posicoes(a,:) - rcm
  END DO

END SUBROUTINE zerar_centroMassas

! ************************************************************
!! Condicionamento: momento linear
!
! Objetivos:
!   Condiciona o momento linear para assumir determinado vetor
!   desejado.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE condicionar_momentoLinear (P, massas, momentos)

  IMPLICIT NONE
  REAL(pf), INTENT(INOUT) :: massas(:), momentos(:,:)
  REAL(pf)                :: P(3)
  REAL(pf), DIMENSION(3)  :: pcm ! analogo ao rcm
  INTEGER                 :: a

  ! Usa o mesmo metodo porque a ideia eh exatamente igual
  pcm = (momentoLinear_total(momentos) - P)/ SUM(massas)
  ! Substitui
  DO a = 1, SIZE(massas)
    momentos(a,:) = momentos(a,:) - massas(a)*pcm
  END DO

END SUBROUTINE condicionar_momentoLinear

! ************************************************************
!! Condicionamento de integrais primeiras
!
! Objetivos:
!   Condiciona vetores ja existentes com valores desejados.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE condicionar_ip (G, massas, posicoes, momentos, H, J, P)
  
  IMPLICIT NONE
  REAL(pf), INTENT(INOUT) :: G, posicoes(:,:), momentos(:,:), massas(:)
  REAL(pf) :: H, J(3), P(3)
  REAL(pf) :: erro_0, erro_1, ENERGIA, LINEAR(3), ANGULAR(3) ! Para medir a taxa de erro
  INTEGER :: i = 0

  ! Zera o centro de massas
  CALL zerar_centroMassas(massas, posicoes)
      
  ! Condiciona o momento linear
  CALL condicionar_momentoLinear(P, massas, momentos)

  ! Condiciona o momento angular
  CALL condicionar_momentoAngular(J, massas, posicoes, momentos)

  ! ! Condiciona a energia total
  CALL condicionar_energiaTotal(H, G, massas, posicoes, momentos)

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
      CALL condicionar_momentoLinear(P, massas, momentos)

      ! Condiciona o momento angular
      CALL condicionar_momentoAngular(J, massas, posicoes, momentos)

      ! ! Condiciona a energia total
      CALL condicionar_energiaTotal(H, G, massas, posicoes, momentos)

      ! Calculo dos erros
      ENERGIA = energia_total(G,massas,posicoes,momentos) - H
      LINEAR = momentoLinear_total(momentos) - P
      ANGULAR = momento_angular_total(posicoes,momentos) - J
      erro_1 = NORM2((/ENERGIA,LINEAR(1),LINEAR(2),LINEAR(3),ANGULAR(1),ANGULAR(2),ANGULAR(3)/))
    END DO
  ENDIF
  WRITE (*,*) ' > condicionamento aplicado ', i, ' vezes para obter o erro ', erro_1
END SUBROUTINE condicionar_ip

! ************************************************************
!! Gera valores condicionados pelas integrais primeiras
!
! Objetivos:
!   Gera vetores aleatorios e os condicionado para valores
!   informados.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE gerar_condicionado_ip (G, N, massas, posicoes, momentos, int_posicoes, int_momentos, int_massas, H, J, P)

  IMPLICIT NONE
  INTEGER, INTENT(IN)     :: N
  REAL(pf), DIMENSION(2), INTENT(IN) :: int_posicoes, int_momentos, int_massas
  REAL(pf), INTENT(INOUT) :: G, H, J(3), P(3)
  REAL(pf), INTENT(INOUT), allocatable :: massas(:), posicoes(:,:), momentos(:,:)

  WRITE (*,'(a)') "GERACAO DAS CONDICOES INICIAIS"
  
  ! Gera os valores
  WRITE (*,'(a)') '  > gerando valores...'
  ALLOCATE(massas (N))
  ALLOCATE(posicoes (N, 3))
  ALLOCATE(momentos (N, 3))
  CALL gerarValores(N, massas, posicoes, momentos, int_posicoes, int_momentos, int_massas)

  ! Condiciona
  WRITE (*,'(a)') '  > condicionando...'
  CALL condicionar_ip(G, massas, posicoes, momentos, H, J, P)

  ! Exibe as integrais primeiras do sistema
  WRITE (*,*) '    * H   =', energia_total(G,massas,posicoes,momentos) 
  WRITE (*,*) '    * Rcm =', centro_massas(massas,posicoes) 
  WRITE (*,*) '    * P   =', momentoLinear_total(momentos) 
  WRITE (*,*) '    * J   =', momento_angular_total(posicoes,momentos) 

  WRITE (*,'(a)') '  > condicoes iniciais geradas!'
  WRITE (*,*)

END SUBROUTINE gerar_condicionado_ip


! ************************************************************
!! Gera valores condicionados pelo modelo de Henon
!
! Objetivos:
!   Gera vetores aleatorios e os condicionado para valores
!   informados.
!
! Modificado:
!   12 de agosto de 2024
!
! Autoria:
!   oap
! 
! DE FATO, NAO USAREI G, H, J E P
SUBROUTINE gerar_condicionado_henon (N, massas, posicoes, momentos, int_posicoes, int_momentos, int_massas)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  REAL(pf), DIMENSION(2), INTENT(IN) :: int_posicoes, int_momentos, int_massas
  REAL(pf), INTENT(INOUT), allocatable :: massas(:), posicoes(:,:), momentos(:,:)
  REAL(pf) :: beta, EP, EC, Qvir ! AARSETH
  integer :: a, m

  WRITE (*,'(a)') "GERACAO DAS CONDICOES INICIAIS (HENON)"
  
  ! Gera os valores
  WRITE (*,'(a)') '  > gerando valores...'
  ALLOCATE(massas (N))
  ALLOCATE(posicoes (N, 3))
  ALLOCATE(momentos (N, 3))
  CALL gerarValores(N, massas, posicoes, momentos, int_posicoes, int_momentos, int_massas)
  
  m = SUM(massas)
  do a = 1, N
    massas(a) = (1.0_pf / N)
  end do
  WRITE (*,*) "M:", SUM(massas)

  ! Condiciona
  WRITE (*,'(a)') '  > condicionando... (HENON)'
  
  ! Zera o centro de massas
  CALL zerar_centroMassas(massas, posicoes)
  
  ! AQUI EU CONSIGO QUE H = -0.25 E V ~ -0.25
  CALL condicionar_energiaTotal(-0.25_pf, 1.0_pf, massas, posicoes, momentos)
    
  ! AGORA EU QUERO QUE V ~ -0.5, OU SEJA, V = 2 * V
  posicoes = posicoes * 0.5_pf
  EP = energia_potencial(1.0_pf, massas, posicoes)
  WRITE (*,*) '    * V   =', EP

  ! AGORA QUERO QUE T/(-V) = 0.5, LOGO v = sqrt(-V/M)
  momentos = SQRT(0.5 / 3) / N

  EC = energia_cinetica(massas, momentos)
  WRITE (*,*) '    * T   =', EC
  WRITE (*,*) '    * H   =', energia_total(1.0_pf,massas,posicoes,momentos) 
  WRITE (*,*) '    * Q   =', EC/ABS(EP)
  WRITE (*,*) '    * R   =', - 1.0_pf * SUM(massas)**2 / (2 * EP)

END SUBROUTINE

END MODULE condicoesIniciais