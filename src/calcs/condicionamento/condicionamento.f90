! *****************************************************************
!! CONDICIONAMENTO DE VALORES INICIAIS
!
! Objetivos:
!   Funcoes para a geracao e condicionamento a partir de restricoes.
! 
! Modificado:
!   05 de maio de 2025
! 
! Autoria:
!   oap
! 
MODULE condicionamento
  USE tipos
  USE mecanica
  USE auxiliares
  USE aleatorio
  USE json_utils_mod
  IMPLICIT NONE
CONTAINS

! ************************************************************
!! Gera valores
!
! Objetivos:
!   Gera valores iniciais aleatorios.
!
! Modificado:
!   02 de maio de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE gerar_valores (N, sorteio, massas, posicoes, momentos, mi)
  TYPE(json_value), POINTER, INTENT(IN) :: sorteio
  INTEGER, INTENT(IN)         :: N
  TYPE(json_value), POINTER   :: sort_mas, sort_pos, sort_mom
  REAL(pf), INTENT(INOUT)     :: posicoes(:,:), momentos(:,:), massas(:)
  LOGICAL, INTENT(IN)         :: mi ! massas iguais

  CALL json % get(sorteio, "massas", sort_mas)
  CALL json % get(sorteio, "posicoes", sort_pos)
  CALL json % get(sorteio, "momentos", sort_mom)

  ! Gera massas
  IF (mi) THEN
    WRITE (*,'(A)') '    * gerando massas (1/N)'
    massas = 1/N
  ELSE
    WRITE (*,'(A)', ADVANCE='no') '    * gerando massas'
    massas = gerar_massas(N, sort_mas)
  ENDIF

  ! Gera as posições
  WRITE (*,'(A)', ADVANCE='no') '    * gerando posicoes'  
  posicoes = gerar_vetores3d(N, sort_pos)

  ! Gera os momentos
  WRITE (*,'(A)', ADVANCE='no') '    * gerando momentos'
  momentos = gerar_vetores3d(N, sort_mom)

END SUBROUTINE gerar_valores

! ************************************************************
!! Gera vetores 3d
!
! Objetivos:
!   Gera vetores 3d utilizados para as posicoes e velocidades.
!
! Modificado:
!   02 de maio de 2025
!
! Autoria:
!   oap
! 
FUNCTION gerar_vetores3d (N, sorteio) RESULT(vetores)

  INTEGER, INTENT(IN)         :: N
  TYPE(json_value), POINTER, INTENT(INOUT) :: sorteio
  REAL(pf), DIMENSION(N,3)    :: vetores
  CHARACTER(:), ALLOCATABLE   :: distribuicao, regiao
  REAL(pf), ALLOCATABLE     :: intervalo(:)
  REAL(pf)                  :: raio, vmin, vmax, distmin

  distribuicao = json_get_string(sorteio, "distribuicao")
  regiao = json_get_string(sorteio, "regiao")
  raio = json_get_float(sorteio, "raio")

  WRITE (*,'(A)') ' ('//distribuicao//')'
  
  ! Intervalo de sorteio
  intervalo = json_get_float_vec(sorteio, "intervalo")
  vmin = intervalo(1)
  vmax = intervalo(2)
  distmin = 0.0_pf
  ! Distancia minima
  IF (SIZE(intervalo) == 3) THEN
    distmin = intervalo(3)
  ENDIF

  SELECT CASE (TRIM(distribuicao))
    ! Uniforme (0,1)
    CASE ("uniforme"); CALL uniforme(vetores, N, distmin, vmin, vmax, regiao, raio)
    ! Normal (0,1)
    CASE ("normal"); CALL normal(vetores, N, distmin, regiao, raio)
    ! Cauchy
    CASE ("cauchy"); CALL cauchy(vetores, N, distmin, regiao, raio)
  END SELECT

END FUNCTION gerar_vetores3d

! ************************************************************
!! Gera vetor de massas
!
! Objetivos:
!   Gera vetor de massas conforme um intervalo informado.
!
! Modificado:
!   02 de maio de 2025
!
! Autoria:
!   oap
! 
FUNCTION gerar_massas (N, sorteio) RESULT(massas)

  INTEGER, INTENT(IN)           :: N
  TYPE(json_value), POINTER     :: sorteio
  REAL(pf), DIMENSION(N)        :: massas
  REAL(pf), ALLOCATABLE         :: intervalo(:)
  CHARACTER(LEN=:), ALLOCATABLE :: distribuicao
  REAL(pf)               :: vet_min(N)

  intervalo = json_get_float_vec(sorteio, "intervalo")
  vet_min = intervalo(1)

  ! AQUI PRECISA APLICAR A DISTRIBUICAO DESEJADA. POR ENQUANTO SO TEM A UNIFORME
  distribuicao = json_get_string(sorteio, "distribuicao")
  WRITE (*,'(A)') ' ('//distribuicao//')'

  IF (distribuicao == "uniforme") THEN
    CALL RANDOM_SEED()
    CALL RANDOM_NUMBER(massas)
    
    ! Agora condiciona no intervalo
    massas = massas * (intervalo(2) - intervalo(1) + 1) + vet_min
  ENDIF

END FUNCTION gerar_massas

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
!! Condicionamento iterativo de integrais primeiras
!
! Objetivos:
!   Condiciona vetores ja existentes com valores desejados de
!   maneira iterativa.
!
! Modificado:
!   02 de maio de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE condicionar_ip_iterativo (G, massas, posicoes, momentos, H, J, P)
  REAL(pf), INTENT(IN) :: G, H, J(3), P(3)
  REAL(pf), INTENT(INOUT) :: massas(:), posicoes(:,:), momentos(:,:)
  REAL(pf) :: erro_0, erro_1, erro_limite = 0.1E-8
  REAL(pf) :: er_energia, er_linear, er_angular
  INTEGER :: N_iter_max = 10, i

  ! Zera o centro de massas
  CALL zerar_centroMassas(massas, posicoes)

  ! Condiciona o momento linear
  CALL condicionar_momentoLinear(P, massas, momentos)

  ! Condiciona o momento angular
  CALL condicionar_momentoAngular(J, massas, posicoes, momentos)

  ! Condiciona a energia total
  CALL condicionar_energiaTotal(H, G, massas, posicoes, momentos)

  ! Calculo dos erros
  er_energia = ABS(energia_total(G,massas, posicoes, momentos) - H)
  er_linear = MAXVAL(ABS(momentoLinear_total(momentos) - P))
  er_angular = MAXVAL(ABS(momento_angular_total(posicoes,momentos) - J))
  erro_1 = MAXVAL((/er_energia,er_linear,er_angular/))

  IF (erro_1 >= erro_limite) THEN
    erro_0 = ABS(erro_1) + 1.0 ! Para entrar no loop
    i = 0
    DO WHILE (ABS(erro_1) >= erro_limite .AND. i <= N_iter_max)
      i = i + 1
      erro_0 = erro_1

      ! Condiciona o momento linear
      CALL condicionar_momentoLinear(P, massas, momentos)

      ! Condiciona o momento angular
      CALL condicionar_momentoAngular(J, massas, posicoes, momentos)

      ! Condiciona a energia total
      CALL condicionar_energiaTotal(H, G, massas, posicoes, momentos)

      ! Calculo dos erros
      er_energia = ABS(energia_total(G,massas, posicoes, momentos) - H)
      er_linear = MAXVAL(ABS(momentoLinear_total(momentos) - P))
      er_angular = MAXVAL(ABS(momento_angular_total(posicoes,momentos) - J))
      erro_1 = MAXVAL((/er_energia,er_linear,er_angular/))
    END DO
  END IF

  WRITE (*,*) ' > condicionamento iterativo aplicado ', i, ' vezes para obter o erro ', erro_1
END SUBROUTINE condicionar_ip_iterativo

! ************************************************************
!! Condicionamento direto de integrais primeiras
!
! Objetivos:
!   Condiciona vetores ja existentes com valores desejados de
!   maneira direta.
!
! Modificado:
!   05 de maio de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE condicionar_ip_direto (G, massas, posicoes, momentos, H, J, P)
  REAL(pf), INTENT(IN) :: G
  REAL(pf), INTENT(IN), OPTIONAL :: H, J(3), P(3)
  REAL(pf) :: ed, pd(3), jd(3)
  REAL(pf), INTENT(INOUT) :: massas(:), posicoes(:,:), momentos(:,:)
  REAL(pf) :: energia, linear(3), angular(3), M_tot_inv
  REAL(pf) :: alpha, potencial ! calculo do alpha
  REAL(pf) :: beta, inercia(3,3), rot(3), rot_(3), S1, S2, K1a(3), K2a(3) ! calculo do beta
  REAL(pf), DIMENSION(:,:), ALLOCATABLE :: n_posicoes, n_momentos
  INTEGER :: i

  ed = 0.0_pf
  pd = 0.0_pf
  jd = 0.0_pf
  IF (PRESENT(H)) ed = H
  IF (PRESENT(J)) jd = J
  IF (PRESENT(P)) pd = P

  ! Zera o centro de massas
  CALL zerar_centroMassas(massas, posicoes)

  ! Calcula as integrais primeiras
  energia = energia_total(G, massas, posicoes, momentos)
  linear = momentoLinear_total(momentos)
  angular = momento_angular_total(posicoes, momentos)
  M_tot_inv = 1.0_pf / SUM(massas)

  ! Calculo do alpha
  potencial = energia_potencial(G, massas, posicoes)
  alpha = 1.0_pf + ed / potencial

  ! Calculo do beta
  inercia = tensor_inercia_geral(massas, posicoes)
  rot  = sistema_linear3(inercia, angular)
  rot_ = sistema_linear3(inercia, jd) * alpha

  S1 = 0.0_pf
  S2 = 0.0_pf

  DO i = 1, SIZE(massas)
    K1a = momentos(i,:) - massas(i) * linear * M_tot_inv
    K1a = K1a - massas(i) * produto_vetorial(posicoes(i,:), rot)
    S1 = S1 + 0.5_pf * DOT_PRODUCT(K1a, K1a) / massas(i)

    K2a = pd * M_tot_inv + produto_vetorial(posicoes(i,:), rot_)
    S2 = S2 + 0.5_pf * DOT_PRODUCT(K2a, K2a) * massas(i)
  END DO

  beta = SQRT((-potencial - S2)/S1)

  WRITE (*,*) '   > coeficientes:'
  WRITE (*,*) '     * alpha =', alpha
  WRITE (*,*) '     * beta  =', beta
  WRITE (*,*) '     * S1    =', S1
  WRITE (*,*) '     * S2    =', S2

  IF (alpha == 0.0 .OR. beta == 0.0) THEN
    ERROR STOP "ERRO: um dos coeficientes alpha ou beta eh nulo."
  END IF

  ! Transforma as coordenadas
  ALLOCATE(n_posicoes(SIZE(massas), 3))
  ALLOCATE(n_momentos(SIZE(massas), 3))

  n_posicoes = posicoes / alpha

  rot = sistema_linear3(inercia, angular - jd*alpha/beta)

  DO i = 1, SIZE(massas)
    n_momentos(i,:) = momentos(i,:) - massas(i) * M_tot_inv * (linear - pd / beta)
    n_momentos(i,:) = n_momentos(i,:) - massas(i) * produto_vetorial(posicoes(i,:), rot)
    n_momentos(i,:) = beta * n_momentos(i,:)
  END DO

  ! Salva
  posicoes = n_posicoes
  momentos = n_momentos

END SUBROUTINE condicionar_ip_direto

! ************************************************************
!! Condicionamento de Aarseth
!
! Objetivos:
!   Condiciona vetores ja existentes atraves do proposto 
!   por (AARSETH, 2003) com as unidades padrao: massa total
!   unitaria, energia total -0.25, relacao do virial atendida.
!
! Modificado:
!   05 de maio de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE condicionar_aarseth (G, massas, posicoes, momentos)
  REAL(pf), INTENT(IN)    :: G
  REAL(pf), INTENT(INOUT) :: massas(:), posicoes(:,:), momentos(:,:)
  REAL(pf) :: energia = -0.25_pf  ! energia padrao
  REAL(pf) :: virial  = 0.5_pf    ! relacao de virial padrao
  REAL(pf) :: ec, ep              ! cinetica e potencial
  REAL(pf) :: Qv, beta

  ! Normaliza as massas
  massas(:) = 1.0_pf/SIZE(massas)

  ! Para comecar, condiciona as integrais primeiras para zero
  CALL condicionar_ip_direto(G, massas, posicoes, momentos)

  ! Metodo de Aarseth
  ep = energia_potencial(G, massas, posicoes)
  ec = energia_cinetica(massas, momentos)

  Qv = SQRT(virial * ABS(ep) / ec)
  beta = (1.0_pf - virial) * ep / energia

  posicoes = posicoes * beta
  momentos = momentos * Qv / SQRT(beta)
END SUBROUTINE condicionar_aarseth

END MODULE condicionamento