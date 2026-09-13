! *****************************************************************
!! CONDICIONAMENTO DE VALORES INICIAIS
!
! Objetivos:
!   Funcoes para a geracao e condicionamento a partir de restricoes.
! 
! Modificado:
!   16 de dezembro de 2025
! 
! Autoria:
!   oap
! 
MODULE condicionamento
  USE iso_fortran_env, only: output_unit
  USE tipos
  USE utilidades
  USE sorteio_mod
  IMPLICIT NONE
CONTAINS

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
!   03 de agosto de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE condicionar_energiaTotal (H, G, massas, posicoes, momentos, eps)
  ! Outra forma de H>0 seria somente aplicar o fator ((H-EP)/EC)**0.5

  IMPLICIT NONE
  REAL(pf), INTENT(IN)    :: H, G, eps
  REAL(pf), INTENT(INOUT) :: posicoes(:,:), momentos(:,:), massas(:)
  REAL(pf)                :: EP, EC, fator
  
  REAL(pf), ALLOCATABLE :: dists2(:), massas2(:), posicoes_fator(:,:)
  REAL(pf) :: erro_newton, pot_fator, pot_der_fator, const, Rab(3)
  INTEGER  :: a, b, i, contador
  INTEGER :: N

  N = SIZE(massas)
  
  ! Calcula as energias
  EP = energia_potencial(G, massas, posicoes, eps)
  EC = energia_cinetica(massas, momentos)

  ! Calcula o fator para bater a energia
  IF (H > 0) THEN
    fator = ((-EP+H)/EC)**0.5
  ELSE
    fator = (-EP/EC)**0.5
  ENDIF

  ! Aplica sobre os momentos
  momentos = fator * momentos

  ! Se a energia desejada nao for positiva, precisa mexer nas posicoes
  IF (H < 0) THEN
    ! Se nao tiver amortecimento, pode aplicar homotetia pois o potencial eh homogeneo
    IF (eps == 0) THEN
      fator = 1.0_pf/(H/EP + 1.0_pf)
      posicoes = fator * posicoes
    
    ! Se tiver amortecimento, precisa usar Newton para encontrar o fator
    ELSE
      CALL cond_pot_amortecido_energia(H, G, massas, posicoes, momentos, eps)
    ENDIF
  ENDIF

END SUBROUTINE condicionar_energiaTotal

! ************************************************************
!! Condicionamento do potencial amortecido [ENERGIA]
!
! Objetivos:
!   Condiciona o potencial amortecido (softening) utilizando 
!   o metodo de Newton, necessario devido a perda de
!   homogeneidade.
! 
! Atencao:
!   Esta rotina nao implementa equilibrio, apenas redimensiona
!   o sistema para obter um V = H - T
!
! Modificado:
!   25 de julho de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE cond_pot_amortecido_energia (H, G, massas, posicoes, momentos, eps, fator_inicial)
  REAL(pf), INTENT(IN) :: H, G, eps, massas(:), momentos(:,:)
  REAL(pf), INTENT(IN), OPTIONAL :: fator_inicial
  REAL(pf), INTENT(INOUT) :: posicoes(:,:)
  REAL(pf), ALLOCATABLE :: dists2(:), massas2(:), posicoes_fator(:,:)
  REAL(pf) :: erro_newton, pot_fator, pot_der_fator, const, Rab(3), fator, EP, EC
  INTEGER :: a, b, i, contador
  INTEGER :: N

  ! Definicoes
  N = SIZE(massas)
  ALLOCATE(dists2(INT(N*(N-1)/2)))
  ALLOCATE(massas2(INT(N*(N-1)/2)))
  ALLOCATE(posicoes_fator(N,3))

  ! Vetores para facilitar o calculo do potencial e sua derivada
  i = 1
  EP = 0.0_pf
  DO a = 2, N
    DO b = 1, a - 1
      Rab = posicoes(a,:) - posicoes(b,:)
      dists2(i) = DOT_PRODUCT(Rab, Rab)
      massas2(i) = massas(a) * massas(b)
      EP = EP - G * massas2(i) / SQRT(dists2(i) + eps*eps)
      i = i + 1
    END DO
  END DO
  EC = energia_cinetica(massas, momentos)

  ! Constante
  const = -eps * (H - EC)

  ! Agora vamos aplicar o metodo de Newton
  erro_newton = 1.0
  IF (PRESENT(fator_inicial)) THEN
    fator = fator_inicial/eps
  ELSE
    fator = 1.0_pf/(eps * (H/EP + 1.0_pf)) ! Chute inicial
  ENDIF
  contador = 0

  DO WHILE (erro_newton > 1E-15 .AND. contador < 50)
    ! Calcula o potencial e sua derivada em relacao ao fator
    pot_fator = 0.0_pf
    pot_der_fator = 0.0_pf

    DO i = 1, SIZE(dists2)
      pot_fator = pot_fator - G*massas2(i)/SQRT(fator*fator*dists2(i) + 1)
      pot_der_fator = pot_der_fator + G*massas2(i)*dists2(i)/SQRT(fator*fator*dists2(i) + 1)**3
    END DO
    pot_der_fator = pot_der_fator * fator

    ! Agora aplica o metodo de Newton
    fator = fator - (pot_fator + const)/pot_der_fator

    ! Vamos verificar o erro
    posicoes_fator = fator * eps * posicoes
    erro_newton = ABS(H - energia_total(G, massas, posicoes_fator, momentos, eps))

    contador = contador + 1
  END DO

  ! Atualiza as posicoes
  posicoes = posicoes_fator

  ! Liberando memoria
  DEALLOCATE(dists2)
  DEALLOCATE(massas2)
  DEALLOCATE(posicoes_fator)
END SUBROUTINE

! ************************************************************
!! Condicionamento do potencial amortecido [EQUILIBRIO]
!
! Objetivos:
!   Condiciona o potencial amortecido (softening) utilizando 
!   o metodo de Newton, necessario devido a perda de
!   homogeneidade.
! 
! Atencao:
!   Esta rotina implementa equilibrio, redimensionando o
!   sistema de modo que 2T + <F,q> = 0. Para condicionar a
!   a energia total, utilize `cond_pot_amortecido_energia`.
!
! Modificado:
!   25 de julho de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE cond_pot_amortecido_equilibrio (H, G, massas, posicoes, eps, fator_inicial)
  REAL(pf), INTENT(IN) :: H, G, eps, massas(:)
  REAL(pf), INTENT(IN), OPTIONAL :: fator_inicial
  REAL(pf), INTENT(INOUT) :: posicoes(:,:)
  REAL(pf), ALLOCATABLE :: posicoes_fator(:,:)
  REAL(pf) :: erro_newton, pot_fator, pot_der_fator, const, Rab(3), fator
  INTEGER :: a, b, i, contador
  INTEGER :: N

  REAL(pf), ALLOCATABLE :: forcas(:,:), forcas_der(:,:)
  REAL(pf) :: pot, pot_der, mult
  REAL(pf) :: Fab(3), dist, f, f_der, denominador

  ! Definicoes
  N = SIZE(massas)
  ALLOCATE(forcas(N,3))
  ALLOCATE(forcas_der(N,3))
  ALLOCATE(posicoes_fator(N,3))

  ! Constante
  const = 2 * H

  ! Agora vamos aplicar o metodo de Newton
  erro_newton = 1.0_pf
  IF (PRESENT(fator_inicial)) THEN
    fator = fator_inicial
  ELSE
    fator = 1.0_pf ! Chute inicial
  ENDIF
  contador = 0

  DO WHILE (erro_newton > 1E-15 .AND. contador < 10)
    forcas = 0.0_pf
    forcas_der = 0.0_pf
    pot = 0.0_pf
    pot_der = 0.0_pf
    DO a = 2, N
      DO b = 1, a - 1
        Rab =  posicoes(b,:) - posicoes(a,:)
        dist = DOT_PRODUCT(Rab, Rab)
        denominador = SQRT(dist + (eps/fator)**2)
        mult = G * massas(a) * massas(b)

        ! Forca
        Fab = mult * Rab / denominador**3
        forcas(a,:) = forcas(a,:) + Fab
        forcas(b,:) = forcas(b,:) - Fab

        ! Derivada da forca
        Fab = Fab / (denominador*denominador)
        forcas_der(a,:) = forcas_der(a,:) + Fab
        forcas_der(b,:) = forcas_der(b,:) - Fab

        ! Potencial
        pot = pot - mult / denominador

        ! Derivada do potencial
        pot_der = pot_der - mult / (denominador**3)
      END DO
    END DO

    f = const - 2 * pot/fator
    f_der = 2 * pot/(fator*fator) - 2*pot_der*eps*eps/(fator**4)
    DO a = 1, N
      ! Funcao
      f = f + DOT_PRODUCT(forcas(a,:), posicoes(a,:)) / fator

      ! Derivada
      f_der = f_der - DOT_PRODUCT(forcas(a,:), posicoes(a,:)) / (fator*fator)
      f_der = f_der + DOT_PRODUCT(forcas_der(a,:), posicoes(a,:)) * 3*eps*eps/(fator**4)
    END DO

    ! Agora aplica o metodo de Newton
    fator = fator - f/f_der

    ! Vamos verificar o erro
    posicoes_fator = fator * posicoes
    forcas = 0.0_pf
    pot = 0.0_pf
    DO a = 2, N
      DO b = 1, a - 1
        Rab = posicoes_fator(b,:) - posicoes_fator(a,:)
        dist = DOT_PRODUCT(Rab, Rab)
        denominador = SQRT(dist + eps**2)
        mult = G * massas(a) * massas(b)

        Fab = mult * Rab / denominador**3
        forcas(a,:) = forcas(a,:) + Fab
        forcas(b,:) = forcas(b,:) - Fab

        pot = pot - mult / denominador
      END DO
    END DO

    erro_newton = const - 2*pot
    DO a = 1, N
      erro_newton = erro_newton + DOT_PRODUCT(forcas(a,:),posicoes_fator(a,:))
    END DO
    erro_newton = ABS(erro_newton)

    contador = contador + 1
  END DO

  ! Atualiza as posicoes
  posicoes = posicoes_fator

  ! Liberando memoria
  DEALLOCATE(forcas)
  DEALLOCATE(forcas_der)
  DEALLOCATE(posicoes_fator)
END SUBROUTINE

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
  pcm = (momento_linear_total(momentos) - P)/ SUM(massas)
  ! Substitui
  DO a = 1, SIZE(massas)
    momentos(a,:) = momentos(a,:) - massas(a)*pcm
  END DO

END SUBROUTINE condicionar_momentoLinear


! ************************************************************
!! Verifica condicoes de Sundman e Delta
!
! Modificado:
!   03 de agosto de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE verificar_condicoes_existencia (G, massas, posicoes, momentos, eps, ed, jd, pd)
  REAL(pf), INTENT(IN) :: G
  REAL(pf), INTENT(IN) :: ed, jd(3), pd(3), eps
  REAL(pf), INTENT(INOUT) :: massas(:), posicoes(:,:), momentos(:,:)
  REAL(pf) :: M_tot_inv, potencial
  REAL(pf) :: rot_(3), sigma_til, inercia(3,3), delta
  REAL(pf) :: momine, sundman

  ! Se nao tiver momentos, nao precisa verificar nada
  IF (NORM2(jd) == 0 .AND. NORM2(pd) == 0) RETURN

  ! Se tiver momentos
  ! Valores para calcular as condicoes
  M_tot_inv = 1.0_pf / SUM(massas)
  potencial = energia_potencial(G, massas, posicoes, eps)
  inercia = tensor_inercia_geral(massas, posicoes)
  rot_ = sistema_linear3(inercia, jd)
  sigma_til = - DOT_PRODUCT(jd, rot_)
  momine = momento_inercia(massas, posicoes)

  ! Verifica condicao de Delta
  delta = potencial**2 - sigma_til * (M_tot_inv * NORM2(pd)**2 - 2.0_pf * ed)  
  IF (delta < 0) THEN
    PRINT *, ' [ATENCAO] Criterio do Delta nao atingido, energia total insuficiente.'
    STOP 0
  ENDIF

  ! Verifica desigualdade de Sundman
  sundman = momine * potencial**2 + 2 * NORM2(jd)**2 * ed
  IF (sundman < 0) THEN
    PRINT *, ' [ATENCAO] Desigualdade de Sundman nao atendida, energia total insuficiente.'
    STOP 0
  ENDIF

END SUBROUTINE verificar_condicoes_existencia

! ************************************************************
!! Condicionamento iterativo de integrais primeiras
!
! Objetivos:
!   Condiciona vetores ja existentes com valores desejados de
!   maneira iterativa.
!
! Modificado:
!   03 de agosto de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE condicionar_ip_iterativo (G, massas, posicoes, momentos, eps, H, mat, mlt, nim)
  REAL(pf), INTENT(IN) :: G, H, eps
  REAL(pf), INTENT(INOUT) :: massas(:), posicoes(:,:), momentos(:,:)
  REAL(pf), INTENT(IN), OPTIONAL :: mat(3), mlt(3)
  REAL(pf) :: J(3), P(3)
  REAL(pf) :: erro_0, erro_1, erro_limite = 0.1E-8
  REAL(pf) :: er_energia, er_linear, er_angular
  INTEGER, INTENT(IN), OPTIONAL :: nim ! N_iter_max
  INTEGER :: N_iter_max = 10, i

  N_iter_max = 10
  IF (PRESENT(nim)) N_iter_max = nim

  J = 0.0_pf
  P = 0.0_pf
  IF (PRESENT(mat)) J = mat
  IF (PRESENT(mlt)) P = mlt

  ! Zera o centro de massas
  CALL zerar_centroMassas(massas, posicoes)

  ! Condiciona o momento linear
  CALL condicionar_momentoLinear(P, massas, momentos)

  ! Condiciona o momento angular
  CALL condicionar_momentoAngular(J, massas, posicoes, momentos)

  ! Condiciona a energia total
  CALL condicionar_energiaTotal(H, G, massas, posicoes, momentos, eps)

  ! Calculo dos erros
  er_energia = ABS(energia_total(G,massas, posicoes, momentos, eps) - H)
  er_linear = MAXVAL(ABS(momento_linear_total(momentos) - P))
  er_angular = MAXVAL(ABS(momento_angular_total(posicoes,momentos) - J))
  erro_1 = MAXVAL((/er_energia,er_linear,er_angular/))

  i = 1
  IF (erro_1 >= erro_limite) THEN
    DO WHILE (ABS(erro_1) >= erro_limite .AND. i < N_iter_max)
      i = i + 1

      ! Condiciona o momento linear
      CALL condicionar_momentoLinear(P, massas, momentos)

      ! Condiciona o momento angular
      CALL condicionar_momentoAngular(J, massas, posicoes, momentos)

      ! Condiciona a energia total
      CALL condicionar_energiaTotal(H, G, massas, posicoes, momentos, eps)

      ! Calculo dos erros
      er_energia = ABS(energia_total(G,massas, posicoes, momentos, eps) - H)
      er_linear = MAXVAL(ABS(momento_linear_total(momentos) - P))
      er_angular = MAXVAL(ABS(momento_angular_total(posicoes,momentos) - J))
      erro_1 = MAXVAL((/er_energia,er_linear,er_angular/))
    END DO
  END IF

  WRITE (output_unit,*)
  WRITE (output_unit,*) ' > condicionamento iterativo aplicado ', i, ' vezes para obter o erro ', erro_1
  WRITE (output_unit,*)
END SUBROUTINE condicionar_ip_iterativo

! ************************************************************
!! Condicionamento direto de integrais primeiras
!
! Objetivos:
!   Condiciona vetores ja existentes com valores desejados de
!   maneira direta.
!
! Modificado:
!   03 de agosto de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE condicionar_ip_direto (G, massas, posicoes, momentos, soft, H, J, P)
  REAL(pf), INTENT(IN) :: G
  REAL(pf), INTENT(IN), OPTIONAL :: H, J(3), P(3), soft
  REAL(pf) :: ed, pd(3), jd(3), eps
  REAL(pf), INTENT(INOUT) :: massas(:), posicoes(:,:), momentos(:,:)
  REAL(pf) :: energia, linear(3), angular(3), M_tot_inv
  REAL(pf) :: alpha, potencial, sigma_til ! calculo do alpha
  REAL(pf) :: beta, inercia(3,3), rot(3), rot_(3), S1, S2 ! calculo do beta
  REAL(pf), DIMENSION(:,:), ALLOCATABLE :: n_posicoes, n_momentos
  INTEGER :: i, sinal

  ed = 0.0_pf
  jd = 0.0_pf
  pd = 0.0_pf
  eps = 0.0_pf
  IF (PRESENT(H)) ed = H
  IF (PRESENT(J)) jd = J
  IF (PRESENT(P)) pd = P
  IF (PRESENT(soft)) eps = soft

  IF (eps .NE. 0) THEN
    WRITE (output_unit,*) ' [ATENCAO] potencial amortecido, o cond. direto nao sera exato!'
  END IF

  ! Zera o centro de massas
  CALL zerar_centroMassas(massas, posicoes)

  ! Verifica as condicoes de existencia
  CALL verificar_condicoes_existencia(G, massas, posicoes, momentos, eps, ed, jd, pd)

  ! Calcula as integrais primeiras e outros valores
  energia = energia_total(G, massas, posicoes, momentos, eps)
  linear = momento_linear_total(momentos)
  angular = momento_angular_total(posicoes, momentos)
  M_tot_inv = 1.0_pf / SUM(massas)
  potencial = energia_potencial(G, massas, posicoes, eps)

  ! Calculo do beta
  inercia = tensor_inercia_geral(massas, posicoes)
  rot  = sistema_linear3(inercia, angular)
  rot_ = sistema_linear3(inercia, jd)

  ! Se nao tiver momentos, ja deixa calculado
  ! Para E >= 0, basta mexer nas velocidades
  ! Para E < 0, eh necessario um incremento para existir beta
  alpha = 1.0_pf
  IF (ed < 0) alpha = alpha + ed / potencial

  ! Se tiver angular (mas nao linear), precisa adicionar
  ! o elemento de rotacao
  IF (NORM2(jd) > 0) THEN
    sigma_til = - DOT_PRODUCT(jd, rot_)
    alpha = - potencial / sigma_til
  
  ! Se nao tiver momento angular mas tiver linear, precisa
  ! adicionar a contribuicao linear para existir beta
  ELSE IF (NORM2(pd) > 0) THEN
    alpha = alpha - 0.5_pf * M_tot_inv * NORM2(pd)**2 / potencial
  ENDIF

  ! Calculando o beta
  S1 = (energia-potencial)-0.5_pf*M_tot_inv*NORM2(linear)**2+0.5_pf*DOT_PRODUCT(angular, rot)
  S2 = (0.5_pf * M_tot_inv * NORM2(pd)**2 - 0.5_pf*alpha*alpha*DOT_PRODUCT(jd, rot_))
  beta = SQRT((ed - alpha * potencial - S2)/S1)

  WRITE (output_unit,*) '   > coeficientes:'
  WRITE (output_unit,*) '     * alpha =', alpha
  WRITE (output_unit,*) '     * beta  =', beta
  WRITE (output_unit,*) '     * S1    =', S1
  WRITE (output_unit,*) '     * S2    =', S2

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
!   25 de julho de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE condicionar_aarseth (G, massas, posicoes, momentos, eps)
  REAL(pf), INTENT(IN)    :: G, eps
  REAL(pf), INTENT(INOUT) :: massas(:), posicoes(:,:), momentos(:,:)
  REAL(pf) :: energia = -0.25_pf  ! energia padrao
  REAL(pf) :: virial  = 0.5_pf    ! relacao de virial padrao
  REAL(pf) :: ec, ep              ! cinetica e potencial
  REAL(pf) :: Qv, beta

  ! Aviso para caso haja amortecedor
  IF (eps .NE. 0) THEN
    WRITE (output_unit,*) ' [ATENCAO] potencial amortecido, nao havera equilibrio inicial!'
    WRITE (output_unit,*) '           para equilibrio, use o "sorteio_aarseth_modificado".'
  END IF
  
  ! Normaliza as massas
  massas(:) = 1.0_pf/SIZE(massas)
  
  ! Para comecar, condiciona as integrais primeiras para zero
  IF (eps == 0) THEN
    CALL condicionar_ip_direto(G, massas, posicoes, momentos)
  ELSE
    CALL condicionar_ip_iterativo(G, massas, posicoes, momentos, eps, energia)
  ENDIF

  ! Metodo de Aarseth
  ep = energia_potencial(G, massas, posicoes, eps)
  ec = energia_cinetica(massas, momentos)

  Qv = SQRT(virial * ABS(ep) / ec)
  beta = (1.0_pf - virial) * ep / energia

  ! Condiciona as velocidades
  momentos = momentos * Qv / SQRT(beta)

  ! Se nao tiver amortecimento, o metodo eh direto
  IF (eps == 0) THEN
    posicoes = posicoes * beta
  ! Se tiver, precisa ser iterativo, usando beta/eps como chute inicial
  ELSE
    CALL cond_pot_amortecido_energia(energia, G, massas, posicoes, momentos, eps, beta)
  ENDIF
END SUBROUTINE condicionar_aarseth

! ************************************************************
!! Condicionamento de Aarseth Modificado
!
! Objetivos:
!   Igual ao condicionamento de Aarseth se o amortecedor for
!   nulo, mas condiciona devidamente se houver amortecimento
!
! Modificado:
!   25 de julho de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE condicionar_aarseth_modificado (G, massas, posicoes, momentos, eps, mat)
  REAL(pf), INTENT(IN)    :: G, eps
  REAL(pf), INTENT(INOUT) :: massas(:), posicoes(:,:), momentos(:,:)
  REAL(pf), INTENT(IN), OPTIONAL :: mat(:) ! momento angular total
  REAL(pf) :: energia = -0.25_pf  ! energia padrao
  REAL(pf) :: ec, ep              ! cinetica e potencial
  REAL(pf) :: Qv, beta, recond_vel
  REAL(pf) :: J(3)

  ! Normaliza as massas
  massas(:) = 1.0_pf/SIZE(massas)
  
  ! Para comecar, condiciona as integrais primeiras
  J = 0.0_pf
  IF (PRESENT(mat)) J = mat 
  IF (eps == 0) THEN
    CALL condicionar_ip_direto(G, massas, posicoes, momentos, 0.0_pf, 0.0_pf, J)
  ELSE
    CALL condicionar_ip_iterativo(G, massas, posicoes, momentos, eps, 0.0_pf, J)
  ENDIF

  ! Metodo de Aarseth
  ep = energia_potencial(G, massas, posicoes, eps)
  ec = energia_cinetica(massas, momentos)

  Qv = SQRT(0.5_pf * ABS(ep) / ec)
  beta = 0.5_pf * ep / energia

  ! Condiciona as velocidades
  momentos = momentos * Qv / SQRT(beta)

  ! Se nao tiver amortecimento, o metodo eh direto
  IF (eps == 0) THEN
    posicoes = posicoes * beta
  ! Se tiver, precisa ser iterativo, usando beta/eps como chute inicial
  ELSE
    CALL cond_pot_amortecido_equilibrio(energia, G, massas, posicoes, eps, beta)
    ep = energia_potencial(G, massas, posicoes, eps)

    ! Recondiciona os momentos lineares
    recond_vel = SQRT(ep/energia - 1.0_pf)
    momentos = momentos * recond_vel
  ENDIF
END SUBROUTINE

END MODULE condicionamento