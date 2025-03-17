! *****************************************************************
!! CONDICIONAMENTO DE VALORES INICIAIS
!
! Objetivos:
!   Funcoes para o condicionamento de particulas.
! 
! Modificado:
!   15 de janeiro de 2025
! 
! Autoria:
!   oap
! 
MODULE condicoesIniciais
  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  USE mecanica
  USE auxiliares
  USE aleatorio
  IMPLICIT NONE
CONTAINS

! ************************************************************
!! Gera vetores 3d
!
! Objetivos:
!   Gera vetores 3d utilizados para as posicoes e velocidades.
!
! Modificado:
!   15 de janeiro de 2025
!
! Autoria:
!   oap
! 
FUNCTION gerar_vetores3d (N, min, max, vid, vire, vira)

  INTEGER, INTENT(INOUT)   :: N
  REAL(pf), INTENT(IN)     :: min, max
  REAL(pf), DIMENSION(N,3) :: gerar_vetores3d
  INTEGER :: i, j
  REAL(pf) :: maior_valor = 1.0_pf
  CHARACTER(20), INTENT(IN) :: vid, vire ! vi_dist, vi_regiao
  REAL(pf), INTENT(IN)  :: vira ! vi_raio

  SELECT CASE (TRIM(vid))
    ! Uniforme (0,1)
    CASE ("uniforme"); CALL uniforme(gerar_vetores3d, N, 0.0_pf, min, max, vire, vira)
    ! Normal (0,1)
    CASE ("normal"); CALL normal(gerar_vetores3d, N, 0.0_pf, min, max, vire, vira)
    ! Cauchy
    CASE ("cauchy"); CALL cauchy(gerar_vetores3d, N, 0.0_pf, min, max, vire, vira)
  END SELECT

END FUNCTION gerar_vetores3d

! ************************************************************
!! Gera vetores 3d intervalado
!
! Objetivos:
!   Gera vetores 3d intervalados, ou seja, com uma distancia minima
!   entre eles.
!
! Modificado:
!   15 de janeiro de 2025
!
! Autoria:
!   oap
! 
FUNCTION gerar_vetores3d_intervalado (N, min, max, distmin, vid, vire, vira)

  INTEGER, INTENT(INOUT)   :: N
  REAL(pf), INTENT(IN)     :: min, max
  REAL(pf), DIMENSION(N,3) :: gerar_vetores3d_intervalado
  INTEGER :: i, j
  REAL(pf), INTENT(IN) :: distmin ! distancia minima
  CHARACTER(20), INTENT(IN) :: vid, vire ! vi_dist, vi_regiao
  REAL(pf), INTENT(IN)  :: vira ! vi_raio

  SELECT CASE (TRIM(vid))
    ! Uniforme (0,1)
    CASE ("uniforme"); CALL uniforme(gerar_vetores3d_intervalado, N, distmin, min, max, vire, vira)
    ! Normal (0,1)
    CASE ("normal"); CALL normal(gerar_vetores3d_intervalado, N, distmin, min, max, vire, vira)
    ! Cauchy
    CASE ("cauchy"); CALL cauchy(gerar_vetores3d_intervalado, N, distmin, min, max, vire, vira)
  END SELECT

END FUNCTION gerar_vetores3d_intervalado

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
  REAL(pf)               :: ajuste(N)

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
SUBROUTINE gerarValores (N, massas, posicoes, momentos, int_posicoes, int_momentos, int_massas, vid, vire, vira)

  IMPLICIT NONE
  INTEGER, INTENT(INOUT)      :: N
  REAL(pf), DIMENSION(2), INTENT(IN) :: int_momentos, int_massas
  REAL(pf), DIMENSION(3), INTENT(IN) :: int_posicoes
  REAL(pf), INTENT(INOUT) :: posicoes(N,3), momentos(N,3), massas(N)
  CHARACTER(20), INTENT(IN) :: vid, vire ! vi_dist, vi_regiao
  REAL(pf), INTENT(IN)  :: vira ! vi_raio

  ! Gera massas
  WRITE (*,*) '    * gerando massas'
  massas = gerar_massas(N, int_massas(1), int_massas(2))

  ! Gera as posições
  WRITE (*,*) '    * gerando posicoes'
  WRITE (*,*) '       * Distribuicao:', vid, vire, vira, int_posicoes(3)
  posicoes = gerar_vetores3d_intervalado(N, int_posicoes(1), int_posicoes(2), int_posicoes(3), vid, vire, vira)

  ! Gera os momentos
  WRITE (*,*) '    * gerando momentos'
  WRITE (*,*) '       * Distribuicao:', vid, vire, vira
  momentos = gerar_vetores3d(N, int_momentos(1), int_momentos(2), vid, 'cubo                ', vira)

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
!   10 de novembro de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE condicionar_ip (G, massas, posicoes, momentos, H, J, P)
  
  IMPLICIT NONE
  REAL(pf), INTENT(IN) :: G
  REAL(pf), INTENT(INOUT) :: posicoes(:,:), momentos(:,:), massas(:)
  REAL(pf) :: H, J(3), P(3)
  REAL(pf) :: erro_0, erro_1, ENERGIA, LINEAR(3), ANGULAR(3) ! Para medir a taxa de erro
  INTEGER :: i = 0
  REAL(pf) :: erro_limite = 0.1E-8

  ! Zera o centro de massas
  CALL zerar_centroMassas(massas, posicoes)
      
  ! Condiciona o momento linear
  CALL condicionar_momentoLinear(P, massas, momentos)

  ! Condiciona o momento angular
  CALL condicionar_momentoAngular(J, massas, posicoes, momentos)

  ! ! Condiciona a energia total
  CALL condicionar_energiaTotal(H, G, massas, posicoes, momentos)

  ! Calculo dos erros
  ENERGIA = ABS(energia_total(G,massas,posicoes,momentos) - H)
  LINEAR  = ABS(momentoLinear_total(momentos) - P)
  ANGULAR = ABS(momento_angular_total(posicoes,momentos) - J)
  erro_1 = MAXVAL((/ENERGIA,LINEAR(1),LINEAR(2),LINEAR(3),ANGULAR(1),ANGULAR(2),ANGULAR(3)/))

  IF (erro_1 >= erro_limite) THEN
    erro_0 = ABS(erro_1) + 1.0 ! Para entrar no loop abaixo
    
    DO WHILE (ABS(erro_1) >= erro_limite .AND. i <= 10)
      i = i + 1
      erro_0 = erro_1

      ! Condiciona o momento linear
      CALL condicionar_momentoLinear(P, massas, momentos)

      ! Condiciona o momento angular
      CALL condicionar_momentoAngular(J, massas, posicoes, momentos)

      ! ! Condiciona a energia total
      CALL condicionar_energiaTotal(H, G, massas, posicoes, momentos)

      ! Calculo dos erros
      ENERGIA = ABS(energia_total(G,massas,posicoes,momentos) - H)
      LINEAR  = ABS(momentoLinear_total(momentos) - P)
      ANGULAR = ABS(momento_angular_total(posicoes,momentos) - J)
      erro_1 = MAXVAL((/ENERGIA,LINEAR(1),LINEAR(2),LINEAR(3),ANGULAR(1),ANGULAR(2),ANGULAR(3)/))
    END DO
  ENDIF

  WRITE (*,*) ' > condicionamento aplicado ', i, ' vezes para obter o erro ', erro_1
END SUBROUTINE condicionar_ip

! ************************************************************
!! Condicionamento de integrais primeiras na forma direta
!
! Objetivos:
!   Condiciona vetores ja existentes com valores desejados na
!   forma direta, sem a necessidade de iteracoes.
!
! Modificado:
!   17 de novembro de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE condicionar_ip_2 (G, m, Rs, Ps, H, J, P)
  
  IMPLICIT NONE
  REAL(pf), INTENT(IN) :: G
  REAL(pf), INTENT(INOUT) :: Rs(:,:), Ps(:,:), m(:)
  REAL(pf) :: H, J(3), P(3), M_tot_inv
  REAL(pf) :: ENERGIA, LINEAR(3), ANGULAR(3) ! Para medir a taxa de erro
  INTEGER :: i, a = 0
  REAL(pf) :: alpha, V0
  REAL(pf) :: beta, S1, S2, tensorInercia(3,3), rot(3), rot_(3), rot_tot(3), K1a(3), K2a(3)
  REAL(pf), ALLOCATABLE :: Rs1(:,:), Ps1(:,:)

  ALLOCATE(Rs1(SIZE(m), 3))
  ALLOCATE(Ps1(SIZE(m), 3))

  ! Zera o centro de massas
  CALL zerar_centroMassas(m, Rs)

  ! Calcula as integrais primeiras
  ENERGIA = energia_total(G, m, Rs, Ps)
  LINEAR = momentoLinear_total(Ps)
  ANGULAR = momento_angular_total(Rs, Ps)
  M_tot_inv = 1.0_pf / SUM(m)
  
  WRITE (*,*) 
  WRITE (*,*) ' > integrais primeiras iniciais'
  WRITE (*,*) '   * E   =', ENERGIA
  WRITE (*,*) '   * J   =', ANGULAR
  WRITE (*,*) '   * P   =', LINEAR
  WRITE (*,*) '   * Qcm =', centro_massas(m, Rs)

  ! Calculo do alpha
  V0 = energia_potencial(G, m, Rs)
  alpha = 1.0_pf + H / V0

  ! Calculo do beta
  tensorInercia = tensor_inercia_geral(m, Rs)
  rot     = sistema_linear3(tensorInercia, ANGULAR)
  rot_    = sistema_linear3(tensorInercia, J) * alpha
  
  S1 = 0.0_pf
  S2 = 0.0_pf

  DO a = 1, SIZE(m)

    K1a = Ps(a,:) - m(a) * LINEAR * M_tot_inv - m(a) * produto_vetorial(Rs(a,:), rot)
    S1 = S1 + 0.5_pf * DOT_PRODUCT(K1a, K1a) / m(a)

    K2a = m(a) * P * M_tot_inv + m(a) * produto_vetorial(Rs(a,:), rot_)
    S2 = S2 + 0.5_pf * DOT_PRODUCT(K2a, K2a) / m(a)

  END DO

  beta = SQRT((-V0 - S2)/S1)

  WRITE (*,*)
  WRITE (*,*) ' > coeficientes:'
  WRITE (*,*) '   * alpha =', alpha
  WRITE (*,*) '   * beta  =', beta
  WRITE (*,*) '   * S1    =', S1
  WRITE (*,*) '   * S2    =', S2

  ! Transforma as coordenadas
  rot_tot = sistema_linear3(tensorInercia, ANGULAR - J * alpha / beta)

  Rs1 = Rs / alpha

  DO a = 1, SIZE(m)
    Ps1(a,:) = beta * (Ps(a,:) - m(a) * M_tot_inv * (LINEAR - P / beta) - m(a) * produto_vetorial(Rs(a,:), rot_tot))
  END DO

  Rs = Rs1
  Ps = Ps1

END SUBROUTINE condicionar_ip_2

! ************************************************************
!! Gera valores condicionados pelas integrais primeiras
!
! Objetivos:
!   Gera vetores aleatorios e os condicionado para valores
!   informados.
!
! Modificado:
!   03 de fevereiro de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE gerar_condicionado_ip (G, N, massas, pos, mom, int_pos, int_mom, int_massas, H, J, P, vid, vire, vira)

  IMPLICIT NONE
  INTEGER, INTENT(INOUT)  :: N
  REAL(pf), DIMENSION(2), INTENT(IN) :: int_mom, int_massas
  REAL(pf), DIMENSION(3), INTENT(IN) :: int_pos
  REAL(pf), INTENT(INOUT) :: G, H, J(3), P(3)
  REAL(pf), INTENT(INOUT), allocatable :: massas(:), pos(:,:), mom(:,:)
  CHARACTER(20), INTENT(IN) :: vid, vire ! vi_dist, vi_regiao
  REAL(pf), INTENT(IN)  :: vira ! vi_raio

  WRITE (*,'(a)') "GERACAO DAS CONDICOES INICIAIS"
  
  ! Gera os valores
  WRITE (*,'(a)') '  > gerando valores...'
  ALLOCATE(massas (N))
  ALLOCATE(pos (N, 3))
  ALLOCATE(mom (N, 3))
  CALL gerarValores(N, massas, pos, mom, int_pos, int_mom, int_massas, vid, vire, vira)

  ! Condiciona
  WRITE (*,'(a)') '  > condicionando...'
  CALL condicionar_ip(G, massas, pos, mom, H, J, P)

  ! Exibe as integrais primeiras do sistema
  WRITE (*,*) '    * H   =', energia_total(G,massas,pos,mom) 
  WRITE (*,*) '    * Rcm =', centro_massas(massas,pos) 
  WRITE (*,*) '    * P   =', momentoLinear_total(mom) 
  WRITE (*,*) '    * J   =', momento_angular_total(pos,mom) 

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
!   03 de fevereiro de 2024
!
! Autoria:
!   oap
!
SUBROUTINE gerar_condicionado_henon (N, massas, posicoes, momentos, int_posicoes, int_momentos, int_massas, vid, vire, vira)
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: N
  REAL(pf), DIMENSION(2), INTENT(IN) :: int_momentos, int_massas
  REAL(pf), DIMENSION(3), INTENT(IN) :: int_posicoes
  REAL(pf), INTENT(INOUT), allocatable :: massas(:), posicoes(:,:), momentos(:,:)
  CHARACTER(20), INTENT(IN) :: vid, vire ! vi_dist, vi_regiao
  REAL(pf), INTENT(IN)  :: vira ! vi_raio
  REAL(pf) :: beta, EP, EC, Qvir ! AARSETH
  REAL(pf) :: P_normas = 0.0_pf ! Soma dos quadrados das normas dos momentos lineares
  REAL(pf) :: E, dist_min, pv ! energia, disstancia minima, proporcao de virial
  integer :: a, b, m

  E = -0.25_pf ! Energia total final
  pv = 0.5_pf  ! Raio de virial

  WRITE (*,'(a)') "GERACAO DAS CONDICOES INICIAIS (HENON)"
  
  ! Gera os valores
  WRITE (*,'(a)') '  > gerando valores...'
  ALLOCATE(massas (N))
  ALLOCATE(posicoes (N, 3))
  ALLOCATE(momentos (N, 3))
  CALL gerarValores(N, massas, posicoes, momentos, int_posicoes, int_momentos, int_massas, vid, vire, vira)

  m = SUM(massas)
  do a = 1, N
    massas(a) = (1.0_pf / N)
  end do
  WRITE (*,*) "M:", SUM(massas)

  ! Condiciona
  WRITE (*,'(a)') '  > condicionando... (HENON)'

  CALL condicionar_ip (1.0_pf, massas, posicoes, momentos, E, (/0.0_pf,0.0_pf,0.0_pf/), (/0.0_pf,0.0_pf,0.0_pf/))

  ! AGORA EU QUERO QUE V ~ -0.5, OU SEJA, V = 2 * V
  EP = energia_potencial(1.0_pf, massas, posicoes)
  ! posicoes = posicoes * EP * 2.0_pf
  posicoes = posicoes * EP / (E * 2.0_pf) 
  EP = energia_potencial(1.0_pf, massas, posicoes)

  ! AGORA QUERO QUE T/(-V) = 0.5, LOGO v = sqrt(-V/M)
  EC = energia_cinetica(massas, momentos)
  momentos = momentos * SQRT(- pv * EP / EC)

  dist_min = 10000.0_pf
  DO a = 2, N
    DO b = 1, a-1
      IF (dist_min .gt. NORM2(posicoes(a,:) - posicoes(b,:))) THEN
        dist_min = NORM2(posicoes(a,:) - posicoes(b,:))
      ENDIF
    END DO
  END DO

  WRITE(*,*) 'Menor distancia:', dist_min

  EC = energia_cinetica(massas, momentos)
  WRITE (*,*)
  WRITE (*,*) ' > novos valores:'
  WRITE (*,*) '    * V   =', EP
  WRITE (*,*) '    * T   =', EC
  WRITE (*,*) '    * H   =', energia_total(1.0_pf,massas,posicoes,momentos) 
  WRITE (*,*) '    * qm  =', dist_min
  WRITE (*,*) '    * I   =', momento_inercia(massas,posicoes)
  WRITE (*,*) '    * Q   =', EC/ABS(EP)
  WRITE (*,*) '    * R   =', - 1.0_pf * SUM(massas)**2 / (2.0_pf * EP)
  WRITE (*,*) '    * P   =', momentoLinear_total(momentos)
  WRITE (*,*) '    * J   =', momento_angular_total(posicoes, momentos)

END SUBROUTINE

! ************************************************************
!! Gera valores condicionados pelo modelo de Henon usando
!  o algoritmo de Aarseth
!
! Objetivos:
!   Gera vetores aleatorios e os condicionado para valores
!   informados.
!
! Modificado:
!   14 de janeiro de 2025
!
! Autoria:
!   oap
!
SUBROUTINE gerar_condicionado_aarseth (N, massas, posicoes, momentos, int_posicoes, int_momentos, int_massas, vid, vire, vira)
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: N
  REAL(pf), DIMENSION(2), INTENT(IN) :: int_momentos, int_massas
  REAL(pf), DIMENSION(3), INTENT(IN) :: int_posicoes
  REAL(pf), INTENT(INOUT), allocatable :: massas(:), posicoes(:,:), momentos(:,:)
  CHARACTER(20), INTENT(IN) :: vid, vire ! vi_dist, vi_regiao
  REAL(pf), INTENT(IN)  :: vira ! vi_raio
  REAL(pf) :: beta, Qv, EP, EC, Qvir ! AARSETH
  REAL(pf) :: P_normas = 0.0_pf ! Soma dos quadrados das normas dos momentos lineares
  REAL(pf) :: E, dist_min, pv ! energia, disstancia minima, proporcao de virial
  integer :: a, b, m

  E = -0.25_pf
  pv = 0.5_pf

  WRITE (*,'(a)') "GERACAO DAS CONDICOES INICIAIS (AARSETH)"
  
  ! Gera os valores
  WRITE (*,'(a)') '  > gerando valores...'
  ALLOCATE(massas (N))
  ALLOCATE(posicoes (N, 3))
  ALLOCATE(momentos (N, 3))
  CALL gerarValores(N, massas, posicoes, momentos, int_posicoes, int_momentos, int_massas, vid, vire, vira)

  m = SUM(massas)
  do a = 1, N
    massas(a) = (1.0_pf / N)
  end do
  WRITE (*,*) "M:", SUM(massas)

  ! Condiciona
  WRITE (*,'(a)') '  > condicionando... (AARSETH)'

  CALL condicionar_ip (1.0_pf, massas, posicoes, momentos, E, (/0.0_pf,0.0_pf,0.0_pf/), (/0.0_pf,0.0_pf,0.0_pf/))

  WRITE (*,*) '    * I   =', momento_inercia(massas,posicoes)

  EP = energia_potencial(1.0_pf, massas, posicoes)
  EC = energia_cinetica(massas, momentos)

  ! agora aplicamos o metodo de Aarseth
  Qvir = 0.5_pf
  Qv = SQRT(Qvir * ABS(EP)/EC)
  momentos = momentos * Qv
  beta = (1.0_pf - Qvir) * EP / E
  posicoes = posicoes * beta
  momentos = momentos / SQRT(beta)

  dist_min = 10000.0_pf
  DO a = 2, N
    DO b = 1, a-1
      IF (dist_min .gt. NORM2(posicoes(a,:) - posicoes(b,:))) THEN
        dist_min = NORM2(posicoes(a,:) - posicoes(b,:))
      ENDIF
    END DO
  END DO

  WRITE(*,*) 'Menor distancia:', dist_min

  EP = energia_potencial(1.0_pf, massas, posicoes)
  EC = energia_cinetica(massas, momentos)
  WRITE (*,*)
  WRITE (*,*) ' > novos valores:'
  WRITE (*,*) '    * V   =', EP
  WRITE (*,*) '    * T   =', EC
  WRITE (*,*) '    * H   =', energia_total(1.0_pf,massas,posicoes,momentos) 
  WRITE (*,*) '    * qm  =', dist_min
  WRITE (*,*) '    * I   =', momento_inercia(massas,posicoes)
  WRITE (*,*) '    * D   =', momento_dilatacao(posicoes,momentos)
  WRITE (*,*) '    * Q   =', EC/ABS(EP)
  WRITE (*,*) '    * R   =', - 1.0_pf * SUM(massas)**2 / (2.0_pf * EP)
  WRITE (*,*) '    * P   =', momentoLinear_total(momentos)
  WRITE (*,*) '    * J   =', momento_angular_total(posicoes, momentos)

END SUBROUTINE

END MODULE condicoesIniciais