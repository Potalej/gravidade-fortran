! *****************************************************************
!! Utilidades
!
! Objetivos:
!   Este arquivo contem helpers de mecanica e auxiliares. O modulo 
!   contem funcionalidades como o calculo de momento angular, 
!   energia, momento de inercia, etc. Ha tambem rotinas e funcoes 
!   auxiliares voltadas para o calculo vetorial, resolucao de 
!   sistemas e afins.
! 
!   Este modulo depende do OpenBLAS/LAPACK devido as rotinas dgesv
!   e dsyev.
! 
! Modificado:
!   02 de fevereiro de 2024 (criado)
!   08 de agosto de 2025 (modificado)
! 
! Autoria:
!   oap
! 
MODULE utilidades
  USE tipos
  IMPLICIT NONE
  EXTERNAL dgesv
  EXTERNAL dsyev

  INTERFACE energia_cinetica
    MODULE PROCEDURE energia_cinetica_vec
    MODULE PROCEDURE energia_cinetica_esc
  END INTERFACE

  INTERFACE energia_potencial
    MODULE PROCEDURE energia_potencial_vec
    MODULE PROCEDURE energia_potencial_esc
  END INTERFACE
  
  INTERFACE energia_total
    MODULE PROCEDURE energia_total_vec
    MODULE PROCEDURE energia_total_esc
  END INTERFACE
  
  INTERFACE momento_inercia
    MODULE PROCEDURE momento_inercia_vec
    MODULE PROCEDURE momento_inercia_esc
  END INTERFACE
CONTAINS

! ************************************************************
!! Produto vetorial
!
! Objetivos:
!   Calcula o produto vetorial entre dois vetores do R3.
!
! Modificado:
!   03 de fevereiro de 2024
!
! Autoria:
!   oap
! 
FUNCTION produto_vetorial (u, v)

  IMPLICIT NONE
  REAL(pf), DIMENSION(:), INTENT(IN) :: u, v
  REAL(pf), DIMENSION(3)             :: produto_vetorial

  produto_vetorial(:) = 0.0_pf

  produto_vetorial(1) =  u(2)*v(3)-v(2)*u(3)
  produto_vetorial(2) = -u(1)*v(3)+v(1)*u(3)
  produto_vetorial(3) =  u(1)*v(2)-v(1)*u(2)

END FUNCTION produto_vetorial

! ************************************************************
!! Determinante de matriz 3x3
!
! Objetivos:
!   Calcula o determinante de uma matriz 3x3.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
REAL FUNCTION determinante (M) RESULT (det)

  REAL(pf), DIMENSION(3,3), INTENT(IN) :: M
  det = M(1,1)*(-M(3,2)*M(2,3)) - M(1,2)*(M(2,1)*M(3,3)-M(2,3)*M(3,1)) + M(1,3)*(M(2,1)*M(3,2)-M(2,2)*M(3,1))
  det = det + M(1,1)*M(2,2)*M(3,3)

END FUNCTION determinante

! ************************************************************
!! Sistema linear de 3 equacoes (dgesv)
!
! Objetivos:
!   Resolve um sistema de 3 equacoes lineares usando o dgesv
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
FUNCTION sistema_linear3 (A, b)

  IMPLICIT NONE
  REAL(pf) :: A(3,3), b(3)
  REAL(pf) :: sistema_linear3(3), matriz(3,3)
  INTEGER  :: PIVOS(3), INFO

  sistema_linear3 = b
  matriz = A

  CALL dgesv(3,1,matriz,3,PIVOS,sistema_linear3,3,INFO)

  IF (INFO < 0) THEN
    WRITE (*, '(a)') 'O ', -INFO, '-ÉSIMO PARAMETRO TEM UM VALOR ILEGAL'
  ELSE IF (INFO > 0) then
    WRITE (*, '(a)') 'MATRIZ SINGULAR! SEM SOLUCAO'
  ENDIF
END FUNCTION sistema_linear3

! ************************************************************
!! Ordenacao de lista indexada do menor ao maior (via bubble)
!
! Objetivos:
!   Indexar uma lista e retornar a lista de indices
!
! Modificado:
!   30 de abril de 2025
!
! Autoria:
!   oap
! 
FUNCTION ordenar_lista_crescente (valores)

  IMPLICIT NONE
  REAL(pf), INTENT(IN)   :: valores(:)
  REAL(pf) :: valor_temp
  INTEGER, ALLOCATABLE :: ordenar_lista_crescente(:), idx(:)
  INTEGER :: N, i, j

  N = SIZE(valores)
  ALLOCATE(idx(N))
  ALLOCATE(ordenar_lista_crescente(N))

  DO i = 1, N
    idx(i) = i
  END DO

  DO i = 1, N-1
    DO j = i+1, N
      IF (valores(idx(j)) < valores(idx(i))) THEN
        valor_temp = idx(i)
        idx(i) = idx(j)
        idx(j) = valor_temp
      ENDIF
    END DO
  END DO

  ordenar_lista_crescente = idx

END FUNCTION ordenar_lista_crescente

! ************************************************************
!! Autovalores de uma matriz quadrada qualquer
!
! Modificado:
!   08 de agosto de 2025
!
! Autoria:
!   oap
! 
FUNCTION autovalores_matriz (matriz, tamanho) RESULT(autovalores)
  INTEGER, INTENT(IN) :: tamanho
  REAL(pf), INTENT(IN), DIMENSION(tamanho,tamanho) :: matriz  
  REAL(pf), DIMENSION(tamanho) :: autovalores

  REAL(pf64), ALLOCATABLE :: workspace(:)
  INTEGER :: info, lwork

  !> Usando o DSYEV para calcular autovalores e autovetores
  lwork = -1
  ALLOCATE(workspace(1))
  call DSYEV('N', 'U', tamanho, matriz, tamanho, autovalores, workspace, lwork, info)
  lwork = INT(workspace(1))
  DEALLOCATE(workspace)
  ALLOCATE(workspace(lwork))

  CALL DSYEV('N', 'U', tamanho, matriz, tamanho, autovalores, workspace, lwork, info)
END FUNCTION autovalores_matriz

! ************************************************************
!! Momento angular individual
!
! Objetivos:
!   Calcula o momento angular de um corpo dado sua posicao e
!   seu momento linear.
!
! Modificado:
!   03 de fevereiro de 2024
!
! Autoria:
!   oap
! 
FUNCTION momento_angular_individual (R, P)
  IMPLICIT NONE
  REAL(pf), DIMENSION(3), INTENT(IN) :: R, P
  REAL(pf), DIMENSION(3)             :: momento_angular_individual
  momento_angular_individual = produto_vetorial(R,P)
END FUNCTION momento_angular_individual

! ************************************************************
!! Momento angular total
!
! Objetivos:
!   Calcula o momento angular de uma lista de corpos dadas as
!   suas posicoes e seus momentos lineares.
!
! Modificado:
!   03 de fevereiro de 2024
!
! Autoria:
!   oap
!
FUNCTION momento_angular_total (Rs, Ps)
  IMPLICIT NONE
  REAL(pf), DIMENSION(:,:), INTENT(IN), TARGET :: Rs, Ps
  REAL(pf), DIMENSION(3)               :: momento_angular_total
  INTEGER                              :: i
  REAL(pf), POINTER :: R_p(:), P_p(:)

  momento_angular_total = (/0.0_pf,0.0_pf,0.0_pf/)
  DO i=1, SIZE(Rs,1)
    R_p => Rs(i,:)
    P_p => Ps(i,:)
    momento_angular_total = momento_angular_total &
                          + momento_angular_individual(R_p, P_p)
  END DO
END FUNCTION momento_angular_total

! ************************************************************
!! Energia cinetica
!
! Objetivos:
!   Calcula a energia cinetica de um conjunto de corpos dadas
!   as suas massas e momentos linares.
!
! Modificado:
!   03 de fevereiro de 2024
!
! Autoria:
!   oap
!
FUNCTION energia_cinetica_vec (m, P) RESULT(ec)
  IMPLICIT NONE
  REAL(pf) :: m(:), P(:,:), ec
  INTEGER  :: i
  ec=0.0_pf
  DO i=1, SIZE(m)
    ec = ec + DOT_PRODUCT(P(i,:), P(i,:))/m(i)
  END DO
  ec = 0.5_pf * ec
END FUNCTION energia_cinetica_vec

FUNCTION energia_cinetica_esc (m, P) RESULT(ec)
  IMPLICIT NONE
  REAL(pf), INTENT(IN) :: m, P(:,:)
  REAL(pf) :: m_inv, ec
  m_inv = 1 / m
  ec = 0.5_pf * sum(P*P) * m_inv
END FUNCTION energia_cinetica_esc

! ************************************************************
!! Energia potencial
!
! Objetivos:
!   Calcula a energia potencial newtoniana dada a constante
!   de gravitacao universal, as massas e as posicoes.
!
!   A formula utilizada eh a seguinte:
!   V = - G sum (i < j) m_i m_j / r_ij
!
! Modificado:
!   21 de julho de 2025
!
! Autoria:
!   oap
!
FUNCTION energia_potencial_vec (G,m,R,eps,dists) RESULT(ep)
  IMPLICIT NONE
  REAL(pf) :: m(:), R(:,:)
  REAL(pf) :: distancia, G, ep, distancia_inv
  REAL(pf), OPTIONAL :: eps
  REAL(pf), OPTIONAL :: dists(:)
  INTEGER  :: i,j,indice
  ep=0.0_pf
  indice = 1
  IF (PRESENT(dists)) THEN
    DO i=2, SIZE(m)
      DO j=1,i-1        
        ! distancia = dists(indice)
        distancia = SQRT(dists(indice)*dists(indice) + eps*eps)
        indice = indice + 1
        distancia_inv = 1.0_pf/distancia
        ep = ep + m(i)*m(j)*distancia_inv
      END DO
    END DO
  ELSE
    DO i=2, SIZE(m)
      DO j=1,i-1
        ! distancia = NORM2(R(i,:)-R(j,:))
        distancia = SQRT(DOT_PRODUCT(R(i,:)-R(j,:),R(i,:)-R(j,:)) + eps*eps)
        distancia_inv = 1.0_pf/distancia
        ep = ep + m(i)*m(j)*distancia_inv
      END DO
    END DO
  ENDIF
  ep = -G*ep
END FUNCTION energia_potencial_vec

FUNCTION energia_potencial_esc (G,m,R,eps,dists) RESULT(ep)
  IMPLICIT NONE
  REAL(pf) :: m2, distancia, ep, distancia_inv
  REAL(pf), INTENT(IN) :: G, m, R(:,:)
  REAL(pf), INTENT(IN), OPTIONAL :: eps
  REAL(pf), INTENT(IN), OPTIONAL :: dists(:)
  INTEGER :: i,j
  m2 = m * m
  ep=0.0_pf
  IF (PRESENT(dists)) THEN
    DO i = 1, SIZE(dists)
      ep = ep + 1.0_pf/SQRT(dists(i)*dists(i) + eps*eps)
    END DO
  ELSE
    DO i=2, SIZE(R,1)
      DO j=1,i-1
        ! distancia = NORM2(R(i,:)-R(j,:))
        distancia = SQRT(DOT_PRODUCT(R(i,:)-R(j,:),R(i,:)-R(j,:)) + eps*eps)
        distancia_inv = 1.0_pf/distancia
        ep = ep + distancia_inv
      END DO
    END DO
  ENDIF
  ep = -G*ep*m2
END FUNCTION energia_potencial_esc

! ************************************************************
!! Energia total do problema de N corpos
!
! Objetivos:
!   Calcula a energia total do problema de N corpos dados o
!   o valor da constante de gravitacao universal, as massas,
!   as posicoes e os momentos lineares.
!
!   A energia total eh dada pela soma da energia cinetica com
!   a energia potencial newtoniana
!
! Modificado:
!   20 de julho de 2025
!
! Autoria:
!   oap
!
FUNCTION energia_total_vec (G, m, R, P, eps, dists) RESULT(e_tot)
  IMPLICIT NONE
  REAL(pf) :: m(:), R(:,:), P(:,:), G
  REAL(pf), OPTIONAL :: eps
  REAL(pf), OPTIONAL :: dists(:)
  REAL(pf) :: e_tot
  IF (.NOT. PRESENT(eps)) eps = 0.0_pf
  IF (PRESENT(dists)) THEN
    e_tot = energia_cinetica(m,P) + energia_potencial(G,m,R,eps,dists)
  ELSE
    e_tot = energia_cinetica(m,P) + energia_potencial(G,m,R,eps)
  ENDIF
END FUNCTION energia_total_vec

FUNCTION energia_total_esc (G, m, R, P, eps, dists) RESULT(e_tot)
  IMPLICIT NONE
  REAL(pf) :: m
  REAL(pf) :: R(:,:), P(:,:), G
  REAL(pf), OPTIONAL :: eps
  REAL(pf), OPTIONAL :: dists(:)
  REAL(pf) :: e_tot
  IF (.NOT. PRESENT(eps)) eps = 0.0_pf
  IF (PRESENT(dists)) THEN
    e_tot = energia_cinetica(m,P) + energia_potencial(G,m,R,eps,dists)
  ELSE
    e_tot = energia_cinetica(m,P) + energia_potencial(G,m,R,eps)
  ENDIF
END FUNCTION energia_total_esc

! ************************************************************
!! Momento de dilatacao
!
! Objetivos:
!   Calcula o momento de dilatacao.
!
! Modificado:
!   05 de janeiro de 2025
!
! Autoria:
!   oap
!
FUNCTION momento_dilatacao (R, P)
  IMPLICIT NONE
  REAL(pf) :: R(:,:), P(:,:)
  REAL(pf) :: momento_dilatacao
  INTEGER  :: i
  momento_dilatacao=0.0_pf
  DO i=1, SIZE(R,1)
    momento_dilatacao = momento_dilatacao + DOT_PRODUCT(R(i,:),P(i,:))
  END DO
END FUNCTION momento_dilatacao

! ************************************************************
!! Momento de inercia
!
! Objetivos:
!   Calcula o momento de inercia.
!
! Modificado:
!   07 de janeiro de 2025
!
! Autoria:
!   oap
!
FUNCTION momento_inercia_vec (m, R) RESULT(momine)
  IMPLICIT NONE
  REAL(pf) :: m(:), R(:,:)
  REAL(pf) :: momine
  INTEGER  :: i
  momine=0.0_pf
  DO i=1, SIZE(R,1)
    momine = momine + m(i) * DOT_PRODUCT(R(i,:),R(i,:))
  END DO
END FUNCTION momento_inercia_vec

FUNCTION momento_inercia_esc (m, R) RESULT(momine)
  IMPLICIT NONE
  REAL(pf) :: m, R(:,:)
  REAL(pf) :: momine
  momine=m*sum(R(:,1)**2 + R(:,2)**2 + R(:,3)**2)
END FUNCTION momento_inercia_esc

! ************************************************************
!! Raio de meia massa (half-mass ratio)
!
! Objetivos:
!   Calcula o raio de meia massa.
!
! Modificado:
!   30 de abril de 2025
!
! Autoria:
!   oap
!
FUNCTION raio_meia_massa (m, R)
  IMPLICIT NONE
  REAL(pf) :: m(:), R(:,:)
  REAL(pf) :: raio_meia_massa, M_met, M_soma, qcm(3)
  INTEGER :: i, j
  INTEGER, ALLOCATABLE :: idx_ord(:)
  REAL(pf), ALLOCATABLE :: raios(:)

  ALLOCATE(idx_ord(SIZE(m)))
  ALLOCATE(raios(SIZE(m)))
  raios = 0.0_pf
  
  ! Calcula os raios
  qcm = centro_massas(m, R)
  DO i=1, SIZE(m)
    raios(i) = NORM2(R(i,:) - qcm)
  END DO

  ! Ordena a lista em ordem crescente
  idx_ord = ordenar_lista_crescente(raios)

  ! Agora vai pegando as massas ate chegar na metade
  M_met = 0.5_pf * SUM(m)
  M_soma = 0.0_pf
  DO i=1, SIZE(m)
    j = idx_ord(i)
    M_soma = M_soma + m(j)
    IF (M_soma .GE. M_met) THEN
      raio_meia_massa = raios(j)
      EXIT
    END IF
  END DO

END FUNCTION raio_meia_massa

! ************************************************************
!! Tempo de relaxacao do raio de meia massa t_rh
!
! Objetivos:
!   Calcula o tempo de relaxacao do raio de meia massa t_rh
!
! Modificado:
!   30 de abril de 2025
!
! Autoria:
!   oap
!
FUNCTION tempo_relaxacao_rh (m, R)
  IMPLICIT NONE
  REAL(pf) :: G=1.0_pf, m(:), R(:,:), tempo_relaxacao_rh
  REAL(pf) :: rh, m_media, log_coulomb
  INTEGER :: N
  
  N = SIZE(m)
  ! Calcula o raio de meia massa
  rh = raio_meia_massa(m, R)
  WRITE(*,*) 'rh:', rh

  ! Media das massas
  m_media = SUM(m) / N

  ! Logaritmo de coulomb
  log_coulomb = LOG(0.4_pf * N)

  ! Aproximacao do valor de t_rh
  tempo_relaxacao_rh = (0.138_pf * N / log_coulomb) * SQRT(rh*rh*rh / G*m_media)

END FUNCTION tempo_relaxacao_rh

! ************************************************************
!! Calculo do termo <F,q> para o virial amortecido
!
! Modificado:
!   25 de julho de 2025
!
! Autoria:
!   oap
!
FUNCTION virial_potencial_amortecido (G, m, R, eps, ep_var) RESULT(f_prod_q)
  REAL(pf), INTENT(IN) :: G, m(:), R(:,:), eps
  REAL(pf), INTENT(OUT), OPTIONAL :: ep_var  ! potencial
  REAL(pf) :: Fab(3), dist2, denominador, f_prod_q, ep
  INTEGER :: a, b
  
  f_prod_q = 0.0_pf
  ep = 0.0_pf
  DO a = 2, SIZE(m)
    DO b = 1, a - 1
      Fab = R(b,:) - R(a,:)
      dist2 = DOT_PRODUCT(Fab, Fab)
      denominador = SQRT(dist2 + eps*eps)

      ! potencial
      ep = ep - G * m(a) * m(b) / denominador

      ! forca
      Fab = G * m(a) * m(b) * Fab / (denominador*denominador*denominador)

      ! virial
      f_prod_q = f_prod_q + DOT_PRODUCT(R(a,:), Fab)
      f_prod_q = f_prod_q - DOT_PRODUCT(R(b,:), Fab)
    END DO
  END DO

  IF (PRESENT(ep_var)) ep_var = ep

END FUNCTION

! ************************************************************
!! Tensor de inercia
!
! Objetivos:
!   Dada a massa e a posicao de um corpo, calcula o tensor de
!   inercia respectivo.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
FUNCTION tensor_inercia (m, R)

  IMPLICIT NONE
  REAL(pf), DIMENSION(3), INTENT(IN) :: R
  REAL(pf), INTENT(IN)               :: m
  REAL(pf), DIMENSION(3,3) :: tensor_inercia
  INTEGER              :: a, b

  DO a = 1, 3
    DO b = 1, 3
      IF (a /= b) THEN
        tensor_inercia(a,b) = m * R(a) * R(b)
      ENDIF
    END DO
  END DO

  tensor_inercia(1,1) = - m * (R(2)**2 + R(3)**2)
  tensor_inercia(2,2) = - m * (R(1)**2 + R(3)**2)
  tensor_inercia(3,3) = - m * (R(1)**2 + R(2)**2)

END FUNCTION tensor_inercia

! ************************************************************
!! Tensor de inercia geral
!
! Objetivos:
!   Dadas as massas e posicoes dos corpos, calcula o tensor de
!   inercia geral, que eh basicamente a soma dos tensores de
!   inercia individuais.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
FUNCTION tensor_inercia_geral (massas, posicoes)

  IMPLICIT NONE
  REAL(pf), DIMENSION(:), INTENT(IN) :: massas
  INTEGER :: a
  REAL(pf), DIMENSION(SIZE(massas),3), TARGET :: posicoes
  REAL(pf), POINTER :: pos_p(:)
  REAL(pf), DIMENSION(3,3) :: tensor_inercia_geral

  tensor_inercia_geral(:,:) = 0.0_pf
  
  DO a = 1, SIZE(massas)
    pos_p => posicoes(a,:)
    tensor_inercia_geral = tensor_inercia_geral + tensor_inercia(massas(a), pos_p)
  END DO   

END FUNCTION tensor_inercia_geral

! ************************************************************
!! Centro de massas
!
! Objetivos:
!   Dadas as massas e posicoes dos corpos, calcula o centro de
!   massas, dado pela media dos corpos ponderada pelas massas.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
FUNCTION centro_massas (massas, posicoes)

  IMPLICIT NONE
  REAL(pf), INTENT(IN)   :: massas(:), posicoes(:,:)
  REAL(pf), DIMENSION(3) :: centro_massas
  INTEGER                :: a

  centro_massas(:) = 0.0_pf

  DO a = 1, SIZE(massas)
    centro_massas = centro_massas + massas(a) * posicoes(a,:)
  END DO
  centro_massas = centro_massas / SUM(massas)

END FUNCTION centro_massas

! ************************************************************
!! Momento linear total
!
! Objetivos:
!   Soma os momentos lineares e retorna tal vetor.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
FUNCTION momento_linear_total (momentos)

  IMPLICIT NONE
  REAL(pf), INTENT(IN)   :: momentos(:,:)
  REAL(pf), DIMENSION(3) :: momento_linear_total
  INTEGER                :: a
  momento_linear_total = 0.0_pf
  DO a = 1, SIZE(momentos,1)
    momento_linear_total = momento_linear_total + momentos(a,:)
  END DO

END FUNCTION momento_linear_total

! ************************************************************
!! Anisotropia do tensor de inércia
!
! Objetivos:
!   Calcula a anisotropia do tensor de inercia, indicando a
!   esfericidade do sistema.
!
! Modificado:
!   30 de abril de 2025
!
! Autoria:
!   oap
! 
FUNCTION anisotropia_tensor_inercia (m, R)

  IMPLICIT NONE
  REAL(pf), INTENT(IN)   :: m(:), R(:,:)
  REAL(pf) :: tensor(3,3), autovalores(3), l1, l2, l3
  INTEGER :: idx(3)
  REAL(pf) :: anisotropia_tensor_inercia 

  !> Calcula o tensor e seus autovalores
  tensor = -tensor_inercia_geral(m, R)
  autovalores = autovalores_matriz(tensor, 3)

  ! Ordena
  idx = ordenar_lista_crescente(autovalores)
  l1 = autovalores(idx(3))
  l2 = autovalores(idx(2))
  l3 = autovalores(idx(1))

  WRITE (*,*) '     * T.I.  =', tensor(1,1), tensor(2,2), tensor(3,3)
  WRITE (*,*) '     * Autovalores (T.I.): ', l3, l2, l1
  
  anisotropia_tensor_inercia = (l2 - l3)/l1

END FUNCTION anisotropia_tensor_inercia

! ************************************************************
!! Anisotropia via velocidades radial e tangencial
!
! Objetivos:
!  Anisotropia via velocidades radial e tangencial
!
! Modificado:
!   30 de abril de 2025
!
! Autoria:
!   oap
! 
FUNCTION anisotropia_velocidades (m, R, P)

  IMPLICIT NONE
  REAL(pf) :: anisotropia_velocidades
  REAL(pf), INTENT(IN) :: m(:), R(:,:), P(:,:)
  REAL(pf) :: P_radial, P_tangente
  REAL(pf) :: media_radial, media_tangente
  INTEGER :: a

  media_radial = 0.0_pf
  media_tangente = 0.0_pf

  DO a=1, SIZE(m)
    P_radial = DOT_PRODUCT(P(a,:), R(a,:)) / NORM2(R(a,:))
    P_tangente = SQRT(DOT_PRODUCT(P(a,:), P(a,:)) - P_radial**2)

    media_radial = media_radial + P_radial * P_radial / m(a)
    media_tangente = media_tangente + P_tangente * P_tangente / m(a)
  END DO

  media_radial = media_radial / SUM(m)
  media_tangente = media_tangente / SUM(m)

  anisotropia_velocidades = 1 - media_tangente / (media_radial + media_radial)

END FUNCTION anisotropia_velocidades

END MODULE utilidades