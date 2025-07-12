! *****************************************************************
!! FUNCOES AUXILIARES GERAIS
!
! Objetivos:
!   Este arquivo contem helpers de varias coisas, como vetores e
!   tal.
!
! Modificado:
!   15 de marco de 2024
! 
! Autoria:
!   oap
! 
MODULE auxiliares
  USE tipos
  IMPLICIT NONE
  EXTERNAL dgesv
  EXTERNAL dsyev
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
FUNCTION momentoLinear_total (momentos)

  IMPLICIT NONE
  REAL(pf), INTENT(IN)   :: momentos(:,:)
  REAL(pf), DIMENSION(3) :: momentoLinear_total
  INTEGER                :: a
  momentoLinear_total = 0.0_pf
  DO a = 1, SIZE(momentos,1)
    momentoLinear_total = momentoLinear_total + momentos(a,:)
  END DO

END FUNCTION momentoLinear_total

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
  INTEGER, parameter :: N = 3
  REAL(pf) :: tensor(N,N), autovalores(N), l1, l2, l3
  REAL(pf64), ALLOCATABLE :: workspace(:)
  INTEGER :: info, lwork, idx(3)
  REAL(pf) :: anisotropia_tensor_inercia 

  tensor = -tensor_inercia_geral(m, R)

  ! Call DSYEV to compute eigenvalues and eigenvectors
  lwork = -1
  ALLOCATE(workspace(1))
  call DSYEV('N', 'U', N, tensor, N, autovalores, workspace, lwork, info)
  lwork = INT(workspace(1))
  DEALLOCATE(workspace)
  ALLOCATE(workspace(lwork))

  CALL dsyev('N', 'U', N, tensor, N, autovalores, workspace, lwork, info)

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

END MODULE auxiliares