! ************************************************************
!! CORRECAO
!
! Calculo da correcao numerica via matriz jacobiana e utilizando
! das propriedades simpleticas do espaco de fases do sistema
! hamiltoniano.
! 
! Funcoes
! 
! = SUBROUTINE corrigir (G, massas, posicoes, momentos, grads, gradsT, vetorCorrecao, corrigiu)
! Funcao central, que aplica a correcao conforme descrito.
! 
! = FUNCTION Gx (G, massas, posicoes, momentos, N)
! Vetor de estados. Contem as dez integrais de movimento no momento pre-correcao.
! 
! = FUNCTION matriz_normal (G, massas, posicoes, momentos, N, grads)
! Matriz dada pelo produto interno dos gradientes de cada integral de movimento.
! Trata-se de uma matriz 4x4.
! 
! = FUNCTION gradiente_energia (G, massas, Rs, Ps, N)
! Calcula o gradiente da energia total do sistema. 
! 
! = FUNCTION gradiente_angularX (posicoes, momentos, N)
! Calcula o gradiente do componente X do momento angular.
! 
! = FUNCTION gradiente_angularY (posicoes, momentos, N)
! Calcula o gradiente do componente Y do momento angular.
! 
! = FUNCTION gradiente_angularZ (posicoes, momentos, N)
! Calcula o gradiente do componente Z do momento angular.
! 
! = FUNCTION gradiente_linearX (N)
! Gradiente do componente X do momento linear total do sistema.
! 
! = FUNCTION gradiente_linearY (N)
! Gradiente do componente Y do momento linear total do sistema.
! 
! = FUNCTION gradiente_linearZ (N)
! Gradiente do componente Z do momento linear total do sistema.
! 
! = FUNCTION gradiente_centroMassasX (massas, N)
! Gradiente do componente X do centro de massas.
! 
! = FUNCTION gradiente_centroMassasY (massas, N)
! Gradiente do componente Y do centro de massas.
! 
! = FUNCTION gradiente_centroMassasZ (massas, N)
! Gradiente do componente Z do centro de massas.
! 
! 
! Modificado:
!   10 de novembro de 2024
! 
! Autoria:
!   oap
! 
MODULE correcao
  USE tipos
  USE auxiliares
  USE mecanica

  IMPLICIT NONE
  PRIVATE
  PUBLIC corrigir, corrigir_apenas_energia

CONTAINS

! ************************************************************
!! Aplicacao
!
! Objetivos:
!   Aplicacao da correcao numerica nas integrais primeiras.
!
! Modificado:
!   10 de novembro de 2024
!
! Autoria:
!   oap
!
SUBROUTINE corrigir (corme, cormnt, G, massas, posicoes, momentos, corrigiu, H, J)
  IMPLICIT NONE
  REAL(pf), INTENT(IN)    :: corme ! CORrecao Margem Erro
  INTEGER,  INTENT(IN)    :: cormnt ! CORrecao Max Num Tentativas
  REAL(pf), INTENT(IN)    :: G, massas(:)
  REAL(pf), INTENT(INOUT) :: posicoes(:,:), momentos(:,:)
  LOGICAL,  INTENT(INOUT) :: corrigiu
  INTEGER                 :: N, a, INFO, b, contador = 0, pivos(4)
  REAL(pf)                :: JJt(4, 4), vetG(4)
  REAL(pf), ALLOCATABLE   :: grads(:,:), gradsT(:,:), vetorCorrecao(:)
  REAL(pf)                :: H, J(3)

  N = SIZE(massas)
  contador = 0

  ALLOCATE(grads(4,6*N))
  ALLOCATE(gradsT(6*N,4))
  ALLOCATE(vetorCorrecao(6*N))
  
  ! enquanto nao tiver aceitado a correcao, roda
  loop: DO WHILE (.TRUE.)

    IF (contador >= cormnt) THEN
      ! WRITE(*,*) "max num tentativas atingido"
      EXIT loop
    ENDIF

    ! Contador de correcoes
    contador = contador + 1

    ! calcula a mariz normal
    JJt = matriz_normal(G, massas, posicoes, momentos, N, grads)

    ! vetor G
    vetG = Gx(G, massas, posicoes, momentos, N) + (/H, J(1), J(2), J(3)/)

    ! resolve o sistema via fatoracao de Cholesky
    CALL dpotrf('L', 4, JJt, 4, INFO)
    CALL dpotrs('L', 4, 1, JJt, 4, vetG, 4, INFO)

    IF (INFO == 0) THEN
      ! aplica a correcao
      DO a = 1, 4
        grads(a,:) = vetG(a)*grads(a,:)
      END DO
      gradsT = transpose(grads)

      DO a = 1, 6*N
        vetorCorrecao(a) = sum(gradsT(a,:))
      END DO

      ! Ve se esta na margem de erro
      IF (ABS(MINVAL(vetorCorrecao)) .le. 1.0E-20_pf .AND. ABS(MAXVAL(vetorCorrecao)) .le. 1.0E-20_pf) THEN
        ! Se estiver, manda embora porque a correcao sera desnecessaria
        ! WRITE(*,*) "dentro da margem aceitavel"
        EXIT
      END IF

      ! aplica a correcao
      DO a = 1, N
        posicoes(a,1) = posicoes(a,1) + vetorCorrecao(6*a-5)
        posicoes(a,2) = posicoes(a,2) + vetorCorrecao(6*a-4)
        posicoes(a,3) = posicoes(a,3) + vetorCorrecao(6*a-3)
        momentos(a,1) = momentos(a,1) + vetorCorrecao(6*a-2)
        momentos(a,2) = momentos(a,2) + vetorCorrecao(6*a-1)
        momentos(a,3) = momentos(a,3) + vetorCorrecao(6*a-0)
      END DO

      corrigiu = .TRUE.

    ELSE
      ! Caso contrario, nao tem como corrigir e so sai
      corrigiu = .FALSE.
      WRITE(*,*) "problema na matriz"
      EXIT loop
    
    ENDIF

  END DO loop
END SUBROUTINE corrigir

! ************************************************************
!! Aplicacao
!
! Objetivos:
!   Aplicacao da correcao numerica da energia total apenas.
!
! Modificado:
!   10 de novembro de 2024
!
! Autoria:
!   oap
!
SUBROUTINE corrigir_apenas_energia (corme, cormnt, G, massas, posicoes, momentos, corrigiu, H0, J, H_atual)
  IMPLICIT NONE
  REAL(pf), INTENT(IN)    :: corme  ! CORrecao Margem Erro
  INTEGER,  INTENT(IN)    :: cormnt ! CORrecao Max Num Tentativas
  REAL(pf), INTENT(IN)    :: G, massas(:)
  REAL(pf), INTENT(INOUT) :: posicoes(:,:), momentos(:,:)
  LOGICAL,  INTENT(INOUT) :: corrigiu
  INTEGER                 :: N, a, b, contador = 0
  REAL(pf)                :: H0, J(3), gradE2, alpha, H_atual
  REAL(pf), ALLOCATABLE   :: gradE(:)

  N = SIZE(massas)
  contador = 0

  aLlocate(gradE(1:6*N))
  
  ! enquanto nao tiver aceitado a correcao, roda
  loop: DO WHILE (.TRUE.)

    IF (contador >= cormnt) THEN
      ! WRITE(*,*) "max num tentativas atingido"
      EXIT loop
    ENDIF

    ! Contador de correcoes
    contador = contador + 1

    gradE=gradiente_energia(G, massas, posicoes, momentos, N)
    gradE2 = DOT_PRODUCT(gradE,gradE)

    alpha = (H0 - H_atual) / gradE2

    ! aplica a correcao
    gradE = gradE * alpha

    ! Ve se esta na margem de erro
    ! IF (MAXVAL(gradE) .le. 1E-8) THEN
    !   ! Se estiver, manda embora porque a correcao sera desnecessaria
    !   ! WRITE(*,*) "dentro da margem aceitavel: ", NORM2(gradE)
    !   EXIT
    ! END IF

    ! aplica a correcao
    DO a = 1, N
      posicoes(a,1) = posicoes(a,1) + gradE(6*a-5)
      posicoes(a,2) = posicoes(a,2) + gradE(6*a-4)
      posicoes(a,3) = posicoes(a,3) + gradE(6*a-3)
      momentos(a,1) = momentos(a,1) + gradE(6*a-2)
      momentos(a,2) = momentos(a,2) + gradE(6*a-1)
      momentos(a,3) = momentos(a,3) + gradE(6*a-0)
    END DO

    corrigiu = .TRUE.

  END DO loop
END SUBROUTINE corrigir_apenas_energia

! ************************************************************
!! Vetor de estado
!
! Objetivos:
!   Calcula o vetor de estados (i.e., vetor com cada integral
!   primeira)
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
FUNCTION Gx (G, massas, posicoes, momentos, N)
  IMPLICIT NONE
  INTEGER  :: N
  REAL(pf) :: G, massas(N), posicoes(N,3), momentos(N,3), Gx(4), J(3), Rcm(3), P(3)

  ! energia total
  Gx(1) = energia_total(G, massas, posicoes, momentos)

  ! momento angular
  J = momento_angular_total(posicoes, momentos)
  Gx(2) = J(1)
  Gx(3) = J(2)
  Gx(4) = J(3)

  Gx = - Gx

END FUNCTION Gx

! ************************************************************
!! Matriz normal
!
! Objetivos:
!   Calcula a matriz produto dos jacobianos (uma matriz normal)
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
FUNCTION matriz_normal (G, massas, posicoes, momentos, N, grads)

  IMPLICIT NONE
  INTEGER  :: N, gi, gj
  REAL(pf) :: G, massas(N), posicoes(N,3), momentos(N,3)
  REAL(pf),INTENT(INOUT) :: grads(4, 6*N)
  REAL(pf) :: matriz_normal(4,4)

  grads(1,:)=gradiente_energia(G, massas, posicoes, momentos, N)
  grads(2,:)=gradiente_angularX(posicoes, momentos, N)
  grads(3,:)=gradiente_angularY(posicoes, momentos, N)
  grads(4,:)=gradiente_angularZ(posicoes, momentos, N)

  DO gi = 1, 4
    DO gj = 1, 4
      matriz_normal(gi, gj) = DOT_PRODUCT(grads(gi,:), grads(gj,:))
    END DO
  END DO

END FUNCTION matriz_normal

! ************************************************************
!! Gradiente: energia total
!
! Objetivos:
!   Vetores de derivadas parciais (vetor gradiente) da energia
!   total.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
FUNCTION gradiente_energia (G, massas, Rs, Ps, N)

  IMPLICIT NONE
  INTEGER  :: N, a, b
  REAL(pf) :: distancia3, distancia(3)
  REAL(pf) :: G, massas(N), Rs(N,3), Ps(N,3), gradiente_energia(1:6*N)

  gradiente_energia(:) = 0.0_pf

  DO a = 1, N
    gradiente_energia(6*a-2) = Ps(a,1)/massas(a)
    gradiente_energia(6*a-1) = Ps(a,2)/massas(a)
    gradiente_energia(6*a-0) = Ps(a,3)/massas(a)
    
    DO b = 1, N
      IF (b /= a) THEN
        distancia = Rs(b,:)-Rs(a,:)
        distancia3 = norm2(distancia)**3
        distancia = (massas(b)/distancia3) * distancia

        gradiente_energia(6*a-5) = gradiente_energia(6*a-5)+distancia(1)
        gradiente_energia(6*a-4) = gradiente_energia(6*a-4)+distancia(2)
        gradiente_energia(6*a-3) = gradiente_energia(6*a-3)+distancia(3)
      ENDIF
    END DO

    gradiente_energia(6*a-5) = gradiente_energia(6*a-5) * (-G*massas(a))
    gradiente_energia(6*a-4) = gradiente_energia(6*a-4) * (-G*massas(a))
    gradiente_energia(6*a-3) = gradiente_energia(6*a-3) * (-G*massas(a))
  END DO
  
END FUNCTION gradiente_energia

! ************************************************************
!! Gradiente: momento angular (x)
!
! Objetivos:
!   Vetores de derivadas parciais (vetor gradiente) do momento
!   angular total (x).
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
FUNCTION gradiente_angularX (posicoes, momentos, N)
  IMPLICIT NONE
  INTEGER  :: N, a
  REAL(pf) :: gradiente_angularX(6*N), posicoes(N,3), momentos(N,3)

  DO a = 1, N
    gradiente_angularX(6*a-5) = 0.0_pf
    gradiente_angularX(6*a-4) = momentos(a,3)
    gradiente_angularX(6*a-3) = - momentos(a,2)
    gradiente_angularX(6*a-2) = 0.0_pf
    gradiente_angularX(6*a-1) = - posicoes(a,3)
    gradiente_angularX(6*a-0) = posicoes(a,2)
  END DO

END FUNCTION gradiente_angularX

! ************************************************************
!! Gradiente: momento angular (y)
!
! Objetivos:
!   Vetores de derivadas parciais (vetor gradiente) do momento
!   angular total (y).
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
FUNCTION gradiente_angularY (posicoes, momentos, N)
  IMPLICIT NONE
  INTEGER  :: N, a
  REAL(pf) :: posicoes(N,3), momentos(N,3), gradiente_angularY(6*N)

  DO a = 1, N
    gradiente_angularY(6*a-5) = - momentos(a,3)
    gradiente_angularY(6*a-4) = 0.0_pf
    gradiente_angularY(6*a-3) = momentos(a,1)
    gradiente_angularY(6*a-2) = posicoes(a,3)
    gradiente_angularY(6*a-1) = 0.0_pf
    gradiente_angularY(6*a-0) = - posicoes(a,1)
  END DO

END FUNCTION gradiente_angularY

! ************************************************************
!! Gradiente: momento angular (z)
!
! Objetivos:
!   Vetores de derivadas parciais (vetor gradiente) do momento
!   angular total (z).
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
FUNCTION gradiente_angularZ (posicoes, momentos, N)
  IMPLICIT NONE
  INTEGER  :: N, a
  REAL(pf) :: posicoes(N,3), momentos(N,3), gradiente_angularZ(6*N)

  DO a = 1, N
    gradiente_angularZ(6*a-5) = momentos(a,2)
    gradiente_angularZ(6*a-4) = - momentos(a,1)
    gradiente_angularZ(6*a-3) = 0.0_pf
    gradiente_angularZ(6*a-2) = - posicoes(a,2)
    gradiente_angularZ(6*a-1) = posicoes(a,1)
    gradiente_angularZ(6*a-0) = 0.0_pf
  END DO

END FUNCTION gradiente_angularZ

! ************************************************************
!! Gradiente: momento linear (x)
!
! Objetivos:
!   Vetores de derivadas parciais (vetor gradiente) do momento
!   linear total (x).
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
FUNCTION gradiente_linearX (N)
  IMPLICIT NONE
  INTEGER :: N, a
  REAL(pf) :: gradiente_linearX(6*N)

  DO a = 1, N
    gradiente_linearX(6*a-5) = 0.0_pf
    gradiente_linearX(6*a-4) = 0.0_pf
    gradiente_linearX(6*a-3) = 0.0_pf
    gradiente_linearX(6*a-2) = 1.0_pf
    gradiente_linearX(6*a-1) = 0.0_pf
    gradiente_linearX(6*a-0) = 0.0_pf
  END DO

END FUNCTION gradiente_linearX

! ************************************************************
!! Gradiente: momento linear (y)
!
! Objetivos:
!   Vetores de derivadas parciais (vetor gradiente) do momento
!   linear total (y).
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
FUNCTION gradiente_linearY (N)
  IMPLICIT NONE
  INTEGER :: N, a
  REAL(pf) :: gradiente_linearY(6*N)

  DO a = 1, N
    gradiente_linearY(6*a-5) = 0.0_pf
    gradiente_linearY(6*a-4) = 0.0_pf
    gradiente_linearY(6*a-3) = 0.0_pf
    gradiente_linearY(6*a-2) = 0.0_pf
    gradiente_linearY(6*a-1) = 1.0_pf
    gradiente_linearY(6*a-0) = 0.0_pf
  END DO

END FUNCTION gradiente_linearY

! ************************************************************
!! Gradiente: momento linear (z)
!
! Objetivos:
!   Vetores de derivadas parciais (vetor gradiente) do momento
!   linear total (z).
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
FUNCTION gradiente_linearZ (N)
  IMPLICIT NONE
  INTEGER :: N, a
  REAL(pf) :: gradiente_linearZ(6*N)

  DO a = 1, N
    gradiente_linearZ(6*a-5) = 0.0_pf
    gradiente_linearZ(6*a-4) = 0.0_pf
    gradiente_linearZ(6*a-3) = 0.0_pf
    gradiente_linearZ(6*a-2) = 0.0_pf
    gradiente_linearZ(6*a-1) = 0.0_pf
    gradiente_linearZ(6*a-0) = 1.0_pf
  END DO

END FUNCTION gradiente_linearZ

! ************************************************************
!! Gradiente: centro de massas (x)
!
! Objetivos:
!   Vetores de derivadas parciais (vetor gradiente) do centro
!   de massas (x).
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
FUNCTION gradiente_centroMassasX (massas, N)
  IMPLICIT NONE
  INTEGER  :: N, a
  REAL(pf) :: M
  REAL(pf) :: massas(N), gradiente_centroMassasX(6*N)

  M = sum(massas)

  DO a = 1, N
    gradiente_centroMassasX(6*a-5) = massas(a)/M
    gradiente_centroMassasX(6*a-4) = 0.0_pf
    gradiente_centroMassasX(6*a-3) = 0.0_pf
    gradiente_centroMassasX(6*a-2) = 0.0_pf
    gradiente_centroMassasX(6*a-1) = 0.0_pf
    gradiente_centroMassasX(6*a-0) = 0.0_pf
  END DO

END FUNCTION gradiente_centroMassasX

! ************************************************************
!! Gradiente: centro de massas (y)
!
! Objetivos:
!   Vetores de derivadas parciais (vetor gradiente) do centro
!   de massas (y).
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
FUNCTION gradiente_centroMassasY (massas, N)
  IMPLICIT NONE
  INTEGER  :: N, a
  REAL(pf) :: M
  REAL(pf) :: massas(N), gradiente_centroMassasY(6*N)

  M = sum(massas)

  DO a = 1, N
    gradiente_centroMassasY(6*a-5) = 0.0_pf
    gradiente_centroMassasY(6*a-4) = massas(a)/M
    gradiente_centroMassasY(6*a-3) = 0.0_pf
    gradiente_centroMassasY(6*a-2) = 0.0_pf
    gradiente_centroMassasY(6*a-1) = 0.0_pf
    gradiente_centroMassasY(6*a-0) = 0.0_pf
  END DO

END FUNCTION gradiente_centroMassasY

! ************************************************************
!! Gradiente: centro de massas (z)
!
! Objetivos:
!   Vetores de derivadas parciais (vetor gradiente) do centro
!   de massas (z).
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
FUNCTION gradiente_centroMassasZ (massas, N)
  IMPLICIT NONE
  INTEGER  :: N, a
  REAL(pf) :: M
  REAL(pf) :: massas(N), gradiente_centroMassasZ(6*N)

  M = SUM(massas)

  DO a = 1, N
    gradiente_centroMassasZ(6*a-5) = 0.0_pf
    gradiente_centroMassasZ(6*a-4) = 0.0_pf
    gradiente_centroMassasZ(6*a-3) = massas(a)/M
    gradiente_centroMassasZ(6*a-2) = 0.0_pf
    gradiente_centroMassasZ(6*a-1) = 0.0_pf
    gradiente_centroMassasZ(6*a-0) = 0.0_pf
  END DO

END FUNCTION gradiente_centroMassasZ

END MODULE correcao