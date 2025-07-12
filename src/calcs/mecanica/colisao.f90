! ************************************************************
!! COLISAO DE PARTICULAS
!
! Objetivos:
!   Aplicacao das colisoes perfeitamente elasticas entre duas
!   particulas.
!
! Modificado:
!   16 de marco de 2025
!
! Autoria:
!   oap
!  
MODULE colisao
  USE tipos
  USE octree
  USE mecanica
  IMPLICIT NONE
  PRIVATE
  PUBLIC verificar_e_colidir

CONTAINS

! ************************************************************
!! Verifica e colide
!
! Objetivos:
!   Determina e aplica colisoes perfeitamente elasticas entre
!   particulas.
!
! Modificado:
!   16 de marco de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE verificar_e_colidir (m, R, P, colmd, paralelo, raios, arvore, modo, dists)
  IMPLICIT NONE
  TYPE(arvore_octo), ALLOCATABLE :: arvore
  REAL(pf) :: m(:), R(:,:), P(:,:), colmd ! maximo de aproximacao
  LOGICAL  :: paralelo
  REAL(pf) :: raios(size(m))
  CHARACTER(LEN=*), INTENT(IN) :: modo
  REAL(pf) :: dists(:)

  ! Aplica colisoes conforme o modo
  IF (TRIM(modo) == 'T' .OR. TRIM(modo) == 'direto') THEN
    CALL verificar_e_colidir_direto(m, R, P, colmd, paralelo, raios, dists)
  ELSE IF (TRIM(modo) == 'octree') THEN
    CALL verificar_e_colidir_octree(m, R, P, paralelo, raios, arvore)
  ELSE
    WRITE (*,*) "ATENCAO! MODO DE COLISAO NAO IDENTIFICADO: ", modo
    CALL ABORT()
  ENDIF
END SUBROUTINE verificar_e_colidir

! ************************************************************
!! Verifica e colide atraves de uma octree
!
! Objetivos:
!   Utilizando uma arvore octenaria, verifica colisoes entre
!   corpos.
!
! Modificado:
!   16 de marco de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE verificar_e_colidir_octree (m, R, P, paralelo, raios, arvore)
  IMPLICIT NONE
  REAL(pf) :: m(:), R(:,:), P(:,:)
  INTEGER :: a, b, i, colididos, colisoes(size(m)-1)
  LOGICAL :: paralelo, colidiram(size(m),size(m))
  REAL(pf) :: raios(size(m))
  TYPE(arvore_octo), ALLOCATABLE :: arvore

  colidiram = .FALSE.

  ! Gera e ordena as arvores binarias
  CALL gerar_octree(arvore, size(m), m, R)

  ! Agora percorre os corpos para detectar colisoes
  DO a=1, size(m)
    colididos = 0
    colisoes = 0
    CALL verificar_colisao_octree(arvore, a, R, raios, colisoes, colididos)

    ! Se tiver tido alguma colisao
    IF (colididos > 0) THEN
      ! DO i=1, colididos
      DO i=1, size(M)
        IF (ANY(colisoes==i))THEN
          b = i
        ELSE
          CYCLE
        ENDIF
        ! b = colisoes(i)
        IF (colidiram(a,b) .OR. colidiram(b,a)) THEN
          CYCLE
        ENDIF

        IF (DOT_PRODUCT(R(b,:) - R(a,:), P(b,:)-P(a,:)) < 0) THEN
          colidiram(a,b) = .TRUE.
          colidiram(b,a) = .TRUE.
          ! WRITE(*,*) 'colidiram octree:', a,b, NORM2(R(a,:)-R(b,:)) ! debug
          CALL colidir(m(a),R(a,:),P(a,:),m(b),R(b,:),P(b,:))
        ENDIF
      END DO
    ENDIF
  END DO

END SUBROUTINE verificar_e_colidir_octree

! ************************************************************
!! Verifica e colide direto
!
! Objetivos:
!   Passa por cada par de particulas avaliando as distancias e
!   taxa de variacao desta, verificando se trata-se de caso de
!   colisao ou nao.
!   Condicoes para colisao:
!     1. fator <= colmd
!     2. <rb - ra, pb - pa> < 0
!
! Modificado:
!   16 de marco de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE verificar_e_colidir_direto (m, R, P, colmd, paralelo, raios, dists)
  IMPLICIT NONE
  REAL(pf) :: m(:), R(:,:), P(:,:), colmd ! maximo de aproximacao
  INTEGER :: a, b
  LOGICAL :: paralelo, colidiram(size(m),size(m))
  REAL(pf) :: m13a, m13b, fator, raios(size(m)), dist
  
  REAL(pf) :: dists(:)
  INTEGER :: indice

  colidiram=.FALSE.

  IF (paralelo) THEN
    !$OMP PARALLEL PRIVATE(a, b, m13a, m13b, fator)
    !$OMP DO
    DO a = 1, SIZE(m)
      m13a = m(a) ** (1.0_pf/3.0_pf)
      DO b = 1, SIZE(m)
        IF (a == b .OR. colidiram(a,b) .OR. colidiram(b,a)) THEN
          CYCLE
        ENDIF
        m13b = m(b) ** (1.0_pf/3.0_pf)
        indice = INT((a-1)*(a-2)/2 + (b-1)) + 1
        fator = dists(indice) / ABS(m13a + m13b)
        IF (fator <= colmd) THEN
          IF (DOT_PRODUCT(R(b,:) - R(a,:), P(b,:)-P(a,:)) < 0) THEN
            colidiram(a,b) = .TRUE.
            colidiram(b,a) = .TRUE.
            CALL colidir (m(a), R(a,:), P(a,:), m(b), R(b,:), P(b,:))
          ENDIF
        ENDIF 
      END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
  ELSE
    DO a = 2, SIZE(m)
      DO b = 1, a-1
        IF (a == b .OR. colidiram(a,b) .OR. colidiram(b,a)) THEN
          CYCLE
        ENDIF
        
        indice = INT((a-1)*(a-2)/2 + (b-1)) + 1
        dist = dists(indice)
        IF (dist <= raios(a) + raios(b)) THEN
          IF (DOT_PRODUCT(R(b,:)-R(a,:), P(b,:)-P(a,:)) < 0) THEN
            colidiram(a,b) = .TRUE.
            colidiram(b,a) = .TRUE.
            CALL colidir (m(a), R(a,:), P(a,:), m(b), R(b,:), P(b,:))
            ! WRITE(*,*) 'colidiram normal:', a,b, NORM2(R(b,:)-R(a,:)) ! debug
          ENDIF
        ENDIF 
      END DO
    END DO
  ENDIF
END SUBROUTINE verificar_e_colidir_direto

! ************************************************************
!! Colisao
!
! Objetivos:
!   Faz o calculo da colisao perfeitamente elastica entre duas
!   particulas A e B.
!
! Modificado:
!   08 de abril de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE colidir (ma, Ra, Pa, mb, Rb, Pb)
  IMPLICIT NONE
  REAL(pf) :: ma, mb, Ra(3), Pa(3), Rb(3), Pb(3)
  REAL(pf) :: Normal(3), Normal_(3), u1, u2, k

  ! vetor normal e normal unitario
  Normal = Rb - Ra
  Normal_ = Normal/norm2(Normal)

  ! calcula a componente normal
  u1 = DOT_PRODUCT(Pa, Normal_) / ma
  u2 = DOT_PRODUCT(Pb, Normal_) / mb

  ! agora calcula o componente de angulo
  k = 2.0_pf * (u2 - u1) * ma * mb / (ma + mb)

  ! por fim, aplica a colisao
  Pa = Pa + k * Normal_
  Pb = Pb - k * Normal_

END SUBROUTINE colidir

END MODULE colisao