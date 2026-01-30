! ************************************************************
!! COLISAO DE PARTICULAS
!
! Objetivos:
!   Aplicacao das colisoes perfeitamente elasticas entre duas
!   particulas.
!
! Modificado:
!   29 de janeiro de 2026
!
! Autoria:
!   oap
!  
MODULE colisao
  USE tipos
  USE octree
  IMPLICIT NONE
  PRIVATE
  PUBLIC verificar_e_colidir

  !> Arvore pre-instanciada para ser utilizada sem realocacao
  TYPE(arvore_octo), ALLOCATABLE :: arvore
CONTAINS

! ************************************************************
!! Verifica e colide
!
! Objetivos:
!   Determina e aplica colisoes perfeitamente elasticas entre
!   particulas.
!
! Modificado:
!   27 de janeiro de 2026
!
! Autoria:
!   oap
! 
SUBROUTINE verificar_e_colidir (m, R, P, paralelo, raios, modo)
  IMPLICIT NONE
  REAL(pf) :: m(:), R(:,:), P(:,:) ! maximo de aproximacao
  LOGICAL  :: paralelo
  REAL(pf) :: raios(size(m))
  CHARACTER(LEN=*), INTENT(IN) :: modo

  ! Aplica colisoes conforme o modo
  IF (TRIM(modo) == 'T' .OR. TRIM(modo) == 'direto') THEN
    CALL verificar_e_colidir_direto(m, R, P, paralelo, raios)
  ELSE IF (TRIM(modo) == 'octree') THEN
    CALL verificar_e_colidir_octree(m, R, P, paralelo, raios)
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
!   24 de julho de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE verificar_e_colidir_octree (m, R, P, paralelo, raios)
  IMPLICIT NONE
  REAL(pf) :: m(:), R(:,:), P(:,:)
  INTEGER :: a, b, i, colididos, colisoes(size(m)-1)
  LOGICAL :: paralelo, colidiram(size(m),size(m))
  REAL(pf) :: raios(size(m))

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
!   29 de janeiro de 2026
!
! Autoria:
!   oap
! 
SUBROUTINE verificar_e_colidir_direto (m, R, P, paralelo, raios)
  IMPLICIT NONE
  REAL(pf) :: m(:), R(:,:), P(:,:)
  INTEGER :: a, b
  LOGICAL :: paralelo, colidiram(INT(size(m)*(size(m)-1)/2))
  REAL(pf) :: raios(size(m)), dist
  INTEGER :: indice

  colidiram=.FALSE.
  DO a = 2, SIZE(m)
    DO b = 1, a-1
      indice = INT((a-1)*(a-2)/2 + (b-1)) + 1
      
      IF (a == b .OR. colidiram(indice)) THEN
        CYCLE
      ENDIF
      
      dist = (R(b,1) - R(a,1))**2 + (R(b,2) - R(a,2))**2 + (R(b,3) - R(a,3))**2
      IF (dist <= (raios(a) + raios(b))**2) THEN
        IF (DOT_PRODUCT(R(b,:)-R(a,:), P(b,:)-P(a,:)) < 0) THEN
          colidiram(indice) = .TRUE.
          CALL colidir (m(a), R(a,:), P(a,:), m(b), R(b,:), P(b,:))
        ENDIF
      ENDIF 
    END DO
  END DO
END SUBROUTINE verificar_e_colidir_direto

! ************************************************************
!! Colisao
!
! Objetivos:
!   Faz o calculo da colisao perfeitamente elastica entre duas
!   particulas A e B.
!
! Modificado:
!   03 de agosto de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE colidir (ma, Ra, Pa, mb, Rb, Pb)
  IMPLICIT NONE
  REAL(pf) :: ma, mb, Ra(3), Pa(3), Rb(3), Pb(3)
  REAL(pf) :: Normal(3), norma2, u1, u2, k

  ! vetor normal e normal unitario
  Normal = Rb - Ra
  norma2 = DOT_PRODUCT(Normal, Normal)

  ! calcula a componente normal
  u1 = DOT_PRODUCT(Pa, Normal) / ma
  u2 = DOT_PRODUCT(Pb, Normal) / mb

  ! agora calcula o componente de angulo
  k = 2.0_pf * (u2 - u1) * ma * mb / ((ma + mb) * norma2)

  ! por fim, aplica a colisao
  Pa = Pa + k * Normal
  Pb = Pb - k * Normal 

END SUBROUTINE colidir

END MODULE colisao