! ************************************************************
!! COLISAO DE PARTICULAS
!
! Objetivos:
!   Aplicacao das colisoes perfeitamente elasticas entre duas
!   particulas.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
! Funcoes relativas a colisao de particulas
! 
! Funcoes
! 
! = SUBROUTINE verificar_e_colidir (m, R, P)
! Passa por cada par de particulas avaliando as distancias e
! a taxa de variacao desta, verificando se trata-se de caso
! de colisao ou nao.
! 
! = subourtine colidir (ma, Ra, Pa, mb, Rb, Pb)
! Faz o calculo da colisao perfeitamente elastica entre
! duas particulas A e B dadas.
! 
MODULE colisao
  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  IMPLICIT NONE
  PRIVATE
  PUBLIC verificar_e_colidir, colidir

CONTAINS

! ************************************************************
!! Verifica e colide
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
!   11 de novembro de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE verificar_e_colidir (m, R, P, colmd, paralelo)
  IMPLICIT NONE
  REAL(pf) :: m(:), R(:,:), P(:,:), colmd ! maximo de aproximacao
  INTEGER :: a, b
  LOGICAL :: paralelo
  REAL(pf) :: m13a, m13b, fator

  IF (paralelo) THEN
    !$OMP PARALLEL PRIVATE(a, b, m13a, m13b, fator)
      !$OMP DO
      DO a = 2, SIZE(m)
        m13a = m(a) ** (1.0_pf/3.0_pf)
        DO b = 1, a-1
          m13b = m(b) ** (1.0_pf/3.0_pf)
          fator = norm2(R(b,:)-R(a,:)) / ABS(m13a + m13b)
          IF (fator <= colmd) THEN
            IF (DOT_PRODUCT(R(b,:) - R(a,:), P(b,:)-P(a,:)) < 0) THEN
              CALL colidir (m(a), R(a,:), P(a,:), m(b), R(b,:), P(b,:))
            ENDIF
          ENDIF 
        END DO
      END DO
      !$OMP END DO
    !$OMP END PARALLEL
  ELSE
    DO a = 2, SIZE(m)
      m13a = m(a) ** (1.0_pf/3.0_pf)
      DO b = 1, a-1
        m13b = m(b) ** (1.0_pf/3.0_pf)
        fator = norm2(R(b,:)-R(a,:)) / ABS(m13a + m13b)
        IF (fator <= colmd) THEN
          ! CALL colidir (m(a), R(a,:), P(a,:), m(b), R(b,:), P(b,:))

          IF (DOT_PRODUCT(R(b,:) - R(a,:), P(b,:)-P(a,:)) < 0) THEN
            CALL colidir (m(a), R(a,:), P(a,:), m(b), R(b,:), P(b,:))
          ENDIF
        ENDIF 
      END DO
    END DO
  ENDIF
END SUBROUTINE verificar_e_colidir

! ************************************************************
!! Colisao
!
! Objetivos:
!   Faz o calculo da colisao perfeitamente elastica entre duas
!   particulas A e B.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE colidir (ma, Ra, Pa, mb, Rb, Pb)
  IMPLICIT NONE
  REAL(pf) :: ma, mb, Ra(3), Pa(3), Rb(3), Pb(3)
  REAL(pf) :: ua(3), ub(3), Normal(3), Normal_(3), u1, u2, ua_p(3), ub_p(3)
  REAL(pf) :: S(3), T(3)

  ! separa as velocidades
  ua = Pa/ma
  ub = Pb/mb

  ! vetor normal e normal unitario
  Normal = Rb - Ra
  Normal_ = Normal/norm2(Normal)

  ! calcula a componente tangente
  u1 = DOT_PRODUCT(ua, Normal_)
  u2 = DOT_PRODUCT(ub, Normal_)
  
  ! Para calcular o plano tangente, tome T = (-n2_, n1_, 0)
  T(1) = - Normal_(2)
  T(2) =   Normal_(1)
  T(3) =   0.0_pf
  ! e S = N_ x T
  S(1) = - Normal_(1) * Normal_(3)
  S(2) = - Normal_(2) * Normal_(3)
  S(3) =   Normal_(1)**2 + Normal_(2)**2
  
  ! calcula as componentes do plano
  ua_p = DOT_PRODUCT(ua, S) * S + DOT_PRODUCT(ua, T) * T
  ub_p = DOT_PRODUCT(ub, S) * S + DOT_PRODUCT(ub, T) * T

  ! obtem as novas velocidades
  Pa = ma * (ua_p + (u1*(ma-mb)+2*mb*u2)/(ma+mb) * Normal_)
  Pb = mb * (ub_p + (u2*(mb-ma)+2*ma*u1)/(ma+mb) * Normal_)

END SUBROUTINE colidir

END MODULE colisao