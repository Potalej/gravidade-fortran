! *****************************************************************
!! ALEATORIO
!
! Objetivos:
!   Este arquivo contem helpers de aleatoriedade. O modulo contem
!   subrotinas de geracao de valores iniciais com diferentes
!   distribuicoes de probabilidade.
! 
! Modificado:
!   15 de janeiro de 2025
! 
! Autoria:
!   oap
! 
MODULE aleatorio
  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  IMPLICIT NONE
CONTAINS

FUNCTION regiao_cubo (vetor3d, raio)

  REAL(pf), INTENT(IN), DIMENSION(3) :: vetor3d
  REAL(pf), INTENT(IN) :: raio
  LOGICAL :: regiao_cubo

  regiao_cubo = (ABS(vetor3d(1)) <= raio .AND. ABS(vetor3d(2)) <= raio .AND. ABS(vetor3d(3)) <= raio)

END FUNCTION regiao_cubo

FUNCTION regiao_esfera (vetor3d, raio)

  REAL(pf), INTENT(IN), DIMENSION(3) :: vetor3d
  REAL(pf), INTENT(IN) :: raio
  LOGICAL :: regiao_esfera

  regiao_esfera = (NORM2(vetor3d) <= raio)

END FUNCTION regiao_esfera

FUNCTION regiao_cobrinha (vetor3d, raio)

  REAL(pf), INTENT(IN), DIMENSION(3) :: vetor3d
  REAL(pf), INTENT(IN) :: raio
  LOGICAL :: regiao_cobrinha

  regiao_cobrinha = (ABS(vetor3d(1)) <= 2 .AND. ABS(vetor3d(3)) <= 2)
  IF (regiao_cobrinha) THEN
    regiao_cobrinha = ((vetor3d(2) <= COS(2*vetor3d(1))+1.0_pf) .AND. (vetor3d(2) >= COS(2*vetor3d(1))))
  ENDIF
  
END FUNCTION regiao_cobrinha

FUNCTION testar_regiao (vetor3d, raio, regiao)

  REAL(pf), INTENT(IN), DIMENSION(3) :: vetor3d
  REAL(pf), INTENT(IN) :: raio
  CHARACTER(LEN=*),INTENT(IN) :: regiao
  LOGICAL :: testar_regiao

  SELECT CASE (TRIM(regiao))
    CASE ('cubo');     testar_regiao = regiao_cubo(vetor3d, raio)
    CASE ('esfera');   testar_regiao = regiao_esfera(vetor3d, raio)
    CASE ('cobrinha'); testar_regiao = regiao_cobrinha(vetor3d, raio)
  END SELECT

END FUNCTION testar_regiao


SUBROUTINE condiciona_esfera (vetor, N, raio)

  REAL(pf), INTENT(INOUT), DIMENSION(N,3) :: vetor
  INTEGER, INTENT(IN)  :: N
  REAL(pf), INTENT(IN) :: raio
  REAL(pf) :: x,y,z,x_,y_,z_
  INTEGER :: a

  DO a = 1,N
    x = vetor(a,1)
    y = vetor(a,2)
    z = vetor(a,3)
    x_ = x * SQRT(1.0_pf - 0.5_pf * y * y - 0.5_pf * z * z + y * y * z * z / 3.0_pf)
    y_ = y * SQRT(1.0_pf - 0.5_pf * z * z - 0.5_pf * x * x + z * z * x * x / 3.0_pf)
    z_ = z * SQRT(1.0_pf - 0.5_pf * x * x - 0.5_pf * y * y + x * x * y * y / 3.0_pf)
    vetor(a,:) = (/x_,y_,z_/)
  END DO

END SUBROUTINE condiciona_esfera

! ************************************************************
!! Gera valores iniciais com distribuicao uniforme
!
! Objetivos:
!   Gera valores iniciais com distribuicao uniforme atraves
!   da subrotina padrao do Fortran.
!
! Modificado:
!   02 de maio de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE uniforme (vetor, N, distmin, vmin, vmax, regiao, raio)

  REAL(pf), INTENT(INOUT), DIMENSION(N,3) :: vetor
  CHARACTER(LEN=*), INTENT(IN) :: regiao
  INTEGER,  INTENT(IN) :: N
  REAL(pf), INTENT(IN) :: raio, distmin, vmin, vmax
  REAL(pf) :: u(3), ajuste(3)
  INTEGER :: i

  ajuste(:) = vmin
  CALL RANDOM_SEED()
  DO i=1, N
    DO
      ! Sorteia e ajusta no intervalo
      CALL RANDOM_NUMBER(u)
      u = raio * (u * (vmax - vmin) + ajuste)

      IF (testar_regiao(u, raio, regiao)) THEN
        IF (i >= 1 .AND. distmin > 0) THEN
          vetor(i,:) = u
          IF (testar_distancias(vetor, N, distmin, i)) THEN
            EXIT
          ENDIF
        ELSE
          vetor(i,:) = u
          EXIT
        ENDIF
      END IF 
    END DO
  END DO

END SUBROUTINE uniforme

! ************************************************************
!! Gera valores iniciais com distribuicao normal(0,1)
!
! Objetivos:
!   Gera valores iniciais com distribuicao normal(0,1) atraves
!   do metodo de Box-Muller. A terceira coordenada eh gerada
!   aplicando novamente o metodo de Box-Muller e extraindo a
!   primeira componente.
!
! Modificado:
!   15 de janeiro de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE normal (vetor, N, distmin, vmin, vmax, regiao, raio)

  INTEGER,  INTENT(IN) :: N
  REAL(pf), INTENT(INOUT), DIMENSION(N,3) :: vetor
  CHARACTER(LEN=*), INTENT(IN)     :: regiao
  REAL(pf), INTENT(IN)     :: vmin, vmax, raio, distmin
  REAL(pf), DIMENSION(N,3) :: ajuste
  REAL(pf) :: x, y, z
  REAL(pf) :: u1, u2, r, theta, max_dist
  REAL(pf) :: PI = 3.14159265358979_pf
  INTEGER  :: a
  ajuste(:,:) = vmin

  CALL RANDOM_SEED()

  max_dist = 0.0_pf

  ! transformacao de Box-Muller
  DO a = 1, N
    DO
      CALL RANDOM_NUMBER(u1)
      CALL RANDOM_NUMBER(u2)
      r = SQRT(-2.0_pf * LOG(u1))
      theta = 2.0_pf * PI * u2
      x = r * COS(theta)
      y = r * SIN(theta)

      CALL RANDOM_NUMBER(u1)
      CALL RANDOM_NUMBER(u2)
      
      r = SQRT(-2.0_pf * LOG(u1))
      theta = 2.0_pf * PI * u2
      z = r * COS(theta)

      IF (testar_regiao((/x,y,z/), raio, regiao)) THEN
        IF (a >= 1 .AND. distmin > 0) THEN
          vetor(a,:) = (/x,y,z/)
          IF (testar_distancias(vetor, N, distmin, a)) THEN
            EXIT
          ENDIF
        ELSE
          vetor(a,:) = (/x,y,z/)
          EXIT
        ENDIF
      END IF 
    END DO
  END DO

END SUBROUTINE normal

! ************************************************************
!! Gera valores iniciais com distribuicao cauchy
!
! Objetivos:
!   Gera valores iniciais com distribuicao cauchy.
!
! Modificado:
!   17 de janeiro de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE cauchy (vetor, N, distmin, vmin, vmax, regiao, raio)

  INTEGER,  INTENT(IN) :: N
  REAL(pf), INTENT(INOUT), DIMENSION(N,3) :: vetor
  REAL(pf), INTENT(IN)     :: vmin, vmax, raio, distmin
  CHARACTER(LEN=*),INTENT(IN) :: regiao
  REAL(pf), DIMENSION(3)   :: u
  REAL(pf) :: PI = 3.14159265358979_pf
  INTEGER :: a
  REAL(pf), DIMENSION(3) :: ajuste
  ajuste(:) = vmin

  CALL RANDOM_SEED()

  DO a=1, N
    DO
      ! Sorteia e ajusta no intervalo
      CALL RANDOM_NUMBER(u)
      u = TAN(PI*u - PI*0.5_pf)

      IF (testar_regiao(u, raio, regiao)) THEN
        IF (a > 1 .AND. distmin > 0) THEN
          vetor(a,:) = u
          IF (testar_distancias(vetor, N, distmin, a)) THEN
            EXIT
          ENDIF
        ELSE
          vetor(a,:) = u
          EXIT
        ENDIF
      END IF 
    END DO
  END DO

END SUBROUTINE cauchy

! ************************************************************
!! Verifica se distancia minima eh atendida.
!
! Objetivos:
!   Percorre a lista de vetores ja gerados para ver se o ultimo
!   gerado esta a uma distancia minima exigida.
!
! Modificado:
!   03 de marco de 2025
!
! Autoria:
!   oap
! 
FUNCTION testar_distancias (vetor, N, distmin, ind)

  LOGICAL :: testar_distancias
  INTEGER, INTENT(IN) :: N, ind
  REAL(pf), INTENT(IN), DIMENSION(N,3) :: vetor
  REAL(pf), INTENT(IN) :: distmin
  INTEGER :: i

  testar_distancias = .TRUE.

  DO i = 1, ind-1
    IF (NORM2(vetor(ind,:) - vetor(i,:)) < distmin) THEN
      testar_distancias = .FALSE.
      EXIT
    ENDIF
  END DO

END FUNCTION

END MODULE aleatorio