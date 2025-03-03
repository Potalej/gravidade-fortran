MODULE funcoes_forca

    USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
    USE OMP_LIB
    IMPLICIT NONE
    PUBLIC

CONTAINS

FUNCTION forcas_par (m, R, G, N, dim, potsoft, potsoft2)
  IMPLICIT NONE
  INTEGER,                     INTENT(IN) :: N, dim
  REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
  REAL(pf), DIMENSION(N),      INTENT(IN) :: m
  REAL(pf),                    INTENT(IN) :: G, potsoft, potsoft2
  
  REAL(pf), DIMENSION(dim) :: Fab, Rab
  INTEGER :: a, b, tid
  REAL(pf) :: distancia, distancia_inv
  REAL(pf), DIMENSION(N, dim) :: forcas_par, forcas_local
  
  forcas_par = 0.0_pf

  !$OMP PARALLEL SHARED(forcas_par) PRIVATE(Fab, Rab, distancia, distancia_inv, a, b, tid)
  !$OMP DO
  DO a = 1, N
    DO b = 1, N
      IF (a==b) THEN
        CYCLE
      ENDIF
      ! distancia entre os corpos
      Rab = R(b,:) - R(a,:)
      distancia = norm2(Rab)
      IF (potsoft .NE. 0) THEN
        distancia = SQRT(distancia*distancia + potsoft2)
      ENDIF
      distancia_inv = 1.0_pf/distancia
      distancia_inv = distancia_inv**3

      ! forca entre os corpos a e b
      Fab = G * m(a) * m(b) * Rab * distancia_inv
      
      ! Adiciona na matriz      
      forcas_par(a,:) = forcas_par(a,:) + Fab
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL

END FUNCTION forcas_par

! ************************************************************
!! Matriz de forcas
!
! Objetivos:
!   Calcula a matriz de forcas a partir das posicoes. Todos os metodos
!   calculam as forcas do mesmo jeito.
!
! Modificado:
!   14 de setembro de 2024
!
! Autoria:
!   oap
! 
FUNCTION forcas_seq (m, R, G, N, dim, potsoft, potsoft2)
  IMPLICIT NONE
  INTEGER,                     INTENT(IN) :: N, dim
  REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
  REAL(pf), DIMENSION(N),      INTENT(IN) :: m
  REAL(pf),                    INTENT(IN) :: G, potsoft, potsoft2

  REAL(pf), DIMENSION(dim) :: Fab, Rab
  INTEGER :: a, b
  REAL(pf) :: distancia, distancia_inv
  REAL(pf), DIMENSION(N, dim) :: forcas_seq

  forcas_seq(:,:) = 0.0_pf

  DO a = 2, N
    DO b = 1, a - 1
      ! distancia entre os corpos
      Rab = R(b,:) - R(a,:)
      distancia = norm2(Rab)
      IF (potsoft .NE. 0) THEN
        distancia = SQRT(distancia*distancia + potsoft2)
      ENDIF
      distancia_inv = 1.0_pf/distancia
      distancia_inv = distancia_inv**3

      ! forca entre os corpos a e b
      Fab = G * m(a) * m(b) * (Rab) * distancia_inv
      ! Adiciona na matriz
      forcas_seq(a,:) = forcas_seq(a,:) + Fab
      forcas_seq(b,:) = forcas_seq(b,:) - Fab
    END DO
  END DO

END FUNCTION forcas_seq

END MODULE funcoes_forca