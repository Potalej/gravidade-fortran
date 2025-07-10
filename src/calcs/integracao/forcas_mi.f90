! ************************************************************
!! Matriz de forcas (massas iguais)
!
! Objetivos:
!   Calcula a matriz de forcas a partir das posicoes. Todos os metodos
!   calculam as forcas do mesmo jeito.
!
! Modificado:
!   03 de junho de 2025
!
! Autoria:
!   oap
! 
MODULE funcoes_forca_mi

  USE tipos
  USE OMP_LIB
  IMPLICIT NONE
  PUBLIC

CONTAINS

! Paralelo
FUNCTION forcas_mi_par (R, G, N, dim, potsoft, potsoft2, distancias) RESULT(forcas)
  IMPLICIT NONE
  INTEGER,                     INTENT(IN) :: N, dim
  REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
  REAL(pf),                    INTENT(IN) :: G, potsoft, potsoft2
  
  REAL(pf), DIMENSION(INT(N*(N-1)/2)) :: distancias
  INTEGER :: indice

  REAL(pf), DIMENSION(dim) :: Fab, Rab
  INTEGER :: a, b, tid
  REAL(pf) :: distancia, distancia3, distancia_inv
  REAL(pf), DIMENSION(N, dim) :: forcas
  
  forcas = 0.0_pf
  indice = 1

  !$OMP PARALLEL SHARED(forcas) PRIVATE(Fab, Rab, distancia, distancia_inv, a, b, tid)
  !$OMP DO
  DO a = 1, N
    DO b = 1, N
      IF (a==b) THEN
        CYCLE
      ENDIF
      ! distancia entre os corpos
      Rab = R(b,:) - R(a,:)
      distancia = norm2(Rab)

      distancias(indice) = distancia
      indice = indice + 1

      IF (potsoft .NE. 0) THEN
        distancia = SQRT(distancia*distancia + potsoft2)
      ENDIF
      distancia3 = distancia * distancia * distancia
      distancia_inv = 1.0_pf/distancia3

      ! forca entre os corpos a e b
      Fab = G * Rab * distancia_inv
      
      ! Adiciona na matriz      
      forcas(a,:) = forcas(a,:) + Fab
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL

END FUNCTION forcas_mi_par

! Sequencial 
FUNCTION forcas_mi_seq (R, G, N, dim, potsoft, potsoft2, distancias) RESULT(forcas)
  IMPLICIT NONE
  INTEGER,                     INTENT(IN) :: N, dim
  REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
  REAL(pf),                    INTENT(IN) :: G, potsoft, potsoft2

  REAL(pf), DIMENSION(INT(N*(N-1)/2)) :: distancias
  INTEGER :: indice

  REAL(pf), DIMENSION(dim) :: Fab, Rab
  INTEGER :: a, b
  REAL(pf) :: distancia, distancia3, distancia_inv
  REAL(pf), DIMENSION(N, dim) :: forcas

  forcas(:,:) = 0.0_pf
  indice = 1

  DO a = 2, N
    DO b = 1, a - 1
      ! distancia entre os corpos
      Rab = R(b,:) - R(a,:)
      distancia = Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3)
      
      distancias(indice) = SQRT(distancia)
      indice = indice + 1
      
      IF (potsoft .NE. 0) THEN
        distancia = distancia + potsoft2
      ENDIF
      distancia3 = distancia * SQRT(distancia)
      distancia_inv = 1.0_pf/distancia3

      ! forca entre os corpos a e b
      Fab = G * Rab * distancia_inv
      ! Adiciona na matriz
      forcas(a,:) = forcas(a,:) + Fab
      forcas(b,:) = forcas(b,:) - Fab
    END DO
  END DO

END FUNCTION forcas_mi_seq

END MODULE funcoes_forca_mi