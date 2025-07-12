! ************************************************************
!! Matriz de forcas (massas iguais)
!
! Objetivos:
!   Calcula a matriz de forcas a partir das posicoes. Todos os metodos
!   calculam as forcas do mesmo jeito.
!
! Modificado:
!   10 de julho de 2025
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
FUNCTION forcas_mi_par (R, G, N, dim, potsoft2, distancias) RESULT(forcas)
  IMPLICIT NONE
  INTEGER,                     INTENT(IN) :: N, dim
  REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
  REAL(pf),                    INTENT(IN) :: G, potsoft2
  
  REAL(pf), DIMENSION(INT(N*(N-1)/2)) :: distancias
  INTEGER :: indice

  REAL(pf), DIMENSION(dim) :: Fab
  INTEGER :: a, b
  REAL(pf) :: distancia
  REAL(pf), DIMENSION(N, dim) :: forcas
  
  forcas = 0.0_pf

  !$OMP PARALLEL SHARED(forcas) PRIVATE(Fab, distancia, a, b, indice)
  !$OMP DO
  DO a = 1, N
    DO b = 1, N
      IF (a==b) THEN
        CYCLE
      ENDIF

      ! Calculo da forca
      Fab = calcular_forca(G, R, a, b, potsoft2, distancia)

      IF (a > b) THEN
        indice = (a-1)*(a-2)/2 + (b-1)
      ELSE
        indice = (b-1)*(b-2)/2 + (a-1)
      ENDIF
      distancias(indice+1) = SQRT(distancia)
      
      ! Adiciona na matriz      
      forcas(a,:) = forcas(a,:) + Fab
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL

END FUNCTION forcas_mi_par

! Sequencial 
FUNCTION forcas_mi_seq (R, G, N, dim, potsoft2, distancias) RESULT(forcas)
  IMPLICIT NONE
  INTEGER,                     INTENT(IN) :: N, dim
  REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
  REAL(pf),                    INTENT(IN) :: G, potsoft2

  REAL(pf), DIMENSION(INT(N*(N-1)/2)) :: distancias
  INTEGER :: indice

  REAL(pf), DIMENSION(dim) :: Fab
  INTEGER :: a, b
  REAL(pf) :: distancia
  REAL(pf), DIMENSION(N, dim) :: forcas

  forcas(:,:) = 0.0_pf
  indice = 1

  DO a = 2, N
    DO b = 1, a - 1
      ! Calculo da forca
      Fab = calcular_forca(G, R, a, b, potsoft2, distancia)

      ! Salva a distancia
      distancias(indice) = SQRT(distancia)
      indice = indice + 1
      
      ! Adiciona na matriz
      forcas(a,:) = forcas(a,:) + Fab
      forcas(b,:) = forcas(b,:) - Fab
    END DO
  END DO

END FUNCTION forcas_mi_seq

! Calculo das forcas em ambos os casos
FUNCTION calcular_forca (G, R, a, b, potsoft2, dist) RESULT(Fab)
  REAL(pf), INTENT(IN) :: G, R(:,:), potsoft2
  INTEGER,  INTENT(IN) :: a, b
  REAL(pf), INTENT(OUT) :: dist
  REAL(pf) :: dist_pot
  REAL(pf) :: distancia3, dist_inv, Fab(3)

  Fab = R(b,:) - R(a,:)
  
  dist = Fab(1)*Fab(1) + Fab(2)*Fab(2) + Fab(3)*Fab(3)
  dist_pot = dist + potsoft2
  distancia3 = dist_pot * SQRT(dist_pot)
  dist_inv = G * 1.0_pf/distancia3

  Fab = Fab * dist_inv
END FUNCTION calcular_forca

END MODULE funcoes_forca_mi