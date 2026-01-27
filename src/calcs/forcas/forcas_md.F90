! ************************************************************
!! Matriz de forcas (massas diferentes)
!
! Objetivos:
!   Calcula a matriz de forcas a partir das posicoes. Todos os metodos
!   calculam as forcas do mesmo jeito.
!
! Modificado:
!   27 de janeiro de 2026
!
! Autoria:
!   oap
! 
MODULE funcoes_forca_md

  USE tipos
  USE OMP_LIB
  IMPLICIT NONE
  PUBLIC

CONTAINS

! Paralelo (CPU)
FUNCTION forcas_par (m, R, G, N, dim, potsoft2) RESULT(forcas)
  IMPLICIT NONE
  INTEGER,                     INTENT(IN) :: N, dim
  REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
  REAL(pf), DIMENSION(N),      INTENT(IN) :: m
  REAL(pf),                    INTENT(IN) :: G, potsoft2
  REAL(pf), DIMENSION(N, dim) :: forcas
  REAL(pf), DIMENSION(dim)    :: Fab
  INTEGER  :: a, b
  
  forcas = 0.0_pf

  !$OMP PARALLEL SHARED(forcas) PRIVATE(Fab, a, b)
  !$OMP DO
  DO a = 1, N
    DO b = 1, N
      IF (a==b) THEN
        CYCLE
      ENDIF
      
      ! Calculo da forca
      Fab = calcular_forca(G, m, R, a, b, potsoft2)
      
      ! Adiciona na matriz      
      forcas(a,1) = forcas(a,1) + Fab(1)
      forcas(a,2) = forcas(a,2) + Fab(2)
      forcas(a,3) = forcas(a,3) + Fab(3)
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL

END FUNCTION forcas_par

! Sequencial
FUNCTION forcas_seq (m, R, G, N, dim, potsoft2) RESULT(forcas)
  IMPLICIT NONE
  INTEGER,                     INTENT(IN) :: N, dim
  REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
  REAL(pf), DIMENSION(N),      INTENT(IN) :: m
  REAL(pf),                    INTENT(IN) :: G, potsoft2
  REAL(pf), DIMENSION(N, dim) :: forcas
  REAL(pf), DIMENSION(dim)    :: Fab
  INTEGER  :: a, b

  forcas(:,:) = 0.0_pf

  DO a = 2, N
    DO b = 1, a - 1
      ! Calculo da forca
      Fab = calcular_forca(G, m, R, a, b, potsoft2)

      ! Adiciona na matriz
      forcas(a,1) = forcas(a,1) + Fab(1)
      forcas(a,2) = forcas(a,2) + Fab(2)
      forcas(a,3) = forcas(a,3) + Fab(3)

      forcas(b,1) = forcas(b,1) - Fab(1)
      forcas(b,2) = forcas(b,2) - Fab(2)
      forcas(b,3) = forcas(b,3) - Fab(3)
    END DO
  END DO

END FUNCTION forcas_seq

! Calculo das forcas para o sequencial e o paralelo (CPU)
FUNCTION calcular_forca (G, m, R, a, b, potsoft2) RESULT(Fab)
  REAL(pf), INTENT(IN)    :: G, m(:), R(:,:), potsoft2
  INTEGER,  INTENT(IN)    :: a, b
  REAL(pf) :: Fab(3)
  REAL(pf) :: dx, dy, dz, r2, rinv, dist

  dx = R(b,1) - R(a,1)
  dy = R(b,2) - R(a,2)
  dz = R(b,3) - R(a,3)

  dist = dx*dx + dy*dy + dz*dz
  r2 = dist + potsoft2
  rinv = G * m(a) * m(b) / (r2 * SQRT(r2))

  Fab(1) = dx * rinv
  Fab(2) = dy * rinv
  Fab(3) = dz * rinv
END FUNCTION calcular_forca

! GPU
#ifdef USAR_GPU
FUNCTION forcas_par_gpu (m, R, G, N, dim, potsoft2) RESULT(forcas)
  INTEGER,                     INTENT(IN) :: N, dim
  REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
  REAL(pf), DIMENSION(N),      INTENT(IN) :: m
  REAL(pf),                    INTENT(IN) :: G, potsoft2
  INTEGER :: i, j, indice, chunk_size
  REAL(pf) :: dx, dy, dz, r2, rinv
  
  ! Tamanho do bloco para melhor localidade (ajuste conforme sua GPU)
  chunk_size = 64
  
  ! 1. Manter dados na GPU entre chamadas
  !$OMP target data map(to: R, m) map(from: forcas)
  
  ! 2. Inicializacao otimizada
  !$OMP target teams distribute parallel do simd num_teams(N/256) thread_limit(256)
  DO i = 1, N
    forcas(i,:) = 0.0_pf
  END DO
  
  ! 3. Calculo das forcas
  !$OMP target teams distribute parallel do private(dx, dy, dz, r2, rinv) &
  !$OMP num_teams(480) thread_limit(256)
  DO i = 1, N
    ! Loop interno sem collapse
    DO j = i+1, N
      dx = R(j,1) - R(i,1)
      dy = R(j,2) - R(i,2)
      dz = R(j,3) - R(i,3)
      
      r2 = dx*dx + dy*dy + dz*dz + potsoft2
      rinv = G * m(i) * m(j) / (r2 * SQRT(r2))

      !$OMP atomic update
      forcas(i,1) = forcas(i,1) + dx * rinv
      !$OMP atomic update
      forcas(i,2) = forcas(i,2) + dy * rinv
      !$OMP atomic update
      forcas(i,3) = forcas(i,3) + dz * rinv
      
      !$OMP atomic update
      forcas(j,1) = forcas(j,1) - dx * rinv
      !$OMP atomic update
      forcas(j,2) = forcas(j,2) - dy * rinv
      !$OMP atomic update
      forcas(j,3) = forcas(j,3) - dz * rinv
    END DO
  END DO
  
  !$OMP end target data
END FUNCTION forcas_par_gpu
#endif

END MODULE funcoes_forca_md