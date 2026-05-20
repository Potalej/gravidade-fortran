! ************************************************************
!! Matriz de forcas (octree)
!
! Objetivos:
!   Calculo das forcas entre os corpos utilizando uma octree e
!   o criterio de Barnes e Hut a partir do parametro `theta`.
!
! Modificado:
!   20 de maio de 2026
!
! Autoria:
!   oap
! 
MODULE funcoes_forca_tree

  USE tipos
  USE octree_mod
  USE OMP_LIB
  IMPLICIT NONE
  PUBLIC

CONTAINS

! Sequencial com massas diferentes
FUNCTION forcas_seq_tree (m, R, G, N, dim, potsoft, theta2) RESULT(forcas)
  IMPLICIT NONE
  INTEGER,                     INTENT(IN) :: N, dim
  REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
  REAL(pf), DIMENSION(N),      INTENT(IN) :: m
  REAL(pf),                    INTENT(IN) :: G, potsoft, theta2
  REAL(pf), DIMENSION(dim)    :: Fab
  REAL(pf), DIMENSION(N, dim) :: forcas
  REAL(pf), DIMENSION(N) :: x, y, z
  INTEGER  :: a, b
  CLASS(OctreeType), ALLOCATABLE :: tree

  forcas(:,:) = 0.0_pf

  x = R(:, 1)
  y = R(:, 2)
  z = R(:, 3)

  ALLOCATE(tree)
  CALL tree % init(m, x, y, z)

  DO a = 1, N
    forcas(a,:) = tree % forces(a, theta2, G, potsoft)
  END DO

  DEALLOCATE(tree)

END FUNCTION forcas_seq_tree

! Paralelo com massas diferentes
FUNCTION forcas_par_tree (m, R, G, N, dim, potsoft, theta2) RESULT(forcas)
  IMPLICIT NONE
  INTEGER,                     INTENT(IN) :: N, dim
  REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
  REAL(pf), DIMENSION(N),      INTENT(IN) :: m
  REAL(pf),                    INTENT(IN) :: G, potsoft, theta2
  REAL(pf), DIMENSION(N, dim) :: forcas
  REAL(pf), DIMENSION(N) :: x, y, z
  INTEGER  :: a, b
  CLASS(OctreeType), ALLOCATABLE :: tree

  forcas(:,:) = 0.0_pf

  x = R(:, 1)
  y = R(:, 2)
  z = R(:, 3)

  ALLOCATE(tree)
  CALL tree % init(m, x, y, z)

  !$OMP PARALLEL SHARED(forcas) PRIVATE(a)
  !$OMP DO
  DO a = 1, N
    forcas(a,:) = tree % forces(a, theta2, G, potsoft)
  END DO
  !$OMP END DO
  !$OMP END PARALLEL

  DEALLOCATE(tree)

END FUNCTION forcas_par_tree

! Sequencial com massas iguais
FUNCTION forcas_mi_seq_tree (R, G, N, dim, potsoft, theta2) RESULT(forcas)
  IMPLICIT NONE
  INTEGER,                     INTENT(IN) :: N, dim
  REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
  REAL(pf),                    INTENT(IN) :: G, potsoft, theta2
  REAL(pf), DIMENSION(dim)    :: Fab
  REAL(pf), DIMENSION(N, dim) :: forcas
  REAL(pf), DIMENSION(N) :: x, y, z, m
  INTEGER  :: a, b
  CLASS(OctreeType), ALLOCATABLE :: tree

  forcas(:,:) = 0.0_pf

  m = 1.0_pf
  x = R(:, 1)
  y = R(:, 2)
  z = R(:, 3)

  ALLOCATE(tree)
  CALL tree % init(m, x, y, z)

  DO a = 1, N
    forcas(a,:) = tree % forces(a, theta2, G, potsoft)
  END DO
  DEALLOCATE(tree)

END FUNCTION forcas_mi_seq_tree

! Paralelo (CPU) com massas iguais
FUNCTION forcas_mi_par_tree (R, G, N, dim, potsoft, theta2) RESULT(forcas)
  IMPLICIT NONE
  INTEGER,                     INTENT(IN) :: N, dim
  REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
  REAL(pf),                    INTENT(IN) :: G, potsoft, theta2
  REAL(pf), DIMENSION(N, dim) :: forcas
  REAL(pf), DIMENSION(N) :: x, y, z, m
  INTEGER  :: a, b
  CLASS(OctreeType), ALLOCATABLE :: tree

  forcas(:,:) = 0.0_pf

  m = 1.0_pf
  x = R(:, 1)
  y = R(:, 2)
  z = R(:, 3)

  ALLOCATE(tree)
  CALL tree % init(m, x, y, z)

  !$OMP PARALLEL SHARED(forcas) PRIVATE(a)
  !$OMP DO
  DO a = 1, N
    forcas(a,:) = tree % forces(a, theta2, G, potsoft)
  END DO
  !$OMP END DO
  !$OMP END PARALLEL

  DEALLOCATE(tree)

END FUNCTION forcas_mi_par_tree

END MODULE funcoes_forca_tree