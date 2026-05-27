! ************************************************************
!! Matriz de forcas (octree)
!
! Objetivos:
!   Calculo das forcas entre os corpos utilizando uma octree e
!   o criterio de Barnes e Hut a partir do parametro `theta`.
!
! Modificado:
!   26 de maio de 2026
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
FUNCTION forcas_seq_tree (m, R, G, N, dim, potsoft, theta2, octree) RESULT(forcas)
  IMPLICIT NONE
  INTEGER,                     INTENT(IN) :: N, dim
  REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
  REAL(pf), DIMENSION(N),      INTENT(IN) :: m
  REAL(pf),                    INTENT(IN) :: G, potsoft, theta2
  REAL(pf), DIMENSION(dim)    :: Fab
  REAL(pf), DIMENSION(N, dim) :: forcas
  INTEGER  :: a, b
  CLASS(OctreeType), INTENT(INOUT) :: octree

  forcas(:,:) = 0.0_pf

  CALL octree % init(m, R)

  DO a = 1, N
    CALL octree % forces(a, theta2, G, potsoft, forcas(a,:))
  END DO

END FUNCTION forcas_seq_tree

! Paralelo com massas diferentes
FUNCTION forcas_par_tree (m, R, G, N, dim, potsoft, theta2, octree) RESULT(forcas)
  IMPLICIT NONE
  INTEGER,                     INTENT(IN) :: N, dim
  REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
  REAL(pf), DIMENSION(N),      INTENT(IN) :: m
  REAL(pf),                    INTENT(IN) :: G, potsoft, theta2
  REAL(pf), DIMENSION(N, dim) :: forcas
  INTEGER  :: a, b
  CLASS(OctreeType), INTENT(INOUT) :: octree

  forcas(:,:) = 0.0_pf

  CALL octree % init(m, R)

  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP SHARED(forcas, octree, theta2, G, potsoft, N) &
  !$OMP PRIVATE(a) SCHEDULE(DYNAMIC)
  DO a = 1, N
    CALL octree % forces(a, theta2, G, potsoft, forcas(a,:))
  END DO
  !$OMP END PARALLEL DO

END FUNCTION forcas_par_tree

END MODULE funcoes_forca_tree