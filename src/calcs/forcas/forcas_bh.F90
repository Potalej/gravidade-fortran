! ************************************************************
!! Matriz de forcas (octree usando Barnes-Hut)
!
! Objetivos:
!   Calculo das forcas entre os corpos utilizando uma octree e
!   o criterio de Barnes e Hut a partir do parametro `theta`.
!
! Modificado:
!   07 de outubro de 2026
!
! Autoria:
!   oap
! 
MODULE funcoes_forca_bh

  USE tipos
  USE octree_mod
  USE OMP_LIB
  IMPLICIT NONE
  PUBLIC

CONTAINS

! Sequencial com massas diferentes
FUNCTION forcas_seq_bh (m, R, G, N, dim, potsoft2, theta2, octree) RESULT(forcas)
  IMPLICIT NONE
  INTEGER,                     INTENT(IN) :: N, dim
  REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
  REAL(pf), DIMENSION(N),      INTENT(IN) :: m
  REAL(pf),                    INTENT(IN) :: G, potsoft2, theta2
  REAL(pf), DIMENSION(dim)    :: Fab
  REAL(pf), DIMENSION(N, dim) :: forcas
  INTEGER  :: a, b
  CLASS(OctreeType), INTENT(INOUT) :: octree

  CALL octree % init(R(:,1),R(:,2),R(:,3))

  DO a = 1, N
    forcas(a,:) = octree % forces(a, theta2, G, potsoft2)
  END DO

END FUNCTION

! Paralelo com massas diferentes
FUNCTION forcas_par_bh (m, R, G, N, dim, potsoft2, theta2, octree) RESULT(forcas)
  IMPLICIT NONE
  INTEGER,                     INTENT(IN) :: N, dim
  REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
  REAL(pf), DIMENSION(N),      INTENT(IN) :: m
  REAL(pf),                    INTENT(IN) :: G, potsoft2, theta2
  REAL(pf), DIMENSION(N, dim) :: forcas
  REAL(pf), DIMENSION(dim, N) :: forcas_dehnen
  INTEGER  :: a, b
  CLASS(OctreeType), INTENT(INOUT) :: octree

  CALL octree % init(R(:,1),R(:,2),R(:,3))

  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP SHARED(forcas, octree, theta2, G, potsoft2, N) &
  !$OMP PRIVATE(a) SCHEDULE(DYNAMIC)
  DO a = 1, N
    forcas(a,:) = octree % forces(a, theta2, G, potsoft2)
  END DO
  !$OMP END PARALLEL DO

END FUNCTION

END MODULE funcoes_forca_bh