! ************************************************************
!! Matriz de forcas (octree usando Barnes-Hut otimizado)
!
! Objetivos:
!   Calculo das forcas entre os corpos utilizando uma octree e
!   o criterio de Barnes e Hut a partir do parametro `theta`.
!   Nesta versao, as listas de interacao sao montadas para
!   grupos de particulas, em vez de individualmente.
!
! Modificado:
!   03 de novembro de 2026
!
! Autoria:
!   oap
! 
MODULE funcoes_forca_bh_groups

  USE tipos
  USE octree_mod
  USE OMP_LIB
  IMPLICIT NONE
  PUBLIC

CONTAINS

! Sequencial com massas diferentes
FUNCTION forcas_seq_bh_groups (m, R, G, N, dim, potsoft2, theta2, octree) RESULT(forcas)
  IMPLICIT NONE
  INTEGER,                     INTENT(IN) :: N, dim
  REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
  REAL(pf), DIMENSION(N),      INTENT(IN) :: m
  REAL(pf),                    INTENT(IN) :: G, potsoft2, theta2
  REAL(pf), DIMENSION(N, dim) :: forcas
  REAL(pf), DIMENSION(dim, N) :: forcas_groups
  CLASS(OctreeType), INTENT(INOUT) :: octree

  CALL octree % init(R(:,1),R(:,2),R(:,3))
  CALL octree % forces_groups(theta2, G, potsoft2, forcas_groups)
  forcas(:,1) = forcas_groups(1,:)
  forcas(:,2) = forcas_groups(2,:)
  forcas(:,3) = forcas_groups(3,:)
END FUNCTION

! Paralelo com massas diferentes
FUNCTION forcas_par_bh_groups (m, R, G, N, dim, potsoft2, theta2, octree) RESULT(forcas)
  IMPLICIT NONE
  INTEGER,                     INTENT(IN) :: N, dim
  REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
  REAL(pf), DIMENSION(N),      INTENT(IN) :: m
  REAL(pf),                    INTENT(IN) :: G, potsoft2, theta2
  REAL(pf), DIMENSION(N, dim) :: forcas
  REAL(pf), DIMENSION(dim, N) :: forcas_groups
  CLASS(OctreeType), INTENT(INOUT) :: octree

  CALL octree % init(R(:,1),R(:,2),R(:,3))
  CALL octree % forces_groups(theta2, G, potsoft2, forcas_groups)
  forcas(:,1) = forcas_groups(1,:)
  forcas(:,2) = forcas_groups(2,:)
  forcas(:,3) = forcas_groups(3,:)
END FUNCTION

END MODULE funcoes_forca_bh_groups