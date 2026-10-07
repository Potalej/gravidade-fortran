! ************************************************************
!! Matriz de forcas (octree usando Dehnen)
!
! Objetivos:
!   Calculo das forcas entre os corpos utilizando uma octree e
!   o algoritmo de Dehnen.
!
! Modificado:
!   07 de outubro de 2026
!
! Autoria:
!   oap
! 
MODULE funcoes_forca_dehnen

  USE tipos
  USE octree_mod
  USE OMP_LIB
  IMPLICIT NONE
  PUBLIC

CONTAINS

! Sequencial com massas diferentes
FUNCTION forcas_seq_dehnen (m, R, G, N, dim, potsoft2, theta2, octree) RESULT(forcas)
  IMPLICIT NONE
  INTEGER,                     INTENT(IN) :: N, dim
  REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
  REAL(pf), DIMENSION(N),      INTENT(IN) :: m
  REAL(pf),                    INTENT(IN) :: G, potsoft2, theta2
  REAL(pf), DIMENSION(dim)    :: Fab
  REAL(pf), DIMENSION(N, dim) :: forcas
  REAL(pf), DIMENSION(dim, N) :: forcas_dehnen
  INTEGER  :: a, b
  CLASS(OctreeType), INTENT(INOUT) :: octree

  CALL octree % init(R(:,1),R(:,2),R(:,3))
  call octree % dehnen_eval(theta2, potsoft2, G, forcas_dehnen)
  
  forcas(:,1) = forcas_dehnen(1,:)
  forcas(:,2) = forcas_dehnen(2,:)
  forcas(:,3) = forcas_dehnen(3,:)

END FUNCTION

! Paralelo
! Por enquanto faz a mesma coisa, nao sei como paralelizar ainda
FUNCTION forcas_par_dehnen (m, R, G, N, dim, potsoft2, theta2, octree) RESULT(forcas)
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
  call octree % dehnen_eval(theta2, potsoft2, G, forcas_dehnen)
  
  forcas(:,1) = forcas_dehnen(1,:)
  forcas(:,2) = forcas_dehnen(2,:)
  forcas(:,3) = forcas_dehnen(3,:)
END FUNCTION

END MODULE funcoes_forca_dehnen