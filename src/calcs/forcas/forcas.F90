! ************************************************************
!! MODULO CENTRALIZADOR DAS FORCAS
!
! Objetivos:
!   Ha tres formas de calcular as forcas: sequencialmente,
!   paralelo (CPU) e paralelo (GPU). Tambem ha duas formas de
!   lidar com as massas: caso com massas diferentes (MD) e com
!   massas iguais (MI), que estao separadas em dois modulos.
!   Alem disso, tambem eh possivel calcular as forcas
!   diretamente ou utilizando uma arvore (octree) atraves do
!   criterio de Barnes-Hut.
! 
!   O que este modulo faz eh auxiliar a classe INTEGRACAO na
!   inicializacao das funcoes de forca, trazendo para ca a
!   tarefa de apontar corretamente os ponteiros das funcoes.
!
! Modificado:
!   26 de maio de 2026
!
! Autoria:
!   oap
! 
MODULE funcoes_forca
  USE tipos
  USE octree_mod

!> Modulos de forca
  USE funcoes_forca_md
  USE funcoes_forca_mi
  USE funcoes_forca_tree

  IMPLICIT NONE
  PUBLIC

!> Interface com o basico das funcoes de forca
  ABSTRACT INTERFACE
    FUNCTION forcas_funcbase (m, R, G, N, dim, potsoft2)
        IMPORT :: pf
        IMPLICIT NONE
        INTEGER,                     INTENT(IN) :: N, dim
        REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
        REAL(pf), DIMENSION(N),      INTENT(IN) :: m
        REAL(pf),                    INTENT(IN) :: G, potsoft2
        REAL(pf), DIMENSION(N, dim) :: forcas_funcbase
    END FUNCTION forcas_funcbase

    FUNCTION forcas_mi_funcbase (R, G, N, dim, potsoft2)
        IMPORT :: pf
        IMPLICIT NONE
        INTEGER,                     INTENT(IN) :: N, dim
        REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
        REAL(pf),                    INTENT(IN) :: G, potsoft2
        REAL(pf), DIMENSION(N, dim) :: forcas_mi_funcbase
    END FUNCTION forcas_mi_funcbase

    FUNCTION forcas_tree_funcbase (m, R, G, N, dim, potsoft, theta2, octree)
        IMPORT :: pf, OctreeType
        IMPLICIT NONE
        INTEGER,                     INTENT(IN) :: N, dim
        REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
        REAL(pf), DIMENSION(N),      INTENT(IN) :: m
        REAL(pf),                    INTENT(IN) :: G, potsoft, theta2
        CLASS(OctreeType),        INTENT(INOUT) :: octree
        REAL(pf), DIMENSION(N, dim) :: forcas_tree_funcbase
    END FUNCTION forcas_tree_funcbase
  END INTERFACE

CONTAINS

! ************************************************************
!! Inicializador
!
! Modificado:
!   26 de maio de 2026
!
! Autoria:
!   oap
! 
SUBROUTINE inicializar_forcas (mi, pcpu, pgpu, tree, f, f_mi, f_t)
  LOGICAL, INTENT(IN) :: mi, pcpu, pgpu, tree
  PROCEDURE(forcas_funcbase), POINTER, INTENT(OUT)    :: f
  PROCEDURE(forcas_mi_funcbase), POINTER, INTENT(OUT) :: f_mi
  PROCEDURE(forcas_tree_funcbase), POINTER, INTENT(OUT) :: f_t

  IF (pgpu) THEN
#ifdef USAR_GPU
    IF (mi) THEN
      f_mi => forcas_mi_par_gpu
    ELSE
      f => forcas_par_gpu
    ENDIF
#else
    WRITE (*,*) "GPU nao compilada. Recompile o programa com -DUSAR_GPU=ON"
    WRITE (*,*) "ou desative a opcao de gpu."
    STOP 0
#endif
  ELSE
    IF (pcpu) THEN
      IF (tree) THEN
        f_t => forcas_par_tree
      ELSE
        IF (mi) THEN
          f_mi => forcas_mi_par
        ELSE
          f => forcas_par
        ENDIF
      ENDIF
    ELSE
      IF (tree) THEN
        f_t => forcas_seq_tree
      ELSE
        IF (mi) THEN
          f_mi => forcas_mi_seq
        ELSE
          f => forcas_seq
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  
END SUBROUTINE inicializar_forcas

END MODULE funcoes_forca