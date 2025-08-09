! ************************************************************
!! MODULO CENTRALIZADOR DAS FORCAS
!
! Objetivos:
!   Ha tres formas de calcular as forcas: sequencialmente,
!   paralelo (CPU) e paralelo (GPU). Tambem ha duas formas de
!   lidar com as massas: caso com massas diferentes (MD) e com
!   massas iguais (MI), que estao separadas em dois modulos.
! 
!   O que este modulo faz eh auxiliar a classe INTEGRACAO na
!   inicializacao das funcoes de forca, trazendo para ca a
!   tarefa de apontar corretamente os ponteiros das funcoes.
!
! Modificado:
!   08 de agosto de 2025
!
! Autoria:
!   oap
! 
MODULE funcoes_forca
  USE tipos

!> Modulos de forca
  USE funcoes_forca_md
  USE funcoes_forca_mi

  IMPLICIT NONE
  PUBLIC

!> Interface com o basico das funcoes de forca
  ABSTRACT INTERFACE
    FUNCTION forcas_funcbase (m, R, G, N, dim, potsoft2, distancias)
        IMPORT :: pf
        IMPLICIT NONE
        INTEGER,                     INTENT(IN) :: N, dim
        REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
        REAL(pf), DIMENSION(N),      INTENT(IN) :: m
        REAL(pf),                    INTENT(IN) :: G, potsoft2
        REAL(pf), DIMENSION(INT(N*(N-1)/2)), INTENT(INOUT) :: distancias
        REAL(pf), DIMENSION(N, dim) :: forcas_funcbase
    END FUNCTION forcas_funcbase

    FUNCTION forcas_mi_funcbase (R, G, N, dim, potsoft2, distancias)
        IMPORT :: pf
        IMPLICIT NONE
        INTEGER,                     INTENT(IN) :: N, dim
        REAL(pf), DIMENSION(N, dim), INTENT(IN) :: R
        REAL(pf),                    INTENT(IN) :: G, potsoft2
        REAL(pf), DIMENSION(INT(N*(N-1)/2)), INTENT(INOUT) :: distancias
        REAL(pf), DIMENSION(N, dim) :: forcas_mi_funcbase
    END FUNCTION forcas_mi_funcbase
  END INTERFACE

CONTAINS

! ************************************************************
!! Inicializador
!
! Modificado:
!   08 de agosto de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE inicializar_forcas (mi, pcpu, pgpu, f, f_mi)
  LOGICAL, INTENT(IN) :: mi, pcpu, pgpu
  PROCEDURE(forcas_funcbase), POINTER, INTENT(OUT)    :: f
  PROCEDURE(forcas_mi_funcbase), POINTER, INTENT(OUT) :: f_mi

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
      IF (mi) THEN
        f_mi => forcas_mi_par
      ELSE
        f => forcas_par
      ENDIF
    ELSE
      IF (mi) THEN
        f_mi => forcas_mi_seq
      ELSE
        f => forcas_seq
      ENDIF
    ENDIF
  ENDIF
  
END SUBROUTINE inicializar_forcas

END MODULE funcoes_forca