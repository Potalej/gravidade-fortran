! *****************************************************************
!! PYTHON
!
! Objetivos:
!   Este arquivo contem helpers em Python para utilizar nas
!   simulacoes, analises e visualizacoes.
! 
! Modificado:
!   02 de fevereiro de 2024
! 
! Autoria:
!   oap
! 

! ************************************************************
!! Rodar plot com matplotlib
!
! Objetivos:
!   Roda o script "plot.py" e exibe um grafico com matplotlib
!
! Modificado:
!   02 de fevereiro de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE rodar_plot (diretorio)
  IMPLICIT NONE
  CHARACTER(LEN=*) :: diretorio
  CALL SYSTEM ('python src/python/plot.py '//TRIM(diretorio))
END SUBROUTINE rodar_plot