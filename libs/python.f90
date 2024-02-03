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
subroutine rodar_plot (diretorio)
  IMPLICIT NONE
  CHARACTER(LEN=*) :: diretorio
  CALL system ('py python/plot.py '//TRIM(diretorio))
end subroutine rodar_plot