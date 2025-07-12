! ************************************************************
!! SIMULACAO: VALORES INICIAIS (VI)
!
! Objetivos:
!   Simulacoes a partir diretamente de valores iniciais.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
MODULE simulacao_vi

  USE tipos
  USE OMP_LIB
  USE simulacao
  USE arquivos_json
  USE json_utils_mod
  IMPLICIT NONE
  PRIVATE
  PUBLIC simular_vi

CONTAINS

! ************************************************************
!! Metodo principal
!
! Objetivos:
!   Faz a simulacao.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE simular_vi (arquivo)
  CHARACTER(LEN=*), INTENT(INOUT) :: arquivo

  TYPE(json_value), POINTER :: infos
  REAL(pf), ALLOCATABLE :: pos3(:), mom3(:)
  REAL(pf), ALLOCATABLE :: massas(:), posicoes(:,:), momentos(:,:)
  INTEGER :: a
  CHARACTER(32) :: a_string

  ! Le o arquivo de valores iniciais
  CALL ler_json(arquivo, infos)

  ! Le os valores iniciais
  massas = json_get_float_vec(infos, 'valores_iniciais.massas')

  ALLOCATE(posicoes(SIZE(massas),3))
  ALLOCATE(momentos(SIZE(massas),3))

  DO a=1, SIZE(massas)
    WRITE(a_string, *) a
    pos3 = json_get_float_vec(infos, 'valores_iniciais.posicoes['//a_string//']')
    mom3 = json_get_float_vec(infos, 'valores_iniciais.momentos['//a_string//']')
    posicoes(a,:) = pos3
    momentos(a,:) = mom3
  END DO

  ! Roda a simulacao
  CALL rodar_simulacao(infos, massas, posicoes, momentos)

END SUBROUTINE simular_vi

END MODULE simulacao_vi