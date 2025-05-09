! ************************************************************
!! SIMULACAO: SORTEIO
!
! Objetivos:
!   Simulacoes a partir do sorteio de valores iniciais.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
MODULE simulacao_sorteio

  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  USE OMP_LIB
  USE simulacao
  USE condicoesIniciais
  USE arquivos
  USE json_module
  USE arquivos_json
  IMPLICIT NONE
  PRIVATE
  PUBLIC simular_sorteio, sorteio_salvar

CONTAINS

! ************************************************************
!! Metodo principal
!
! Objetivos:
!   Aplica o sorteio e faz a simulacao.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE simular_sorteio (arquivo)
  CHARACTER(LEN=*), INTENT(IN) :: arquivo
  TYPE(json_value), POINTER :: infos
  CHARACTER(LEN=:), ALLOCATABLE :: modo
  ! Vetores
  REAL(pf), allocatable :: massas(:), posicoes(:,:), momentos(:,:)

  ! Le o arquivo de configuracoes
  CALL ler_json(arquivo, infos)
  modo = json_get_string(infos, "modo")

  ! Faz o condicionamento
  CALL condicionar(infos, massas, posicoes, momentos, modo)

  ! Roda a simulacao no intervalo [t0, tf]
  CALL rodar_simulacao (infos, massas, posicoes, momentos)

END SUBROUTINE simular_sorteio

! ************************************************************
!! Sorteia e salva
!
! Objetivos:
!   Sorteia os valores iniciais e salva com um preset de 
!   valores iniciais.
!
! Modificado:
!   01 de maio de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE sorteio_salvar (arquivo_in)
  IMPLICIT NONE
  CHARACTER(LEN=*) :: arquivo_in
  REAL(pf), ALLOCATABLE :: massas(:), posicoes(:,:), momentos(:,:)
  TYPE(json_value), POINTER :: infos
  CHARACTER(LEN=:), ALLOCATABLE :: modo

  ! Le o arquivo de configuracoes
  CALL ler_json(arquivo_in, infos)
  modo = json_get_string(infos, "modo")

  ! Faz o condicionamento
  CALL condicionar(infos, massas, posicoes, momentos, modo)

  ! Verifica se o diretorio de saida existe
  CALL diretorio_vi()

  ! Agora salva
  CALL salvar_auto_vi(infos, massas, posicoes, momentos)
END SUBROUTINE sorteio_salvar

END module