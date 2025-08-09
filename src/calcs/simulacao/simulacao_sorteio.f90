! ************************************************************
!! SIMULACAO: SORTEIO
!
! Objetivos:
!   Simulacoes a partir do sorteio de valores iniciais.
!
! Modificado:
!   20 de julho de 2025
!
! Autoria:
!   oap
! 
MODULE simulacao_sorteio

  USE tipos
  USE simulador_mod
  USE condicoes_iniciais
  USE arquivos_mod
  USE json_module
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
!   20 de julho de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE simular_sorteio (arquivo)
  CHARACTER(LEN=*), INTENT(IN) :: arquivo
  TYPE(json_value), POINTER :: infos
  CHARACTER(LEN=:), ALLOCATABLE :: modo
  REAL(pf) :: eps
  ! Vetores
  REAL(pf), allocatable :: massas(:), posicoes(:,:), momentos(:,:)

  ! Le o arquivo de configuracoes
  CALL ler_json(arquivo, infos)
  modo = json_get_string(infos, "modo")
  eps = json_get_float(infos, 'integracao.amortecedor')

  ! Faz o condicionamento
  CALL condicionar(infos, massas, posicoes, momentos, modo, eps)

  ! Roda a simulacao no intervalo [t0, tf]
  CALL rodar_simulacao(infos, massas, posicoes, momentos)

END SUBROUTINE simular_sorteio

! ************************************************************
!! Sorteia e salva
!
! Objetivos:
!   Sorteia os valores iniciais e salva com um preset de 
!   valores iniciais.
!
! Modificado:
!   20 de julho de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE sorteio_salvar (arquivo_in)
  IMPLICIT NONE
  CHARACTER(LEN=*) :: arquivo_in
  REAL(pf) :: eps
  REAL(pf), ALLOCATABLE :: massas(:), posicoes(:,:), momentos(:,:)
  TYPE(json_value), POINTER :: infos
  CHARACTER(LEN=:), ALLOCATABLE :: modo

  ! Le o arquivo de configuracoes
  CALL ler_json(arquivo_in, infos)
  modo = json_get_string(infos, "modo")
  eps = json_get_float(infos, 'integracao.amortecedor')

  ! Faz o condicionamento
  CALL condicionar(infos, massas, posicoes, momentos, modo, eps)

  ! Verifica se o diretorio de saida existe
  CALL diretorio_vi()

  ! Agora salva
  CALL salvar_auto_vi(infos, massas, posicoes, momentos)
END SUBROUTINE sorteio_salvar

END module