! ************************************************************
!! SIMULACAO: VALORES INICIAIS (VI)
!
! Objetivos:
!   Simulacoes a partir diretamente de valores iniciais.
!
! Modificado:
!   13 de setembro de 2026
!
! Autoria:
!   oap
! 
MODULE simulacao_vi

  USE tipos
  USE simulador_mod
  USE arquivos_mod
  IMPLICIT NONE
  PRIVATE
  PUBLIC simular_vi

CONTAINS

! ************************************************************
!! Validacao dos parametros do arquivo de entrada de sorteio
!
! Modificado:
!   13 de setembro de 2026
!
! Autoria:
!   oap
!
SUBROUTINE parametros (infos)
  TYPE(json_value), POINTER, INTENT(IN) :: infos
  CHARACTER(LEN=40) :: pars(26)

  ! gerais
  pars(1) = "modo"
  pars(2) = "nome"
  pars(3) = "N"
  pars(4) = "G"
  pars(5) = "paralelo"
  pars(6) = "gpu"
  pars(7) = "massas_iguais"
  pars(8) = "exibir"

  ! integracao
  pars(9)  = "integracao.metodo"
  pars(10) = "integracao.timestep"
  pars(11) = "integracao.amortecedor"
  pars(12) = "integracao.t0"
  pars(13) = "integracao.tf"
  pars(14) = "integracao.checkpoints"
  pars(15) = "integracao.tree"
  pars(16) = "integracao.theta"

  ! colisoes
  pars(17) = "colisoes.colidir"
  pars(18) = "colisoes.metodo"
  pars(19) = "colisoes.densidade"
  pars(20) = "colisoes.permitir_choque_inicial"

  ! correcao
  pars(21) = "correcao.corrigir"
  pars(22) = "correcao.margem_erro"
  pars(23) = "correcao.max_num_tentativas"

  ! valores iniciais
  pars(24) = "valores_iniciais.massas"
  pars(25) = "valores_iniciais.posicoes"
  pars(26) = "valores_iniciais.momentos"

  CALL validar_json(infos, pars)
END SUBROUTINE

! ************************************************************
!! Metodo principal
!
! Objetivos:
!   Faz a simulacao.
!
! Modificado:
!   13 de setembro de 2026
!
! Autoria:
!   oap
! 
SUBROUTINE simular_vi (arquivo, out_dir, out_ext)
  CHARACTER(LEN=*), INTENT(IN) :: arquivo
  CHARACTER(LEN=*), INTENT(IN) :: out_dir, out_ext

  TYPE(json_value), POINTER :: infos
  REAL(pf), ALLOCATABLE :: pos3(:), mom3(:)
  REAL(pf), ALLOCATABLE :: massas(:), posicoes(:,:), momentos(:,:)
  INTEGER :: a
  CHARACTER(32) :: a_string

  ! Le o arquivo de valores iniciais
  CALL ler_json(arquivo, infos)

  ! valida o arquivo
  CALL parametros(infos)

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
  CALL rodar_simulacao(out_dir, out_ext, infos, massas, posicoes, momentos)

END SUBROUTINE simular_vi

END MODULE simulacao_vi