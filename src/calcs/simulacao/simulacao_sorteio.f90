! ************************************************************
!! SIMULACAO: SORTEIO
!
! Objetivos:
!   Simulacoes a partir do sorteio de valores iniciais.
!
! Modificado:
!   05 de janeiro de 2026
!
! Autoria:
!   oap
! 
MODULE simulacao_sorteio

  USE tipos
  USE simulador_mod
  USE condicoes_iniciais
  USE sorteio_mod
  USE arquivos_mod
  IMPLICIT NONE
  PRIVATE
  PUBLIC simular_sorteio, sorteio_salvar

CONTAINS

! ************************************************************
!! Geracao de configuracoes a partir do json_value
!
! Objetivos:
!   Converte json_value em objetos sorteio e sorteio_vetores.
!
! Modificado:
!   16 de dezembro de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE configuracoes_sorteio (dados, configs)
  TYPE(json_value), POINTER, INTENT(IN) :: dados
  TYPE(sorteio), POINTER, INTENT(INOUT) :: configs
  TYPE(json_value), POINTER :: sorteio_json
  LOGICAL :: encontrado

  ! Configuracoes basicas
  CALL json % get(dados, "N", configs % N)
  configs % G = json_get_float(dados, "G")
  configs % amortecedor = json_get_float(dados, 'integracao.amortecedor')
  configs % modo = json_get_string(dados, "modo")

  ! Parte de sorteio
  CALL json % get (dados, "sorteio", sorteio_json)

  ! Integrais primeiras desejadas
  configs % ed = json_get_float(sorteio_json, "integrais.energia_total")
  configs % jd = json_get_float(sorteio_json, "integrais.angular_total")
  configs % pd = json_get_float(sorteio_json, "integrais.linear_total")

  ! Configuracoes das massas 
  configs % massas % distribuicao = json_get_string(sorteio_json, "massas.distribuicao")
  configs % massas % regiao = json_get_string(sorteio_json, "massas.regiao")
  configs % massas % intervalo = json_get_float_vec(sorteio_json, "massas.intervalo")
  CALL json % get(sorteio_json, "massas.normalizadas", configs % massas % normalizado, encontrado)
  IF (.NOT. encontrado) configs % massas % normalizado = .FALSE.

  ! Configuracoes das posicoes
  configs % posicoes % distribuicao = json_get_string(sorteio_json, "posicoes.distribuicao")
  configs % posicoes % regiao = json_get_string(sorteio_json, "posicoes.regiao")
  configs % posicoes % intervalo = json_get_float_vec(sorteio_json, "posicoes.intervalo")
  
  ! Configuracoes dos momentos
  configs % momentos % distribuicao = json_get_string(sorteio_json, "momentos.distribuicao")
  configs % momentos % regiao = json_get_string(sorteio_json, "momentos.regiao")
  configs % momentos % intervalo = json_get_float_vec(sorteio_json, "momentos.intervalo")
END SUBROUTINE

! ************************************************************
!! Metodo principal
!
! Objetivos:
!   Aplica o sorteio e faz a simulacao.
!
! Modificado:
!   05 de janeiro de 2026
!
! Autoria:
!   oap
! 
SUBROUTINE simular_sorteio (arquivo, out_dir, out_ext)
  CHARACTER(LEN=*), INTENT(IN) :: arquivo
  CHARACTER(LEN=*), INTENT(IN) :: out_dir, out_ext
  TYPE(json_value), POINTER :: infos
  ! Vetores
  REAL(pf), allocatable :: massas(:), posicoes(:,:), momentos(:,:)
  TYPE(sorteio), POINTER :: sorteio_infos

  ALLOCATE(sorteio_infos)

  ! Le o arquivo de configuracoes
  CALL ler_json(arquivo, infos)

  ! Gera os valores iniciais e faz o seu condicionamento
  CALL configuracoes_sorteio(infos, sorteio_infos)
  CALL gerar_condicionar(sorteio_infos, massas, posicoes, momentos)

  ! Roda a simulacao no intervalo [t0, tf]
  CALL rodar_simulacao(out_dir, out_ext, infos, massas, posicoes, momentos)

END SUBROUTINE simular_sorteio

! ************************************************************
!! Sorteia e salva
!
! Objetivos:
!   Sorteia os valores iniciais e salva com um preset de 
!   valores iniciais.
!
! Modificado:
!   16 de dezembro de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE sorteio_salvar (arquivo_in, out_dir)
  IMPLICIT NONE
  CHARACTER(LEN=*) :: arquivo_in
  CHARACTER(LEN=*), INTENT(IN) :: out_dir
  REAL(pf), ALLOCATABLE :: massas(:), posicoes(:,:), momentos(:,:)
  TYPE(json_value), POINTER :: infos
  TYPE(sorteio), POINTER :: sorteio_infos

  ALLOCATE(sorteio_infos)

  ! Le o arquivo de configuracoes
  CALL ler_json(arquivo_in, infos)

  ! Gera os valores iniciais e faz o seu condicionamento
  CALL configuracoes_sorteio(infos, sorteio_infos)
  CALL gerar_condicionar(sorteio_infos, massas, posicoes, momentos)

  ! Verifica se o diretorio de saida existe
  CALL diretorio_vi(out_dir)

  ! Agora salva
  CALL salvar_auto_vi(out_dir, infos, massas, posicoes, momentos)
END SUBROUTINE sorteio_salvar

END module