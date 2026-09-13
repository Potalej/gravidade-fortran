! ************************************************************
!! SIMULACAO: SORTEIO
!
! Objetivos:
!   Simulacoes a partir do sorteio de valores iniciais.
!
! Modificado:
!   13 de setembro de 2026
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
!! Validacao dos parametros do arquivo de entrada de sorteio
!
! Modificado:
!   13 de setembro de 2026
!
! Autoria:
!   oap
! 
SUBROUTINE parametros (infos, integracao)
  TYPE(json_value), POINTER, INTENT(IN) :: infos
  LOGICAL, INTENT(IN) :: integracao
  CHARACTER(LEN=40)   :: pars(39)

  pars(:) = ""

  ! gerais
  pars(1) = "modo"
  pars(2) = "" ! nome
  pars(3) = "N"
  pars(4) = "G"

  IF (integracao) THEN
    pars(5) = "paralelo"
    pars(6) = "gpu"
    pars(7) = "massas_iguais"
    pars(8) = "exibir"
  ENDIF

  ! sorteio
  pars(9)  = "sorteio.integrais"
  pars(10) = "sorteio.integrais.energia_total"
  pars(11) = "sorteio.integrais.angular_total"
  pars(12) = "sorteio.integrais.linear_total"
  pars(13) = "sorteio.massas"
  pars(14) = "sorteio.massas.normalizadas"
  pars(15) = "sorteio.massas.intervalo"
  pars(16) = "sorteio.massas.distribuicao"
  pars(17) = "sorteio.posicoes"
  pars(18) = "sorteio.posicoes.intervalo"
  pars(19) = "sorteio.posicoes.distribuicao"
  pars(20) = "sorteio.posicoes.regiao"
  pars(21) = "sorteio.momentos"
  pars(22) = "sorteio.momentos.intervalo"
  pars(23) = "sorteio.momentos.distribuicao"
  pars(24) = "sorteio.momentos.regiao"

  IF (integracao) THEN
    ! integracao
    pars(25) = "integracao.metodo"
    pars(26) = "integracao.timestep"
    pars(27) = "integracao.amortecedor"
    pars(28) = "integracao.t0"
    pars(29) = "integracao.tf"
    pars(30) = "integracao.checkpoints"
    pars(31) = "integracao.tree"
    pars(32) = "integracao.theta"

    ! colisoes
    pars(33) = "colisoes.colidir"
    pars(34) = "colisoes.metodo"
    pars(35) = "colisoes.densidade"
    pars(36) = "colisoes.permitir_choque_inicial"

    ! correcao
    pars(37) = "correcao.corrigir"
    pars(38) = "correcao.margem_erro"
    pars(39) = "correcao.max_num_tentativas"
  ENDIF

  CALL validar_json(infos, pars)
END SUBROUTINE

! ************************************************************
!! Geracao de configuracoes a partir do json_value
!
! Objetivos:
!   Converte json_value em objetos sorteio e sorteio_vetores.
!
! Modificado:
!   29 de maio de 2026
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
!   13 de setembro de 2026
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

  ! Valida o arquivo de entrada
  CALL parametros(infos, .TRUE.)

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
!   13 de setembro de 2026
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

  ! Valida o arquivo de entrada
  CALL parametros(infos, .FALSE.)

  ! Gera os valores iniciais e faz o seu condicionamento
  CALL configuracoes_sorteio(infos, sorteio_infos)
  CALL gerar_condicionar(sorteio_infos, massas, posicoes, momentos)

  ! Verifica se o diretorio de saida existe
  CALL diretorio_vi(out_dir)

  ! Agora salva
  CALL salvar_auto_vi(out_dir, infos, massas, posicoes, momentos)
END SUBROUTINE sorteio_salvar

END module