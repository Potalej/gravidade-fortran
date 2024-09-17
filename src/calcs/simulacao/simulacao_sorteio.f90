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
  USE leitura
  USE arquivos
  IMPLICIT NONE
  PRIVATE
  PUBLIC simular_sorteio, sorteio_salvar

  ! Instanciamento da classe
  TYPE(simular) :: Simulador
  TYPE(preset_config) :: configs

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
  CHARACTER(256), INTENT(INOUT) :: arquivo

  ! Tempo de execucao
  REAL :: t0, tf
  ! Vetores
  REAL(pf), allocatable :: massas(:), posicoes(:,:), momentos(:,:)
  ! Quantidade total de passos
  INTEGER :: qntd_total_passos
  CHARACTER(8) :: datahoje
  CHARACTER(3) :: numero
  CHARACTER(16) :: nome_arq
  CHARACTER(12) :: nome_sorteio
  INTEGER :: i = 1
  LOGICAL :: arquivo_existe = .TRUE.

  ! Le o arquivo de configuracoes
  CALL configs % config(arquivo)

  ! Verifica o modo de sorteio
  IF (configs%modo == "sorteio_henon") THEN
    ! Gera conforme as configuracoes de Henon
    CALL gerar_condicionado_henon(configs%N, &
      massas,    &
      posicoes,  &
      momentos,  &
      configs%int_posicoes, & ! Intervalo de posicoes
      configs%int_momentos, & ! Intervalo de momentos
      configs%int_massas)     ! Intervalo de massas 
  ELSE
    ! Gera os valores a partir das integrais primeiras
    CALL gerar_condicionado_ip(configs%G, &
        configs%N, &
        massas,    &
        posicoes,  &
        momentos,  &
        configs%int_posicoes, & ! Intervalo de posicoes
        configs%int_momentos, & ! Intervalo de momentos
        configs%int_massas,   & ! Intervalo de massas 
        configs%Etot,         & ! Energia total
        configs%Jtot,         & ! Momento angular total
        configs%Ptot )          ! Momento angular total
  ENDIF

  CALL diretorio_out()

  ! Se o instante inicial for negativo, entao vai rodar ao contrario
  IF (configs%t0 < 0) THEN
    IF (configs%tf == 0) THEN
      ! Roda apenas o passado
      WRITE (*,*) " > Intervalo [", configs%t0, ",", configs%tf, "]"
      CALL rodar(-configs%timestep,massas,posicoes,momentos,configs%t0,0)
    ELSE IF (configs%tf > 0) THEN
      ! Roda o passado e o futuro
      WRITE (*,*) " > Intervalo [", configs%t0, ",", 0, "]"
      CALL rodar(-configs%timestep,massas,posicoes,momentos,configs%t0,0)
      WRITE (*,*) " > Intervalo [", 0, ",", configs%tf, "]"
      CALL rodar(configs%timestep,massas,posicoes,momentos,0,configs%tf)
    ENDIF
  ! Se for positivo, apenas roda normal
  ELSE
    ! Roda apenas o futuro
    WRITE (*,*) " > Intervalo [", 0, ",", configs%tf, "]"
    CALL rodar(configs%timestep,massas,posicoes,momentos,0,configs%tf)
  ENDIF

END SUBROUTINE simular_sorteio

! ************************************************************
!! Roda simulacao
!
! Objetivos:
!   Roda uma simulacao com as condicoes informadas.
!
! Modificado:
!   26 de maio de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE rodar (timestep, massas, posicoes, momentos, tempo_inicial, tempo_final)
  REAL(pf), allocatable :: massas(:), posicoes(:,:), momentos(:,:)
  REAL(pf), INTENT(IN)  :: timestep
  REAL(pf)              :: t0, tf
  INTEGER               :: qntd_total_passos
  INTEGER, INTENT(IN)   :: tempo_inicial, tempo_final

  qntd_total_passos = (configs%tf - configs%t0) / configs%timestep
  WRITE (*,*)

  ! timer
  t0 = omp_get_wtime()
  
  CALL Simulador%Iniciar(configs, massas, posicoes, momentos, timestep)

  CALL Simulador%rodar(configs % tf - configs % t0)

  tf = omp_get_wtime()
  WRITE (*,*) ' * tempo ', configs%integrador, ': ', tf - t0
  WRITE (*,*) ' * tempor por passo: ', (tf-t0) / qntd_total_passos
END SUBROUTINE rodar

! ************************************************************
!! Sorteia e salva
!
! Objetivos:
!   Sorteia os valores iniciais e salva com um preset de 
!   valores iniciais.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE sorteio_salvar (dir)
  IMPLICIT NONE
  CHARACTER(LEN=*) :: dir
  REAL(pf), ALLOCATABLE :: massas(:), posicoes(:,:), momentos(:,:)
  CHARACTER(8) :: datahoje
  CHARACTER(3) :: numero
  CHARACTER(16) :: nome_arq
  CHARACTER(12) :: nome_sorteio
  INTEGER :: i = 1
  LOGICAL :: arquivo_existe = .TRUE.

  ! em string
  CALL DATE_AND_TIME(datahoje)

  ! Le o arquivo de configuracoes
  CALL configs % config(dir)

  ! Verifica o modo de sorteio
  IF (configs%modo == "sorteio_henon") THEN
    ! Gera conforme as configuracoes de Henon
    CALL gerar_condicionado_henon(configs%N, &
      massas,    &
      posicoes,  &
      momentos,  &
      configs%int_posicoes, & ! Intervalo de posicoes
      configs%int_momentos, & ! Intervalo de momentos
      configs%int_massas)     ! Intervalo de massas 
  ELSE
    ! Gera os valores a partir das integrais primeiras
    CALL gerar_condicionado_ip(configs%G, &
        configs%N, &
        massas,    &
        posicoes,  &
        momentos,  &
        configs%int_posicoes, & ! Intervalo de posicoes
        configs%int_momentos, & ! Intervalo de momentos
        configs%int_massas,   & ! Intervalo de massas 
        configs%Etot,         & ! Energia total
        configs%Jtot,         & ! Momento angular total
        configs%Ptot )          ! Momento angular total
  ENDIF

  ! Gera o nome
  DO WHILE (arquivo_existe)
    WRITE(numero, '(I3.3)') i
    i = i + 1

    ! cria nome 
    nome_sorteio = TRIM(datahoje)//"_"//TRIM(numero)
    nome_arq = nome_sorteio // ".txt"

    ! verifica se existe
    INQUIRE(file="./out/auto_vi/"//nome_arq, exist=arquivo_existe)
  END DO

  ! Salva o preset gerado
  CALL salvar_sorteio('out/','auto_vi/', nome_arq,  &
    "Sorteio_"//nome_sorteio, &
    configs % G,          &
    massas,               &
    posicoes,             &
    momentos,             &
    configs % t0,         &
    configs % tf,         &
    configs % timestep,   &
    configs % potsoft,   &
    configs % integrador, &
    configs % corretor,   &
    configs % corretor_margem_erro, &
    configs % corretor_max_num_tentativas, &
    configs % colisoes,   &
    configs % colisoes_max_distancia,   &
    configs % passos_antes_salvar, &
    configs % paralelo)
END SUBROUTINE sorteio_salvar

END module