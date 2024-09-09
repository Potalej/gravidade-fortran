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

  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  USE OMP_LIB
  USE simulacao
  USE leitura
  IMPLICIT NONE
  PRIVATE
  PUBLIC simular_vi

  ! Instanciamento da classe
  TYPE(simular) :: Sim_rk4, Sim_verlet, Sim_corrigir
  TYPE(preset_config) :: configs

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

  ! Tempo de execucao
  REAL :: t0, tf
  ! Quantidade total de passos
  INTEGER :: qntd_total_passos

  ! Le o arquivo de configuracoes
  CALL configs % valores_iniciais(arquivo)

  ! Se o instante inicial for negativo, entao vai rodar ao contrario
  IF (configs%t0 < 0) THEN
    IF (configs%tf == 0) THEN
      ! Roda apenas o passado
      WRITE (*,*) " * Intervalo [", configs%t0, ",", configs%tf, "]"
      CALL rodar(-configs%timestep,configs%massas,configs%R,configs%P)
    ELSE IF (configs%tf > 0) THEN
      ! Roda o passado e o futuro
      WRITE (*,*) " * Intervalo [", configs%t0, ",", 0, "]"
      CALL rodar(-configs%timestep,configs%massas,configs%R,configs%P)
      WRITE (*,*) " * Intervalo [", 0, ",", configs%tf, "]"
      CALL rodar(configs%timestep,configs%massas,configs%R,configs%P)
    ENDIF
  ! Se for positivo, apenas roda normal
  ELSE
    ! Roda apenas o futuro
    WRITE (*,*) " * Intervalo [", 0, ",", configs%tf, "]"
    CALL rodar(configs%timestep,configs%massas,configs%R,configs%P)
  ENDIF

END SUBROUTINE simular_vi

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
SUBROUTINE rodar (timestep, massas, posicoes, momentos)
  REAL(pf), allocatable :: massas(:), posicoes(:,:), momentos(:,:)
  REAL(pf), INTENT(IN)  :: timestep
  REAL(pf)              :: t0, tf
  INTEGER               :: qntd_total_passos

  qntd_total_passos = (configs%tf - configs%t0) / configs%timestep
  ! timer
  t0 = omp_get_wtime()

  WRITE (*,*)

  SELECT CASE (configs%integrador)
    CASE ("verlet")
      Sim_verlet % corrigir = configs%corretor
      Sim_verlet % corrigir_margem_erro = configs%corretor_margem_erro
      Sim_verlet % corrigir_max_num_tentativas = configs%corretor_max_num_tentativas
      
      Sim_verlet % colidir  = configs%colisoes
      Sim_verlet % colisoes_max_distancia = configs%colisoes_max_distancia
      
      CALL Sim_verlet%Iniciar(configs%G, massas, posicoes, momentos, &
        timestep, configs%potsoft, configs%passos_antes_salvar,configs%integrador)
      CALL Sim_verlet%rodar_verlet(configs % tf - configs % t0)
    CASE ("rk4")
      Sim_rk4 % corrigir = configs%corretor
      Sim_rk4 % corrigir_margem_erro = configs%corretor_margem_erro
      Sim_rk4 % corrigir_max_num_tentativas = configs%corretor_max_num_tentativas

      Sim_rk4 % colidir  = configs%colisoes
      Sim_rk4 % colisoes_max_distancia = configs%colisoes_max_distancia

      CALL Sim_rk4%Iniciar(configs%G, massas, posicoes, momentos, &
        timestep, configs%potsoft, configs%passos_antes_salvar,configs%integrador)
      CALL Sim_rk4%rodar_rk4(configs % tf - configs % t0)
  END SELECT

  tf = omp_get_wtime()
  WRITE (*,*) ' * tempo ', configs%integrador, ': ', tf - t0
  WRITE (*,*) ' * tempor por passo: ', (tf-t0) / qntd_total_passos
END SUBROUTINE rodar

END MODULE simulacao_vi