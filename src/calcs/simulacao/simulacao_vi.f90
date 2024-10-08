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
  TYPE(simular) :: Simulador
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

  CALL Simulador%Iniciar(configs, massas, posicoes, momentos, timestep)

  CALL Simulador%rodar(configs % tf - configs % t0)

  tf = omp_get_wtime()
  WRITE (*,*) ' * tempo ', configs%integrador, ': ', tf - t0
  WRITE (*,*) ' * tempor por passo: ', (tf-t0) / qntd_total_passos
END SUBROUTINE rodar

END MODULE simulacao_vi