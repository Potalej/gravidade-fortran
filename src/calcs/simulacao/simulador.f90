! ************************************************************
!! Simulador
!
! Objetivos:
!   Rotinas para invocar o simulador mais facilmente.
!
! Modificado:
!   11 de novembro de 2025
!
! Autoria:
!   oap
! 
MODULE simulador_mod
  USE tipos
  USE json_utils_mod
  USE OMP_LIB
  USE simulacao

  IMPLICIT NONE
  PRIVATE
  PUBLIC rodar_simulacao

CONTAINS

SUBROUTINE rodar_simulacao (od, oe, infos, massas, posicoes, momentos)
  CHARACTER(LEN=*), INTENT(IN) :: od, oe ! out_dir, out_ext
  TYPE(json_value), POINTER :: infos
  REAL(pf), INTENT(IN) :: massas(:), posicoes(:,:), momentos(:,:)
  INTEGER  :: t0, tf
  REAL(pf) :: timestep

  CALL json % get(infos, 'integracao.t0', t0)
  CALL json % get(infos, 'integracao.tf', tf)
  timestep = json_get_float(infos, 'integracao.timestep')

  IF (t0 < 0) THEN
    IF (tf == 0) THEN
      ! Roda apenas o passado
      WRITE (*,*) " > Intervalo [", t0, ",", tf, "]"
      CALL rodar_simulacao_intervalo(od, oe, infos, massas, posicoes, momentos, -timestep, t0, 0)

    ELSE IF (tf > 0) THEN     
      ! Roda o passado e o futuro
      WRITE (*,*) " > Intervalo [", t0, ",", 0, "]"
      CALL rodar_simulacao_intervalo(od, oe, infos, massas, posicoes, momentos, -timestep, t0, 0)
      
      WRITE (*,*) " > Intervalo [", 0, ",", tf, "]"
      CALL rodar_simulacao_intervalo(od, oe, infos, massas, posicoes, momentos, timestep, 0, tf)

    ENDIF
  
  ! Se for positivo, apenas roda normal
  ELSE
    ! Roda apenas o futuro
    WRITE (*,*) " > Intervalo [", 0, ",", tf, "]"
    CALL rodar_simulacao_intervalo(od, oe, infos, massas, posicoes, momentos, timestep, 0, tf)
  ENDIF
END SUBROUTINE rodar_simulacao

SUBROUTINE rodar_simulacao_intervalo (od, oe, infos, massas, posicoes, momentos, timestep, t0, tf)
  CHARACTER(LEN=*), INTENT(IN):: od, oe ! out_dir, out_ext
  TYPE(json_value), POINTER :: infos
  REAL(pf), INTENT(IN) :: massas(:), posicoes(:,:), momentos(:,:)
  REAL(pf), INTENT(IN) :: timestep
  INTEGER, INTENT(IN)  :: t0, tf
  REAL(pf) :: timer_0, timer_1
  INTEGER :: qntd_total_passos
  CHARACTER(LEN=:), ALLOCATABLE :: metodo
  
  !> Instanciamento do simulador
  TYPE(simular) :: simulador

  qntd_total_passos = ABS((tf - t0)/timestep)

  ! Timer
  timer_0 = omp_get_wtime()
  
  ! Inicializando a classe de simulacao
  CALL simulador%iniciar(infos, massas, posicoes, momentos, timestep, od, oe)

  ! Agora roda
  CALL simulador%rodar(tf - t0)

  metodo = json_get_string(infos, 'integracao.metodo')

  timer_1 = omp_get_wtime()
  WRITE (*,*) ' * tempo ', metodo, ': ', timer_1 - timer_0
  WRITE (*,*) ' * tempo por passo: ', (timer_1 - timer_0) / qntd_total_passos
END SUBROUTINE rodar_simulacao_intervalo

END MODULE simulador_mod