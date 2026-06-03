! ************************************************************
!! Modulo de testes
!
! Objetivos:
!   Contem alguns testes que sao rodados e exibidos na tela,
!   so para garantir que o programa esta funcionando.
!
! Modificado:
!   03 de junho de 2026
!
! Autoria:
!   oap
!  
MODULE testes_mod

!> Versionamento
  USE version
  USE iso_fortran_env, only: output_unit
!> Tipos  
  USE tipos
  USE simulacao
  USE arquivos_mod
  USE utilidades
  USE testes_auxiliares_mod
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC testar

CONTAINS

SUBROUTINE print_espacador (arquivo)
  INTEGER :: arquivo
  WRITE (arquivo, *)
  WRITE (arquivo, *) "==============================================="
  WRITE (arquivo, *)
END SUBROUTINE

SUBROUTINE print_titulo (titulo, arquivo)
  INTEGER :: arquivo
  CHARACTER(LEN=*), INTENT(IN) :: titulo

  WRITE (arquivo, *)
  WRITE (arquivo, *) "==============================================="
  WRITE (arquivo, *) "# " // TRIM(titulo)
  WRITE (arquivo, *) "==============================================="
  WRITE (arquivo, *)
END SUBROUTINE

SUBROUTINE testar ()
  INTEGER :: arquivo = output_unit
  CHARACTER(LEN=20) :: data_hoje

  IF (usar_gpu) THEN
    WRITE (arquivo,*) 'v', version_string, '_', precisao, '_GPU (', build_date, ' ', build_time, ')'
  ELSE
    WRITE (arquivo,*) 'v', version_string, '_', precisao, ' (', build_date, ' ', build_time, ')'
  ENDIF
  CALL data_hora_string(data_hoje)
  WRITE (arquivo,*) data_hoje

  CALL limpar_out()

  !# Teste dos metodos
  !> Validacao: Todos devem rodar corretamente
  CALL print_titulo("TESTE 1: USO DOS INTEGRADORES",arquivo)
  WRITE(arquivo,*) "> Descricao:"
  WRITE(arquivo,*) "  Todos os integradores devem rodar corretamente usando massas iguais"
  WRITE(arquivo,*) "  e usando massas diferentes. Se rodou, entao esta ok."
  WRITE(arquivo,*)
  CALL TESTE_metodos(arquivo)
  CALL print_espacador(arquivo)

  CALL limpar_out()
  
  !# Teste de ordem de convergencia
  !> Validacao: Todos devem ter a ordem de convergencia esperada (se tiver pontos flutuantes suficientes)
  CALL print_titulo("TESTE 2: ORDEM DE CONVERGENCIA DOS INTEGRADORES",arquivo)
  WRITE(arquivo,*) "> Descricao:"
  WRITE(arquivo,*) "  A ultima coluna deve convergir para a ordem do metodo. Para os metodos"
  WRITE(arquivo,*) "  de ordem alta pode ser necessario usar precisao quadrupla."
  CALL TESTE_ordem_convergencia(arquivo)
  CALL print_espacador(arquivo)

  CALL limpar_out()
END SUBROUTINE

! Para testar o uso dos metodos
SUBROUTINE TESTE_metodos (arquivo)
  INTEGER :: arquivo
  TYPE(json_value), POINTER :: infos_md, infos_mi
  TYPE(simular),    POINTER :: simulador
  CHARACTER(LEN=:), ALLOCATABLE :: integradores(:)
  TYPE(informacoes_type) :: it_md, it_mi
  INTEGER  :: i, N
  REAL(pf) :: cron0, cronf

  ! Le os arquivo de valores iniciais
  CALL ler_json("testes/presets/vi_metodos_md.json", infos_md, .FALSE.)
  CALL ler_json("testes/presets/vi_metodos_mi.json", infos_mi, .FALSE.)
  CALL lista_integradores(integradores)

  WRITE (arquivo, *) "> Informacoes dos presets:"
  CALL get_informacoes(infos_md, it_md)
  CALL get_informacoes(infos_mi, it_mi)
  
  ! WRITE(*, '(3X,A,2X,E12.4)')
  WRITE (arquivo,'(3X,A,I8)') "- N = ", it_md % N
  WRITE (arquivo,'(3X,A,F12.1,A,F12.1,A)') "- intervalo: [", it_md % t0, ",", it_md % tf, "]"
  WRITE (arquivo,'(3X,A,E12.5)')           "- timestep: ", it_md % timestep
  WRITE (arquivo,'(3X,A,E12.5)')           "- amortecedor: ", it_md % amortecedor

  DO i = 1, SIZE(integradores)
    WRITE (arquivo,*)
    WRITE(arquivo, *) '------------------------------'
    WRITE (arquivo,*)
    WRITE(arquivo, *) '[metodo: ', TRIM(integradores(i)), ']'

    ! Atualiza
    CALL atualizar_json_integracao(infos_md, "metodo", integradores(i))
    CALL atualizar_json_integracao(infos_mi, "metodo", integradores(i))
    
    ! Simulacao (MASSAS DIFERENTES)
    ALLOCATE(simulador)
    cron0 = OMP_GET_WTIME()
    CALL simulador%iniciar(infos_md, it_md % m, it_md % qs0, it_md % ps0, it_md % timestep, &
                            "testes", ".bin", .FALSE.)
    CALL simulador%rodar(INT(it_md % tf - it_md % t0))
    cronf = OMP_GET_WTIME()
    WRITE(arquivo, '(3X,A,2X,E12.4,A)') 'E0 - E', simulador % errene, "[MD]"
    WRITE(arquivo, '(3X,A,2X,F12.3,A)') 'Tempo: ', cronf - cron0, "[MD]"

    NULLIFY(simulador % integrador)
    DEALLOCATE(simulador)

    ! Simulacao (MASSAS IGUAIS)
    ALLOCATE(simulador)
    cron0 = OMP_GET_WTIME()
    CALL simulador%iniciar(infos_mi, it_mi % m, it_mi % qs0, it_mi % ps0, it_mi % timestep, &
                            "testes", ".bin", .FALSE.)
    CALL simulador%rodar(INT(it_mi % tf - it_mi % t0))
    cronf = OMP_GET_WTIME()
    WRITE(arquivo, '(3X,A,2X,E12.4,A)') 'E0 - E', simulador % errene, "[MI]"
    WRITE(arquivo, '(3X,A,2X,F12.3,A)') 'Tempo: ', cronf - cron0, "[MI]"

    NULLIFY(simulador % integrador)
    DEALLOCATE(simulador)
  END DO
END SUBROUTINE

! Para testar a ordem de convergencia
SUBROUTINE TESTE_ordem_convergencia (arquivo)
  INTEGER :: arquivo
  TYPE(json_value), POINTER :: infos
  TYPE(simular), POINTER :: simulador
  CHARACTER(LEN=:), ALLOCATABLE :: integradores(:)
  REAL(pf) :: erro_num, erro_den, dt_local
  INTEGER :: i, teste, corpo, qntd_testes = 6
  CHARACTER(32) :: corpo_string

  TYPE(informacoes_type) :: it
  REAL(pf), ALLOCATABLE :: q(:,:,:), p(:,:,:)

  ! le o arquivo de valores iniciais
  CALL ler_json("testes/presets/vi_ordem.json", infos, .FALSE.)

  ! lista de integradores
  CALL lista_integradores(integradores)

  ! informacoes
  CALL get_informacoes(infos, it)

  ALLOCATE(q(qntd_testes, it % N, 3))
  ALLOCATE(p(qntd_testes, it % N, 3))

  DO i = 1, SIZE(integradores)

    CALL atualizar_json_integracao(infos, "metodo", integradores(i))

    WRITE(arquivo,*)
    WRITE(arquivo, *) '------------------------------'
    WRITE(arquivo,*)
    WRITE(arquivo, *) '[metodo: ', TRIM(integradores(i)), ']'

    dt_local = it % timestep

    DO teste = 1, qntd_testes
      ALLOCATE(simulador)

      ! captura as posicoes e momentos iniciais
      q(teste,:,:) = it % qs0
      p(teste,:,:) = it % ps0

      ! timestep
      dt_local = 0.5_pf * dt_local

      ! agora roda a simulacao
      CALL simulador%iniciar(infos, it%m, q(teste,:,:), p(teste,:,:), dt_local, &
                            "testes", ".bin", .FALSE.)
      CALL simulador%rodar(INT(it % tf - it % t0))
      
      q(teste,:,:) = simulador % R
      p(teste,:,:) = simulador % P

      NULLIFY(simulador % integrador)
      DEALLOCATE(simulador)
    END DO

    ! tendo rodado todos os testes, agora calcula os erros
    DO teste = 2, qntd_testes - 1
      erro_num = NORM2(q(teste-1,:,:) - q(teste,:,:))
      erro_den = NORM2(q(teste,:,:) - q(teste+1,:,:))
      WRITE(arquivo, *) erro_num, erro_den, LOG(erro_num / erro_den)/LOG(2.0_pf)
    END DO

  END DO

  DEALLOCATE(q, p)
END SUBROUTINE

END MODULE