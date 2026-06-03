MODULE testes_auxiliares_mod

USE tipos
USE arquivos_mod
USE utilidades

IMPLICIT NONE
PRIVATE
PUBLIC lista_integradores, &
       limpar_out, &
       atualizar_json_integracao, &
       atualizar_json_sorteio, &
       informacoes_type, &
       get_informacoes

TYPE :: informacoes_type
    ! integracao
    REAL(pf) :: t0, tf, timestep, amortecedor
    ! corpos
    INTEGER :: N
    REAL(pf), ALLOCATABLE :: m(:), qs0(:,:), ps0(:,:), qs(:,:), ps(:,:)
    ! integrais
    REAL(pf) :: E0, J0(3), P0(3), qcm0(3)
    REAL(pf) :: E, J(3), P(3), qcm(3)
END TYPE

CONTAINS

! Apaga o diretorio "testes/data
SUBROUTINE limpar_out ()
    CALL EXECUTE_COMMAND_LINE("rm -rf testes/data", .TRUE.)
END SUBROUTINE

! Lista de integradores
SUBROUTINE lista_integradores (integradores)
  CHARACTER(LEN=:), ALLOCATABLE, INTENT(INOUT) :: integradores(:)
  ALLOCATE(CHARACTER(LEN=20) :: integradores(21))
  ! runge-kuttas
  integradores(1)  = "euler_exp"
  integradores(2)  = "euler_imp"
  integradores(3)  = "rungekutta2"
  integradores(4)  = "rungekutta3"
  integradores(5)  = "rungekutta4"
  ! simpleticos - basicos
  integradores(6)  = "euler_simp"
  integradores(7)  = "verlet"
  ! simpleticos - metodos de Ruth
  integradores(8)  = "ruth3"
  integradores(9)  = "ruth4"
  ! simpleticos - composicao de metodos de Euler
  integradores(10) = "ecp4s5"
  integradores(11) = "ecp4s6"
  ! simpleticos - Runge-Kutta-Nystrom
  integradores(12) = "rkn551"
  integradores(13) = "rkn671"
  ! simpleticos - composicao de metodos de Verlet
  integradores(14) = "svcp6s9"
  integradores(15) = "svcp8s15"
  integradores(16) = "svcp8s17"
  integradores(17) = "svcp10s35"
  ! multipasso lineares
  integradores(18) = "ab2"
  integradores(19) = "ab3"
  integradores(20) = "ab4"
  integradores(21) = "ab5"
END SUBROUTINE lista_integradores

! Atualiza um JSON (integracao)
SUBROUTINE atualizar_json_integracao (infos, chave, valor)
  TYPE(json_value), POINTER, INTENT(INOUT) :: infos
  CHARACTER(LEN=*), INTENT(IN) :: chave, valor
  TYPE(json_value), POINTER    :: subobj, integracao
  LOGICAL :: encontrado

  CALL json % get(infos, TRIM("integracao."//chave), subobj, encontrado)
  IF (encontrado) CALL json % remove(subobj)

  CALL json % get(infos, "integracao", integracao, encontrado)
!   CALL json % get(infos, "valores_iniciais", subobj, encontrado)
  CALL json % add(integracao, TRIM(chave), TRIM(valor))
END SUBROUTINE

! Atualiza um JSON (sorteio)
SUBROUTINE atualizar_json_sorteio (infos, chave, valor)
  TYPE(json_value), POINTER, INTENT(INOUT) :: infos
  CHARACTER(LEN=*), INTENT(IN) :: chave, valor
  TYPE(json_value), POINTER :: subobj, sorteio
  TYPE(json_value), POINTER :: p_m, p_q, p_p
  LOGICAL :: encontrado

  CALL json % get(infos, TRIM("sorteio."//chave), subobj, encontrado)
  IF (encontrado) CALL json % remove(subobj)

  CALL json % get(infos, "sorteio", sorteio, encontrado)
  CALL json % add(sorteio, TRIM(chave), TRIM(valor))
END SUBROUTINE

! Monta um vetor de informacoes
SUBROUTINE get_informacoes (ijson, it)
  TYPE(json_value), POINTER, INTENT(IN) :: ijson
  TYPE(informacoes_type), INTENT(INOUT) :: it
  REAL(pf), ALLOCATABLE :: q3(:), p3(:)
  INTEGER :: p
  CHARACTER(32) :: p_string
  
  CALL json % get(ijson, 'integracao.t0', it % t0)
  CALL json % get(ijson, 'integracao.tf', it % tf)
  it % timestep = json_get_float(ijson, 'integracao.timestep')
  CALL json % get(ijson, 'integracao.amortecedor', it % amortecedor)
  CALL json % get(ijson, 'N', it % N)
  it % m = json_get_float_vec(ijson, 'valores_iniciais.massas')

  ALLOCATE(it % qs0(it % N,3))
  ALLOCATE(it % ps0(it % N,3))
  ALLOCATE(it % qs(it % N,3))
  ALLOCATE(it % ps(it % N,3))
  DO p = 1, it % N
    WRITE(p_string, *) p
    q3 = json_get_float_vec(ijson, 'valores_iniciais.posicoes['//p_string//']')
    p3 = json_get_float_vec(ijson, 'valores_iniciais.momentos['//p_string//']')
    it % qs0(p,:) = q3
    it % ps0(p,:) = p3
  END DO
  it % qs = it % qs0
  it % ps = it % ps0
END SUBROUTINE

END MODULE