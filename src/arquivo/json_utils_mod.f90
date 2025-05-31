MODULE json_utils_mod
    USE json_module, only: json_core, json_value, json_ck
    USE tipos
    IMPLICIT NONE

    TYPE(json_core) :: json
    PUBLIC
CONTAINS

FUNCTION json_get_string (dados, chave) RESULT(texto)
    TYPE(json_core) :: json
    TYPE(json_value), POINTER, INTENT(IN) :: dados
    CHARACTER(LEN=*), INTENT(IN) :: chave
    CHARACTER(LEN=:, KIND=json_ck), ALLOCATABLE :: texto_ck
    CHARACTER(LEN=:), ALLOCATABLE :: texto
    INTEGER :: i

    CALL json % get(dados, chave, texto_ck)

    ALLOCATE(CHARACTER(LEN=LEN(texto_ck)) :: texto)
    DO i = 1, LEN(texto)
        texto(i:i) = ACHAR(IACHAR(texto_ck(i:i)))
    END DO
END FUNCTION json_get_string

FUNCTION json_get_float (dados, chave) RESULT(valor)
    TYPE(json_core) :: json
    TYPE(json_value), POINTER, INTENT(IN) :: dados
    CHARACTER(LEN=*), INTENT(IN) :: chave
    REAL(pf64) :: valor64
    REAL(pf)   :: valor

    CALL json % get(dados, chave, valor64)
    valor = REAL(valor64, KIND=pf)
END FUNCTION json_get_float

FUNCTION json_get_float_vec (dados, chave) RESULT(valor)
    TYPE(json_core) :: json
    TYPE(json_value), POINTER, INTENT(IN) :: dados
    CHARACTER(LEN=*), INTENT(IN) :: chave
    REAL(pf64), ALLOCATABLE :: valor64(:)
    REAL(pf), ALLOCATABLE   :: valor(:)

    CALL json % get(dados, chave, valor64)
    valor = REAL(valor64, KIND=pf)
END FUNCTION json_get_float_vec

END MODULE json_utils_mod