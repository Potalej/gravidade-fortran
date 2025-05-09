MODULE json_utils_mod
    USE json_module, only: json_core, json_value, json_ck
    IMPLICIT NONE

    TYPE(json_core) :: json
    PUBLIC :: json, json_core, json_value, json_get_string
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

END MODULE json_utils_mod