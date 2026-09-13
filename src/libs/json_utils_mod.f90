MODULE json_utils_mod
    USE json_module, only: json_core, json_value, json_ck
    USE tipos
    IMPLICIT NONE

    TYPE(json_core) :: json
    PUBLIC
CONTAINS

LOGICAL FUNCTION json_chave_existe (dados, chave)
    TYPE(json_core) :: json
    TYPE(json_value), POINTER, INTENT(IN) :: dados
    CHARACTER(LEN=*), INTENT(IN) :: chave

    CALL json % info(dados, chave, found=json_chave_existe)
END FUNCTION

FUNCTION json_get_string (dados, chave, encontrado_par) RESULT(text)
    TYPE(json_core) :: json
    TYPE(json_value), POINTER, INTENT(IN) :: dados
    CHARACTER(LEN=*), INTENT(IN) :: chave
    CHARACTER(LEN=:, KIND=json_ck), ALLOCATABLE :: text_ck
    CHARACTER(LEN=:), ALLOCATABLE :: text
    LOGICAL, INTENT(INOUT), OPTIONAL :: encontrado_par
    LOGICAL :: encontrado
    INTEGER :: i

    CALL json % get(dados, chave, text_ck, encontrado)
    IF (PRESENT(encontrado_par)) encontrado_par = encontrado

    IF (.NOT. encontrado) RETURN

    ALLOCATE(CHARACTER(LEN=LEN(text_ck)) :: text)
    DO i = 1, LEN(text)
        text(i:i) = ACHAR(IACHAR(text_ck(i:i)))
    END DO
END FUNCTION json_get_string

FUNCTION json_get_int (dados, chave, encontrado_par) RESULT(valor)
    TYPE(json_core) :: json
    TYPE(json_value), POINTER, INTENT(IN) :: dados
    CHARACTER(LEN=*), INTENT(IN) :: chave
    LOGICAL, INTENT(INOUT), OPTIONAL :: encontrado_par
    INTEGER :: valor
    LOGICAL :: encontrado
    CALL json % get(dados, chave, valor, encontrado)
    IF (PRESENT(encontrado_par)) encontrado_par = encontrado
END FUNCTION

FUNCTION json_get_logical (dados, chave, encontrado_par) RESULT(valor)
    TYPE(json_core) :: json
    TYPE(json_value), POINTER, INTENT(IN) :: dados
    CHARACTER(LEN=*), INTENT(IN) :: chave
    LOGICAL, INTENT(INOUT), OPTIONAL :: encontrado_par
    LOGICAL :: valor
    LOGICAL :: encontrado
    CALL json % get(dados, chave, valor, encontrado)
    IF (PRESENT(encontrado_par)) encontrado_par = encontrado
END FUNCTION

FUNCTION json_get_float (dados, chave, encontrado_par) RESULT(valor)
    TYPE(json_core) :: json
    TYPE(json_value), POINTER, INTENT(IN) :: dados
    CHARACTER(LEN=*), INTENT(IN) :: chave
    LOGICAL, INTENT(INOUT), OPTIONAL :: encontrado_par
    REAL(pf)   :: valor
    LOGICAL :: encontrado
    CALL json % get(dados, chave, valor, encontrado)
    IF (PRESENT(encontrado_par)) encontrado_par = encontrado
END FUNCTION json_get_float

FUNCTION json_get_float_vec (dados, chave, encontrado_par) RESULT(valor)
    TYPE(json_core) :: json
    TYPE(json_value), POINTER, INTENT(IN) :: dados
    CHARACTER(LEN=*), INTENT(IN) :: chave
    LOGICAL, INTENT(INOUT), OPTIONAL :: encontrado_par
    REAL(pf), ALLOCATABLE   :: valor(:)
    LOGICAL :: encontrado
    CALL json % get(dados, chave, valor, encontrado)
    IF (PRESENT(encontrado_par)) encontrado_par = encontrado
END FUNCTION json_get_float_vec

FUNCTION json_get_float_matrix (dados, chave, linhas, colunas, encontrado_par) RESULT(valor)
    TYPE(json_core) :: json
    TYPE(json_value), POINTER, INTENT(IN) :: dados
    CHARACTER(LEN=*), INTENT(IN) :: chave
    INTEGER, INTENT(IN) :: linhas, colunas
    LOGICAL, INTENT(INOUT), OPTIONAL :: encontrado_par
    INTEGER       :: linha
    CHARACTER(32) :: linha_str
    REAL(pf), ALLOCATABLE :: valor_linha(:), valor(:,:)
    LOGICAL :: encontrado

    ALLOCATE(valor(linhas, colunas))
    ALLOCATE(valor_linha(colunas))
    DO linha = 1, linhas
        WRITE(linha_str, *) linha
        valor_linha = json_get_float_vec(dados, chave//"["//linha_str//"]", encontrado)
        IF (.NOT. encontrado) THEN
            IF (PRESENT(encontrado_par)) encontrado_par = encontrado
            RETURN
        ENDIF
        valor(linha,:) = valor_linha
    END DO
    DEALLOCATE(valor_linha)
END FUNCTION json_get_float_matrix

SUBROUTINE json_clone (entrada, saida)
    TYPE(json_core) :: json
    TYPE(json_value), POINTER, INTENT(IN) :: entrada
    TYPE(json_value), POINTER, INTENT(OUT) :: saida
    call json % clone(entrada, saida)
END SUBROUTINE json_clone

END MODULE json_utils_mod