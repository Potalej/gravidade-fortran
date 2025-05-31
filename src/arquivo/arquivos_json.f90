! *****************************************************************
!! ARQUIVOS DE SORTEIOS
!
! Objetivos:
!   Rotinas para manipular arquivos de sorteios. Especificamente, a
!   leitura dos arquivos de pre-valores para condicionamento.
!   
! Modificado:
!   01 de maio de 2025
! 
! Autoria:
!   oap
!
MODULE arquivos_json
    USE tipos
    USE json_utils_mod

    IMPLICIT NONE
    CHARACTER(6)  :: DIR_OUT  = "./out/"
    CHARACTER(11) :: DIR_DATA = "./out/data/"
    CHARACTER(23) :: DIR_AVI  = "./out/valores_iniciais/"

CONTAINS

! ************************************************************
!! Leitura de valores de sorteio
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE ler_json (arquivo, dados)
    CHARACTER(LEN=*), INTENT(IN) :: arquivo
    TYPE(json_value), POINTER, INTENT(INOUT) :: dados

    WRITE (*,'(A)') 'PRESET: '//arquivo
    CALL json % parse(file=arquivo, p=dados)
END SUBROUTINE ler_json

! ************************************************************
!! Gera um nome de arquivo em um diretorio e com uma extensao
!
! Modificado:
!   05 de maio de 2025
!
! Autoria:
!   oap
! 
FUNCTION gerar_nome_arquivo (diretorio, extensao) RESULT(arquivo_nome)
    CHARACTER(LEN=*) :: diretorio, extensao
    CHARACTER(8)  :: data_hoje
    CHARACTER(3)  :: numero
    CHARACTER(12) :: data_numero
    CHARACTER(:), ALLOCATABLE :: arquivo_nome
    INTEGER       :: i
    LOGICAL       :: arquivo_existe

    ALLOCATE(CHARACTER(13+LEN(extensao)) :: arquivo_nome)

    CALL DATE_AND_TIME(data_hoje)
    arquivo_existe = .TRUE.
    i = 1
    DO WHILE (arquivo_existe)
        WRITE(numero, '(I3.3)') i
        i = i + 1

        data_numero = TRIM(data_hoje)//"_"//TRIM(numero)
        arquivo_nome = data_numero//"."//extensao

        ! Verifica se existe
        INQUIRE(file=diretorio//arquivo_nome, exist=arquivo_existe)
    END DO
END FUNCTION gerar_nome_arquivo

! ************************************************************
!! Salva arquivo de valores iniciais no diretorio auto_vi
!
! Modificado:
!   05 de maio de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE salvar_auto_vi (infos_sorteio, massas, posicoes, momentos)
    TYPE(json_value), POINTER   :: infos_sorteio
    REAL(pf)                    :: massas(:), posicoes(:,:), momentos(:,:)

    CALL salvar_vi_json(DIR_AVI, infos_sorteio, massas, posicoes, momentos, .TRUE.)
END SUBROUTINE salvar_auto_vi

! ************************************************************
!! Salva arquivo de valores iniciais em algum diretorio
!
! Modificado:
!   05 de maio de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE salvar_vi_json (diretorio, infos_sorteio, massas, posicoes, momentos, gerar_nome)
    CHARACTER(LEN=*), INTENT(IN) :: diretorio
    TYPE(json_value), POINTER :: infos_sorteio, dummy
    REAL(pf)                  :: massas(:), posicoes(:,:), momentos(:,:)
    REAL(pf64), ALLOCATABLE   :: massas64(:), posicoes64(:,:), momentos64(:,:)
    TYPE(json_core)           :: json
    TYPE(json_value), POINTER :: modo, vi, ap, am
    LOGICAL, OPTIONAL :: gerar_nome
    INTEGER :: i
    LOGICAL :: arquivo_existe, vi_existe
    CHARACTER(8) :: datahoje
    CHARACTER(3) :: numero
    CHARACTER(17) :: nome_arq
    CHARACTER(12) :: nome_sorteio

    ! Verificando nome para arquivo de saida
    IF (PRESENT(gerar_nome) .AND. gerar_nome) THEN
        nome_arq = gerar_nome_arquivo(diretorio, "json")
        WRITE(*,*) 'nome do arquivo:' // nome_arq
    ELSE
        nome_arq = ".json"
    ENDIF

    ! Altera o modo se for o caso
    CALL json % get(infos_sorteio, 'modo', modo)
    CALL json % remove(modo)
    CALL json % add(infos_sorteio, 'modo', 'vi')

    ! Por fim, adiciona os valores iniciais na nova estrutura se necessario
    CALL json % get(infos_sorteio, 'valores_iniciais', dummy, vi_existe)

    IF (.NOT. vi_existe) THEN
        CALL json % create_object(vi, 'valores_iniciais')
        CALL json % add(infos_sorteio, vi)

        ! Adiciona as massas
        massas64 = REAL(massas, KIND=pf64)
        CALL json % add(vi, 'massas', massas64)
        
        ! Adiciona as posicoes e os momentos
        CALL json % create_array(ap, 'posicoes')
        CALL json % create_array(am, 'momentos')
        CALL json % add(vi, ap)
        CALL json % add(vi, am)
        
        posicoes64 = REAL(posicoes, KIND=pf64)
        momentos64 = REAL(momentos, KIND=pf64)

        DO i=1, SIZE(massas)
            CALL json % add(ap, 'posicoes', posicoes64(i,:))
            CALL json % add(am, 'momentos', momentos64(i,:))
        END DO   
    ENDIF

    ! Salva o arquivo
    CALL json % print(infos_sorteio, TRIM(diretorio//nome_arq))
    
    ! Limpa a memoria
    NULLIFY(modo, vi, ap, am)
    CALL json % destroy(infos_sorteio)
END SUBROUTINE salvar_vi_json

END MODULE arquivos_json