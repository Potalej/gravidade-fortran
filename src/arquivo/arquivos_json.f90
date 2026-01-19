! *****************************************************************
!! ARQUIVOS DE SORTEIOS
!
! Objetivos:
!   Rotinas para manipular arquivos de sorteios. Especificamente, a
!   leitura dos arquivos de pre-valores para condicionamento.
!   
! Criado:
!   01 de maio de 2025
!
! Modificado:
!   18 de janeiro de 2026
! 
! Autoria:
!   oap
!
MODULE arquivos_json
    USE tipos
    USE json_utils_mod

    IMPLICIT NONE
    CHARACTER(5)  :: DIR_OUT  = "./out"
    CHARACTER(18) :: DIR_AVI  = "/valores_iniciais/"

CONTAINS

! ************************************************************
!! Leitura de valores de sorteio
!
! Modificado:
!   18 de janeiro de 2026
!
! Autoria:
!   oap
! 
SUBROUTINE ler_json (arquivo, dados, p_exibir)
    CHARACTER(LEN=*), INTENT(IN) :: arquivo
    TYPE(json_value), POINTER, INTENT(INOUT) :: dados
    LOGICAL, OPTIONAL :: p_exibir
    LOGICAL :: exibir

    exibir = MERGE(p_exibir, .TRUE., PRESENT(p_exibir))

    IF (exibir) WRITE (*,'(A)') 'PRESET: '//arquivo
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
!   10 de novembro de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE salvar_auto_vi (out_dir, infos_sorteio, massas, posicoes, momentos)
    TYPE(json_value), POINTER   :: infos_sorteio
    REAL(pf)                    :: massas(:), posicoes(:,:), momentos(:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: out_dir

    IF (PRESENT(out_dir)) THEN
        CALL salvar_vi_json(TRIM(out_dir)//DIR_AVI, &
                            infos_sorteio, massas, posicoes, momentos, .TRUE.)
    ELSE
        CALL salvar_vi_json(DIR_OUT//DIR_AVI, &
                            infos_sorteio, massas, posicoes, momentos, .TRUE.)
    ENDIF
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
    TYPE(json_core)           :: json
    TYPE(json_value), POINTER :: modo, vi, ap, am
    LOGICAL, OPTIONAL :: gerar_nome
    INTEGER :: i
    LOGICAL :: vi_existe
    CHARACTER(17) :: nome_arq

    ! Verificando nome para arquivo de saida
    IF (PRESENT(gerar_nome) .AND. gerar_nome) THEN
        nome_arq = gerar_nome_arquivo(diretorio, "json")
        WRITE (*,*)
        WRITE(*,*) ' > nome do arquivo: ' // nome_arq
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
        CALL json % add(vi, 'massas', massas)
        
        ! Adiciona as posicoes e os momentos
        CALL json % create_array(ap, 'posicoes')
        CALL json % create_array(am, 'momentos')
        CALL json % add(vi, ap)
        CALL json % add(vi, am)
        
        DO i=1, SIZE(massas)
            CALL json % add(ap, 'posicoes', posicoes(i,:))
            CALL json % add(am, 'momentos', momentos(i,:))
        END DO   
    ENDIF

    ! Salva o arquivo
    CALL json % print(infos_sorteio, TRIM(diretorio//nome_arq))
END SUBROUTINE salvar_vi_json

END MODULE arquivos_json