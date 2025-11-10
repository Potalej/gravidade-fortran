! ************************************************************
!! GRAVIDADE-FORTRAN
!
! Objetivos:
!   Arquivo central do gravidade-fortran, sendo o que chama
!   todas as outras funcoes.
!
! Modificado:
!   10 de novembro de 2025
!
! Autoria:
!   oap
! 
PROGRAM main
  USE version
  USE tipos
  USE integradores, only: listar_integradores_disponiveis
  USE simulacao_sorteio
  USE simulacao_vi
  USE pyplot
  IMPLICIT NONE

  CHARACTER(256)  :: arq ! Arquivo
  CHARACTER(15)   :: acao ! Modo em que sera operado
  CHARACTER(15)   :: out_comando
  CHARACTER(256)  :: out_dir

  IF (command_argument_count() == 0) THEN
    STOP 'Nenhum parametro informado! Para mais informacoes, use --help'
  ENDIF

  ! O primeiro argumento deve ser a acao
  CALL get_command_argument(1, acao)

  ! O segundo deve ser o arquivo
  CALL get_command_argument(2, arq)

  ! Se tiver quatro, o terceiro deve ser o "pasta-saida" e o quarto o valor
  IF (command_argument_count() == 4) THEN
    CALL get_command_argument(3, out_comando)
    IF (TRIM(out_comando) /= "-ps" .AND. TRIM(out_comando) /= "--pasta-saida") THEN
      STOP 'O terceiro parametro foi inserido incorretamente. Para mais informacoes, use --help'
    ENDIF

    CALL get_command_argument(4, out_dir)
    IF (out_dir == "") out_dir = 'out'

  ! Se nao, usa o padrao
  ELSE
    out_dir = 'out'
  ENDIF

  CALL cabecalho

  ! Acoes
  SELECT CASE (acao)
    ! Helper
    CASE ('-h', '--help')
      CALL help()
      STOP
    
    ! Usa um preset para gerar valores e simular
    CASE ('-s', '--sorteio')
      CALL simular_sorteio(arq, TRIM(out_dir))

    ! Usa um preset para gerar valores e salva-los, sem simular
    CASE ('-sv', '--sorteio-salvar')
      CALL sorteio_salvar(arq, TRIM(out_dir))
    
    ! Valores iniciais
    CASE ('-vi', '--valores-iniciais')
      CALL simular_vi(arq, TRIM(out_dir))

    ! Para visualizar as trajetorias
    CASE ('-e', '--exibir')
      CALL rodar_exibir(arq)

    ! Para visualizar informacoes
    CASE ('-d', '--dados')
      CALL rodar_dados(arq)

    ! Default
    CASE default
      WRITE (*,*) 'Entrada vazia!!'
      STOP
  END select

CONTAINS

  ! Ajudador
  SUBROUTINE help
    IF (usar_gpu) THEN
      WRITE (*,*) 'v', version_string, '_', precisao, '_GPU (', build_date, ' ', build_time, ')'
    ELSE
      WRITE (*,*) 'v', version_string, '_', precisao, ' (', build_date, ' ', build_time, ')'
    ENDIF
    WRITE (*,*)
    WRITE (*,*) 'SINOPSE:'
    WRITE (*,*) '   # ./gravidade [--help|-h]'
    WRITE (*,*) '   # ./gravidade [OPCAO] [ARQUIVO]'
    WRITE (*,*)
    WRITE (*,*) 'OPCOES:'
    WRITE (*,*) '   -h, --help'
    WRITE (*,*) '         Exibe a ajuda.'
    WRITE (*,*)
    WRITE (*,*) '   (-s, --sorteio) [arquivo]' 
    WRITE (*,*) '         Utiliza um preset informado para gerar valores iniciais aleatorios'
    WRITE (*,*)
    WRITE (*,*) '   (-vi, --valores-iniciais) [arquivo]' 
    WRITE (*,*) '         Utiliza os valores iniciais contidos em um arquivo informado para simular.'
    WRITE (*,*)
    WRITE (*,*) '   (-e, --exibir) [arquivo]' 
    WRITE (*,*) '         Gera um grafico com as trajetorias dos corpos no arquivo informado.'
    WRITE (*,*)
    WRITE (*,*) 'PARAMETROS:'
    WRITE (*,*) '   (-ps, --pasta-saida) [diretorio]'
    WRITE (*,*) '         Cria o diretorio `diretorio` se nao existir e o utiliza como pasta de saida.'
    WRITE (*,*)
    WRITE (*,*) 'INTEGRADORES DISPONIVEIS:'
    
    
    CALL listar_integradores_disponiveis()
  END SUBROUTINE help

  ! Imprime o cabecalho do programa
  SUBROUTINE cabecalho
    WRITE(*,*)
    WRITE(*,*) "|=========================================|"
    WRITE(*,*) "|                    _    _         _     |"
    WRITE(*,*) "|  __ _ _ _ __ ___ _(_)__| |__ _ __| |___ |"
    WRITE(*,*) "| / _` | '_/ _` \ V / / _` / _` / _` / -_)|"
    WRITE(*,*) "| \__, |_| \__,_|\_/|_\__,_\__,_\__,_\___||"
    WRITE(*,*) "| |___/                                   |"
    WRITE(*,*) "|=========================================|"
    WRITE(*,*)
  END SUBROUTINE cabecalho
     
  SUBROUTINE rodar_exibir (dir)
    IMPLICIT NONE
    CHARACTER(len=*) :: dir
    REAL(pf), allocatable:: R(:,:,:), P(:,:,:), massas(:)
    REAL(pf) :: G, h

    ! Abre o arquivo salvo para leitura
    CALL ler_csv(dir, h, G, massas, R, P)
    
    ! Exibe as trajetorias
    CALL plotar_trajetorias(R,1,2)

  END SUBROUTINE rodar_exibir

  SUBROUTINE rodar_dados (dir)
    IMPLICIT NONE
    CHARACTER(len=*) :: dir
    REAL(pf), allocatable:: R(:,:,:), P(:,:,:), massas(:)
    REAL(pf) :: G, h

    ! Abre o arquivo salvo para leitura
    CALL ler_csv(dir, h, G, massas, R, P)
    
    ! Exibe as trajetorias
    CALL plotar_trajetorias(R,1,2)
  END SUBROUTINE

END program