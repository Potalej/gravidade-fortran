! ************************************************************
!! GRAVIDADE-FORTRAN
!
! Objetivos:
!   Arquivo central do gravidade-fortran, sendo o que chama
!   todas as outras funcoes.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
PROGRAM main
  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  USE OMP_LIB
  USE simulacao_sorteio
  USE simulacao_vi
  USE leitura
  USE condicoesIniciais
  USE plot
  IMPLICIT NONE

  CHARACTER(256)  :: arq ! Arquivo
  CHARACTER(15)   :: acao ! Modo em que sera operado
  INTEGER :: i

  IF (command_argument_count() == 0) THEN
    STOP 'Nenhum parametro informado! Para mais informacoes, use --help'
  ENDIF

  ! O primeiro argumento deve ser a acao
  CALL get_command_argument(1, acao)

  ! O segundo deve ser o arquivo
  CALL get_command_argument(2, arq)

  CALL cabecalho

  ! Acoes
  SELECT CASE (acao)
    ! Helper
    CASE ('-h', '--help')
      CALL help()
      STOP
    
    ! Usa um preset para gerar valores e simular
    CASE ('-s', '--sorteio')
      CALL simular_sorteio(arq)

    ! Usa um preset para gerar valores e salva-los, sem simular
    CASE ('-sv', '--sorteio-salvar')
      CALL sorteio_salvar(arq)
    
    ! Valores iniciais
    CASE ('-vi', '--valores-iniciais')
      CALL simular_vi(arq)

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