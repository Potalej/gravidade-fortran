program main
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use OMP_LIB
  use simulacao_sorteio
  use simulacao_vi
  use leitura
  use condicoesArtigo
  use plot
  IMPLICIT NONE

  CHARACTER(256)  :: arq ! Arquivo
  CHARACTER(15)   :: acao ! Modo em que sera operado
  INTEGER :: i

  if (command_argument_count() == 0) then
    STOP 'Nenhum parametro informado! Para mais informacoes, use --help'
  endif

  ! O primeiro argumento deve ser a acao
  CALL get_command_argument(1, acao)

  ! O segundo deve ser o arquivo
  CALL get_command_argument(2, arq)

  ! Acoes
  SELECT CASE (acao)
    ! Helper
    CASE ('-h', '--help')
      CALL help()
      STOP
    
    ! Usa um preset para gerar valores
    CASE ('-s', '--sorteio')
      CALL cabecalho
      CALL simular_sorteio(arq)
    
    ! Valores iniciais
    CASE ('-vi', '--valores-iniciais')
      CALL cabecalho
      CALL simular_vi(arq)

    ! Para visualizar as trajetorias
    CASE ('-e', '--exibir')
      CALL cabecalho
      CALL rodar_exibir(arq)

    ! Default
    CASE default
      WRITE (*,*) 'Entrada vazia!!'
      STOP
  end select

contains

  ! Ajudador
  subroutine help
    WRITE (*,*)
    WRITE (*,*) 'Uso:'
    WRITE (*,*) '   # ./gravidade [--help|-h]'
    WRITE (*,*) '   # ./gravidade preset=\"preset.txt\"'
    WRITE (*,*)
    WRITE (*,*) 'Opcoes:'
    WRITE (*,*) '   -h --help    Exibe a ajuda.'
    WRITE (*,*) '   preset       Define um preset para gerar condicoes iniciais' 
    WRITE (*,*)
  end subroutine help

  ! Imprime o cabecalho do programa
  subroutine cabecalho
    write(*,*)
    write(*,*) "|===================================|"
    write(*,*) "| █▀▀ █▀█ ▄▀█ █░█ █ █▀▄ ▄▀█ █▀▄ █▀▀ |"
    write(*,*) "| █▄█ █▀▄ █▀█ ▀▄▀ █ █▄▀ █▀█ █▄▀ ██▄ |"
    write(*,*) "|===================================|"
    write(*,*)
  end subroutine cabecalho
     
  subroutine rodar_exibir (dir)
    IMPLICIT NONE
    CHARACTER(len=*) :: dir
    real(pf), allocatable:: R(:,:,:), P(:,:,:), massas(:)
    real(pf) :: G, h

    ! Abre o arquivo salvo para leitura
    call ler_csv(dir, h, G, massas, R, P)
    
    ! Exibe as trajetorias
    call plotar_trajetorias(R,1,2)

  end subroutine rodar_exibir

end program