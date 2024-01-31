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

    ! Default
    CASE default
      WRITE (*,*) 'Entrada vazia!!'
      STOP
  end select

contains

  ! Ajudador
  subroutine help
    WRITE (*,*) 'SINOPSE:'
    WRITE (*,*) '   # ./gravidade [--help|-h]'
    WRITE (*,*) '   # ./gravidade [OPCAO] [ARQUIVO]'
    WRITE (*,*)
    WRITE (*,*) 'OPCOES:'
    WRITE (*,*) '   -h, --help'
    WRITE (*,*) '         Exibe a ajuda.'
    WRITE (*,*)
    WRITE (*,*) '   -s, --sorteio' 
    WRITE (*,*) '         Utiliza um preset informado para gerar valores iniciais aleatorios'
    WRITE (*,*)
    WRITE (*,*) '   -vi, --valores-iniciais' 
    WRITE (*,*) '         Utiliza os valores iniciais contidos em um arquivo informado para simular.'
    WRITE (*,*)
    WRITE (*,*) '   -e, --exibir' 
    WRITE (*,*) '         Gera um grafico com as trajetorias dos corpos no arquivo informado.'
    WRITE (*,*)
  end subroutine help

  ! Imprime o cabecalho do programa
  subroutine cabecalho
    write(*,*)
    write(*,*) "|=========================================|"
    write(*,*) "|                    _    _         _     |"
    write(*,*) "|  __ _ _ _ __ ___ _(_)__| |__ _ __| |___ |"
    write(*,*) "| / _` | '_/ _` \ V / / _` / _` / _` / -_)|"
    write(*,*) "| \__, |_| \__,_|\_/|_\__,_\__,_\__,_\___||"
    write(*,*) "| |___/                                   |"
    write(*,*) "|=========================================|"
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