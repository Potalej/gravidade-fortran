program main
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use OMP_LIB
  use simulacao_sorteio
  use simulacao_vi
  use leitura
  use condicoesArtigo
  use plot
  implicit none

  character(256)  :: preset  ! Configuracao de preset
  character(256)  :: vi      ! Valores iniciais
  character(1024) :: comando ! Comando todo
  character(256)  :: exibir  ! Exibir trajetorias

  type(preset_config) :: preset_configs

  call ler_cli      ! Le os argumentos
  call formatar_cli ! Formata os argumentos
  call rodar        ! Roda o programa
contains

  subroutine ler_cli
    integer :: comprimento
    integer :: io, io2
    character(200) :: arg
    integer :: i

    if (command_argument_count() == 0) then
      ERROR STOP 'Nenhum parametro informado! Para mais informacoes, use --help'
    endif

    ! Tratamento de casos de parametros
    do i=1, command_argument_count()
      call get_command_argument(i, arg)

      select case (arg)
        case ('-h', '--help')
          call help()
          stop
      end select
    end do

    ! Salva os comandos
    comando = ""
    call get_command(command = comando, status=io)
    if (io==0) then
      call get_command_argument(0, length=comprimento, status=io2)  
      if (io2==0) then
        comando = "&cmd "//adjustl(trim(comando(comprimento+1:)))//" /"
      else
        comando = "&cmd "//adjustl(trim(comando))//" /"
      end if
    else
      write(*,*) io,"Erro capturando a linha de comando."
    end if
  end subroutine

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
     
  subroutine formatar_cli
    character(256) :: msg
    namelist /cmd/ preset
    namelist /cmd/ vi
    namelist /cmd/ exibir
    integer :: io

    if (len_trim(comando)>0) then
      msg = ''
      read(comando,nml = cmd,iostat = io,iomsg = msg)
      if (io/=0) then
        error stop "Erro analisando a linha de comando ou cmd.conf " // msg
      end if
    end if
  end subroutine formatar_cli

  subroutine rodar
    IMPLICIT NONE
    namelist /cmd/ preset
    namelist /cmd/ vi
    namelist /cmd/ exibir
    LOGICAL :: sorteio = .false.          ! PROVISORIO
    LOGICAL :: valores_iniciais = .false. ! PROVISORIO
    LOGICAL :: exibir_grafico = .true.    ! PROVISORIO

    if (preset == "" .AND. vi == "") then
      WRITE (*,*) 'Entrada vazia!'

    else if (sorteio) then
      call cabecalho
      call simular_sorteio(preset)

    else if (valores_iniciais) then
      call cabecalho
      call simular_vi(vi)

    else if (exibir_grafico) then
      call cabecalho
      call rodar_exibir(exibir)
    endif

  end subroutine rodar

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