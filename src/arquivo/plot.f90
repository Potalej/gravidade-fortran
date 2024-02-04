! *****************************************************************
!! PLOT 
!
! Objetivos:
!   Este arquivo utiliza da biblioteca MATPLOTLIB do Python para 
!   plotar alguns graficos simples, tendo em vista confirmar 
!   visualmente os calculos e os scripts antes de tentar rodar 
!   simulacoes mais pesadas.
! 
! Modificado:
!   02 de fevereiro de 2024
! 
! Autoria:
!   oap
! 
module plot
  use, intrinsic :: iso_fortran_env, only: pf=>real64

  use arquivos
contains    

  ! ************************************************************
  !! salvar_dados_xy
  !
  ! Objetivos:
  !   Salva valores X e Y em coluna em um arquivo .TXT
  !
  ! Modificado:
  !   02 de fevereiro de 2024
  !
  ! Autoria:
  !   oap
  ! 
  subroutine salvar_dados_xy (arquivo,dim,x,y)

    IMPLICIT NONE
    CHARACTER (LEN=*) :: arquivo
    INTEGER (kind=4) u
    INTEGER (kind=4) ios
    INTEGER (kind=4) i
    INTEGER  :: dim
    REAL(pf) :: x(dim), y(dim)

    CALL capturar_unidade(u)

    IF (u == 0) THEN
      WRITE(*, '(a)') 'PLOTAR_XY - ERRO FATAL!'
      WRITE(*, '(a)') '  Nao foi possivel encontrar uma unidade de FORTRAN vazia.'
    END IF

    OPEN( unit=u, file=arquivo, status='replace', iostat=ios )

    IF (ios /= 0) THEN
      WRITE(*, '(a)') 'PLOTAR_XY - ERRO FATAL!'
      WRITE(*, '(a)') '  Nao foi possivel abrir o arquivo de saida.'
    END IF

    DO i=1, dim
      WRITE(u, *) x(i), y(i)
    END DO

    CLOSE(unit=u)

  end subroutine salvar_dados_xy

  ! ************************************************************
  !! plotar_trajetorias
  !
  ! Objetivos:
  !   Plotar em 2d as trajetorias de um problema de N-corpos. Os
  !   eixos utilizados para o plot devem ser informados em absc
  !   (abscissas) e orde (ordenadas), que por padrao sao 1 e 2,
  !   respectivamente.
  !
  ! Modificado:
  !   22 de janeiro de 2023
  !
  ! Autoria:
  !   oap
  ! 
  subroutine plotar_trajetorias (R, absc, orde)

    IMPLICIT NONE
    INTEGER                   :: absc, orde
    REAL(pf),INTENT(INOUT)    :: R(:,:,:)
    INTEGER                   :: i, escala
    INTEGER (kind=4) ierror
    CHARACTER(11)              :: dir = './out/plot/'
    CHARACTER(19)             :: novo_dir
    CHARACTER(28)             :: path_comando
    CHARACTER(4)              :: i_int
    CHARACTER(:), ALLOCATABLE :: data_dir
    LOGICAL                   :: existe=.true.

    ! Verifica se o diretorio de plots existe, senao cria
    call diretorio_out()
    inquire(file=dir,exist=existe)
    if (.not. existe) then
      call criar_dir('plot','out')
    end if

    ! Cria um nome para o diretorio a partir da data
    call nome_data(dir, novo_dir)
    ! Cria o diretorio para os plots
    call criar_dir(novo_dir,dir)

    ! Captura a escala da quantidade de corpos (10**x)
    escala = FLOOR(LOG10(FLOAT(size(R,2)))) + 1
    ! Aloca o espaco para o nome do arquivo
    ALLOCATE(CHARACTER(41 + escala) :: data_dir)
    do i = 1, size(R,2) ! quantidade de corpos
      WRITE(i_int,'(I4.4)') i      
      data_dir = dir//TRIM(novo_dir) // '/corpo_' // TRIM(i_int) // '.txt'
      WRITE(*,*) data_dir
      CALL salvar_dados_xy(data_dir, size(R,1), R(:,i,absc), R(:,i,orde))
    end do
    
    ! Plota
    CALL rodar_plot(novo_dir)

  end subroutine plotar_trajetorias

  subroutine nome_data (dir, nome)

    IMPLICIT NONE
    CHARACTER(LEN=*) :: dir       ! diretorio
    CHARACTER(LEN=*) :: nome      ! nome final
    LOGICAL      :: existe=.true. ! se o arquivo existe
    CHARACTER(9) :: data_hoje     ! captura a data
    CHARACTER(3) :: num           ! (1-999)
    INTEGER      :: i = 1         ! iterador

    ! Data em string
    call date_and_time(data_hoje)
    
    do while (existe)
      WRITE (num, '(I3.3)') i

      ! Criacao do nome
      nome = TRIM(trim(data_hoje) // '_' // TRIM(num))
      
      ! Verifica se existe
      inquire(file=dir//nome,exist=existe)
      
      i = i+1
    end do

  end subroutine nome_data

end module plot