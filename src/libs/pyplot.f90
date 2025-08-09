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
!   02 de fevereiro de 2024 (criado)
!   08 de agosto de 2025 (modificado)
! 
! Autoria:
!   oap
! 
MODULE pyplot
  USE tipos
  USE arquivos_mod
CONTAINS    

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
SUBROUTINE salvar_dados_xy (arquivo,dim,x,y)

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
  ENDIF

  OPEN( unit=u, file=arquivo, status='replace', iostat=ios )

  IF (ios /= 0) THEN
    WRITE(*, '(a)') 'PLOTAR_XY - ERRO FATAL!'
    WRITE(*, '(a)') '  Nao foi possivel abrir o arquivo de saida.'
  ENDIF

  DO i=1, dim
    WRITE(u, *) x(i), y(i)
  END DO

  CLOSE(unit=u)

END SUBROUTINE salvar_dados_xy

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
SUBROUTINE plotar_trajetorias (R, absc, orde)

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
  LOGICAL                   :: existe=.TRUE.

  ! Verifica se o diretorio de plots existe, senao cria
  CALL diretorio_out()
  INQUIRE(file=dir,exist=existe)
  IF (.not. existe) THEN
    CALL criar_dir('plot','out')
  ENDIF

  ! Cria um nome para o diretorio a partir da data
  CALL nome_data(dir, novo_dir)
  ! Cria o diretorio para os plots
  CALL criar_dir(novo_dir,dir)

  ! Captura a escala da quantidade de corpos (10**x)
  escala = FLOOR(LOG10(FLOAT(SIZE(R,2)))) + 1
  ! Aloca o espaco para o nome do arquivo
  ALLOCATE(CHARACTER(41 + escala) :: data_dir)
  do i = 1, SIZE(R,2) ! quantidade de corpos
    WRITE(i_int,'(I4.4)') i      
    data_dir = dir//TRIM(novo_dir) // '/corpo_' // TRIM(i_int) // '.txt'
    WRITE(*,*) data_dir
    CALL salvar_dados_xy(data_dir, SIZE(R,1), R(:,i,absc), R(:,i,orde))
  END DO
  
  ! Plota
  CALL rodar_plot(novo_dir)

END SUBROUTINE plotar_trajetorias

! ************************************************************
!! Cria um nome pela data
!
! Objetivos:
!   Cria o nome de arqiuvo a partir da data.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE nome_data (dir, nome)

  IMPLICIT NONE
  CHARACTER(LEN=*) :: dir       ! diretorio
  CHARACTER(LEN=*) :: nome      ! nome final
  LOGICAL      :: existe=.TRUE. ! se o arquivo existe
  CHARACTER(9) :: data_hoje     ! captura a data
  CHARACTER(3) :: num           ! (1-999)
  INTEGER      :: i = 1         ! iterador

  ! Data em string
  CALL DATE_AND_TIME(data_hoje)
  
  DO WHILE (existe)
    WRITE (num, '(I3.3)') i

    ! Criacao do nome
    nome = TRIM(TRIM(data_hoje) // '_' // TRIM(num))
    
    ! Verifica se existe
    INQUIRE(file=dir//nome,exist=existe)
    
    i = i+1
  END DO

END SUBROUTINE nome_data

! ************************************************************
!! Rodar plot com matplotlib
!
! Objetivos:
!   Roda o script "plot.py" e exibe um grafico com matplotlib
!
! Modificado:
!   02 de fevereiro de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE rodar_plot (diretorio)
  IMPLICIT NONE
  CHARACTER(LEN=*) :: diretorio
  CALL SYSTEM ('python src/python/plot.py '//TRIM(diretorio))
END SUBROUTINE rodar_plot

END MODULE pyplot