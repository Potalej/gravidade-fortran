! *****************************************************************
!! PLOT 
!
! Objetivos:
!   Este arquivo utiliza da biblioteca GNUFOR para plotar alguns
!   graficos simples, tendo em vista confirmar visualmente os
!   calculos e os scripts antes de tentar rodar simulacoes mais
!   pesadas.
! 
! Modificado:
!   19 de janeiro de 2023
! 
! Autoria:
!   oap
! 
module plot
  use, intrinsic :: iso_fortran_env, only: pf=>real64

  use arquivos
contains    

  ! ************************************************************
  !! plotar_xy
  !
  ! Objetivos:
  !   Plotar graficos 2d, informados a dimensao dos arrays e os
  !   valores dos eixos.
  !
  ! Modificado:
  !   22 de janeiro de 2023
  !
  ! Autoria:
  !   oap
  ! 
  subroutine plotar_xy (dim,x,y)

    implicit none
    INTEGER (kind=4) ierror
    INTEGER  :: dim            
    REAL(pf) :: x(dim), y(dim)

    ! TODO: Adicionar tratamento de erros com o ierror
    call write_xy_data('data_temp.txt',size(x,1),x,y,ierror)
    call write_xy_plot('comando_temp.txt','data_temp.txt',ierror)
    call run_gnuplot('comando_temp.txt')

  end subroutine plotar_xy

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

    implicit none
    INTEGER :: absc, orde
    REAL(pf),INTENT(INOUT)    :: R(:,:,:)
    INTEGER                   :: i, escala
    INTEGER (kind=4) ierror
    CHARACTER                 :: dir = 'data_temp_'
    CHARACTER(19)             :: novo_dir
    CHARACTER(28)             :: path_comando
    CHARACTER(4)              :: i_int
    CHARACTER(:), ALLOCATABLE :: data_dir
    LOGICAL                   :: existe=.true.

    ! Verifica se o diretorio de plots existe, senao cria
    inquire(file='./plot',exist=existe)
    if (.not. existe) then
      call criar_dir('./plot')
    end if

    ! Cria um nome para o diretorio a partir da data
    call nome_data('./plot/', novo_dir)
    print *, novo_dir
    path_comando = novo_dir // '/plot.txt'
    ! Cria o diretorio para os plots
    call criar_dir(novo_dir)

    ! Captura a escala da quantidade de corpos (10**x)
    escala = FLOOR(LOG10(FLOAT(size(R,2)))) + 1
    ! Aloca o espaco para o nome do arquivo
    ALLOCATE(CHARACTER(24 + escala) :: data_dir)
    ! ALLOCATE(CHARACTER(escala) :: i_int)
    print *, "N:", size(R,2)
    do i = 1, size(R,2) ! quantidade de corpos
      write(i_int,'(I4.4)') i      
      data_dir = novo_dir // '/corpo_' // TRIM(i_int) // '.txt'
      call write_xy_data(data_dir,size(R,1),R(:,i,absc), R(:,i,orde),ierror)
    end do
    call escrever_comando_GNU(path_comando, novo_dir, novo_dir)
    call run_gnuplot(path_comando)

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
      nome = TRIM(dir // trim(data_hoje) // '_' // TRIM(num))
      
      ! Verifica se existe
      inquire(file=nome,exist=existe)
      
      i = i+1
    end do

  end subroutine nome_data
  
  subroutine escrever_comando_GNU (nome,titulo,dir)

    IMPLICIT NONE
    CHARACTER(LEN=*) :: nome, titulo, dir
    INTEGER          :: u=13

    OPEN(u,file=nome,status="new")
    WRITE(u,*) 'set title "' // titulo // '"'
    WRITE(u,*) 'set xlabel "x"'
    WRITE(u,*) 'set ylabel "y"'
    WRITE(u,*) 'FILES = system("ls '//dir//'/corpo_*.txt")'
    WRITE(u,*) 'plot for [data in FILES] data u 1:2 w lines'
    WRITE(u,*) 'pause -1'
    WRITE(u,*) 'q'
    CLOSE(U)

  end subroutine escrever_comando_GNU

end module plot