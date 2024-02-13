! Arquivo
! 
! Para o manejo das questoes voltadas para arquivos .CSV, como a
! abertura, a leitura e a escrita destes.
! 
! O objeto `arquivo` tem as propriedades de criacao, escrita, fechamento,
! criacao de nome e criacao de formato.
! 
! = subroutine criar (self, idarq, qntdCorpos, dimensao)
! cria um arquivo .CSV que se configura para fazer a formatacao de uma
! determinada quantidade de particulas e uma determinada quantidade de 
! dimensoes. O `idarq` é utilizado como identificar unico do arquivo.
! 
! = subroutine escrever (self, array)
! dado um real :: array(2,3,3), cuja ideia eh ser (/Rk, Pk/), os vetores
! sao salvos no CSV.
! 
! = subroutine fechar (self)
! fecha o arquivo quando termina a simulacao.
! 
! = subroutine criarFormato (self, qntdCorpos, dimensao)
! para uma determinada quantidade de particulas e dimensoes, eh criada a
! formatacao para transformar o array em uma string corretamente. Sendo
! N := qntdCorpos e D := dimensao, por exemplo, deve ser gerada a seguinte
! formatacao: '(2 (N ( D( F15.7, :, "," ) )))'
! 
! = subroutine nomeArquivo (self)
! cria o nome do arquivo baseado no dia corrente e contando a partir do 10,
! ou seja, nao havendo arquivos CSV do mesmo dia eh criado um AAAAMMDD_10.csv,
! e caso haja eh criado AAAAMMDD_11.csv, AAAAMMDD_12.csv, etc.
! 

module arquivos
  use, intrinsic :: iso_fortran_env, only: pf=>real64

  ! OPENMP
  use OMP_LIB

  implicit none
  private
  public arquivo, ler_csv, criar_dir, salvar_sorteio, espacosVazios, capturar_unidade, diretorio_out
  public arquivo, ler_csv, criar_dir, salvar_sorteio, espacosVazios, capturar_unidade, diretorio_out

  ! classe de arquivo
  type :: arquivo

  ! id do arquivo
  integer :: idarq, qntdCorpos_int, dimensao_int
  ! nome do arquivo, qntd de corpos, formato e dimensao
  character(:), allocatable :: nomearq, formato, formatoMassas, qntdCorpos, dimensao
  ! extensão
  character(4) :: extensao = '.csv'
  ! diretório padrão (fora da pasta build)
  character(11) :: dir = "./out/data/"
  
  contains
    procedure :: criar, escrever, fechar, nomeArquivo, criarFormato, escrever_massas, escrever_cabecalho, &
                 diretorio_data
    procedure :: criar, escrever, fechar, nomeArquivo, criarFormato, escrever_massas, escrever_cabecalho, &
                 diretorio_data
  end type

contains

  subroutine diretorio_out ()
    implicit none
    logical :: existe = .true.
    ! verifica se existe o diretorio out
    inquire(file="./out", exist=existe)
    if (.NOT. existe) then
      call criar_dir("./out")
    endif
  end subroutine diretorio_out

  subroutine diretorio_data (self)
    implicit none
    class(arquivo), intent(inout) :: self
    logical :: existe = .true.
    call diretorio_out()
    ! verifica se existe o diretorio padrao
    inquire(file=trim(self % dir), exist=existe)
    if (.NOT. existe) then
      call criar_dir(trim(self % dir))
    endif
  end subroutine diretorio_data

  ! para criacao do nome do arquivo
  subroutine nomeArquivo (self)

    implicit none
    class(arquivo), intent(inout) :: self

    ! para iterar e nao repetir arquivo
    integer :: i = 1
    character(3) :: numero
    logical :: existe

    ! para capturar a data
    character(8) :: datahoje

    ! verifica se existe o diretorio padrao
    call self % diretorio_data()

    ! Por padrao, existe
    existe = .true.
 
    ! em string
    call date_and_time(datahoje)

    do while (existe)
      write(numero, '(I3.3)') i
      i = i + 1

      ! cria nome 
      self % nomearq = self % dir//trim(datahoje)//"_"//trim(numero)//self % extensao

      ! verifica se existe
      inquire(file=trim(self % nomearq), exist=existe)
    end do

  end subroutine nomeArquivo

  ! formatacao do arquivo
  subroutine criarFormato (self, qntdCorpos, dimensao)

    implicit none
    class(arquivo), intent(inout) :: self
    integer, intent(in)           :: qntdCorpos, dimensao
    character                     :: formato
    integer                       :: i = 1

    ! salva a quantidade de corpos e dimensao
    self % qntdCorpos = espacosVazios(qntdCorpos)
    self % dimensao = espacosVazios(dimensao)

    self % qntdCorpos_int = qntdCorpos
    self % dimensao_int = dimensao

    self % formato = '(2(' // self % dimensao // '(' // self % qntdCorpos // '(F25.13, :, ","))))'
    self % formatoMassas = '(' // self % qntdCorpos // '(F25.7, :, ","))'


  end subroutine criarFormato

  ! criacao do arquivo
  subroutine criar (self, idarq, qntdCorpos, dimensao)

    implicit none
    class(arquivo), intent(inout) :: self
    integer, intent(in)           :: idarq, qntdCorpos, dimensao
    character(len=18)             :: nomearq
    logical                       :: existe ! para verificar se ja existe ou nao

    WRITE (*, '(a)') 'CRIAR ARQUIVO PARA SALVAR PLOT:'

    ! cria formatacao
    call self % criarFormato(qntdCorpos, dimensao)
    WRITE (*, '(a)') '  > formato : ' // self % formato

    ! criacao do nome do arquivo
    call self % nomeArquivo()
    WRITE (*,'(a)') '  > arquivo de saída: ' // self % nomearq

    ! agora cria o arquivo
    self % idarq = idarq
    open(idarq, file = self % nomearq, status='new')
    WRITE (*,'(a)') '  > arquivo criado!'

    WRITE(*,*)

  end subroutine criar

  ! escrever massas
  subroutine escrever_massas (self, massas)

    implicit none
    class(arquivo), intent(in) :: self
    real(pf), intent(in)       :: massas(:)

    ! salva
    write (self % idarq, self % formatoMassas) massas

  end subroutine escrever_massas

  ! Escreve o cabecalho (contendo h, G, massas)
  subroutine escrever_cabecalho (self, h, G, massas)

    implicit none
    class(arquivo), intent(in) :: self
    real(pf), intent(in)       :: massas(:)
    real(pf)                   :: h, G

    ! Salva h
    write (self % idarq, "(F25.7, :, ',')") h

    ! Salva G
    write (self % idarq, "(F25.7, :, ',')") G

    ! Salva as massas
    call self%escrever_massas(massas)

  end subroutine escrever_cabecalho

  ! escrever no arquivo
  subroutine escrever (self, array)

    implicit none
    class(arquivo), intent(in) :: self
    real(pf), intent(in)       :: array(2,self % qntdCorpos_int,self % dimensao_int)

    ! salva 
    write (self % idarq, self % formato) array

  end subroutine escrever

  ! fechar arquivo
  subroutine fechar (self)

    implicit none
    class(arquivo), intent(in) :: self
    
    close(self % idarq)

  end subroutine fechar

  ! ler arquivo CSV
  subroutine ler_csv (nome, h, G, massas, R, P)

    character(len=*), intent(in)         :: nome
    real(pf), intent(inout)              :: G, h
    real(pf), allocatable, intent(inout) :: R(:,:,:), P(:,:,:), massas(:)
    character(len=10000)                 :: massas_string
    integer :: iu, i, qntdLinhas = 0, io, qntdCorpos = 0
    real(pf) :: t0, tf

    WRITE(*, '(a)') "LER_CSV:"
    WRITE(*, '(a)') "  > arquivo: " // trim(nome)

    open(newunit=iu,file=nome,status='old',action='read')

    ! Captura o tamanho do passo
    read(iu, *) h

    ! Captura a gravidade
    read(iu, *) G

    ! captura a string de massas
    read(iu,'(A)') massas_string    
    
    ! captura a quantidade de corpos a partir da quantidade de virgulas
    do i = 1, len(massas_string)
      if (massas_string(i:i) == ',') then
        qntdCorpos = qntdCorpos + 1
      end if
    end do
    qntdCorpos = qntdCorpos + 1 ! numero de virgulas = N - 1

    ! captura o numero de linhas do CSV
    do while (.true.) 
      read(iu,*, iostat=io)
      if (io /= 0) exit
      qntdLinhas = qntdLinhas+1
    end do
    
    rewind(iu)
    ! Captura o tamanho do passo
    read(iu, *) h

    ! Captura a gravidade
    read(iu, *) G    
    
    ! captura as massas
    allocate(massas(qntdCorpos))
    read(iu, *) massas

    ! aloca os tamanhos
    allocate(R(qntdLinhas,qntdCorpos,3))
    allocate(P(qntdLinhas,qntdCorpos,3))

    ! captura as posicoes e momentos
    t0 = omp_get_wtime()
    do i = 1, qntdLinhas-1
      read(iu,*) R(i,:,:),P(i,:,:)
    end do
    tf = omp_get_wtime()
    
    close(iu)

    WRITE (*,'(a,F10.4,a)') "  > tempo de leitura: ", tf-t0, "s"
    WRITE (*,*)

  end subroutine ler_csv

  ! Criacao de diretorio
  subroutine criar_dir (dir, onde)

    IMPLICIT NONE
    CHARACTER(LEN=*) :: dir
    CHARACTER(LEN=*),OPTIONAL :: onde
    CHARACTER(LEN=len(dir)) :: res
    CHARACTER(:), ALLOCATABLE :: comando
    INTEGER :: i
    res = dir
    ! Remove o "./" se tiver
    do i = 1, len(dir)
      if (dir(i:i) == "/" .OR. dir(i:i) == ".") then
        res(i:i+1) = " "
      end if
    end do

    ALLOCATE(CHARACTER(3+LEN(onde)+10+LEN(res)) :: comando)
    comando = "cd "//onde//" && mkdir "// trim(res)

    call SYSTEM(comando)

  end subroutine criar_dir

  ! Escreve um preset sorteado como um preset de valores iniciais
  subroutine salvar_sorteio (onde, subdir, arquivo, nome, G, massas, R, P, t0, tf, timestep, metodo, corretor, colisoes)

    IMPLICIT NONE
    CHARACTER(LEN=*)      :: onde, subdir, arquivo, metodo, nome
    CHARACTER(LEN=256)    :: dir_arquivo 
    CHARACTER(LEN=3)      :: num_arquivo
    LOGICAL               :: corretor, colisoes, diretorio_existe, arquivo_existe
    REAL(pf)              :: G, t0, tf, timestep
    REAL(pf),allocatable  :: massas(:), R(:,:), P(:,:)
    INTEGER               :: u, i, arq_i

    WRITE(*,'(a)') 'SALVAR SORTEIO:'

    ! Verifica se o diretorio desejado existe
    inquire(file=onde//trim(subdir), exist=diretorio_existe)
    if (.NOT. diretorio_existe) then
      call criar_dir (subdir, onde)
    end if

    ! Agora verifica se o arquivo ja existe
    dir_arquivo = TRIM(onde//subdir) // TRIM(arquivo)
    inquire(file=TRIM(dir_arquivo), exist=arquivo_existe)
    if (arquivo_existe) then
      arq_i = 1
      DO WHILE (arquivo_existe)
        WRITE(num_arquivo, '(I3.3)') arq_i
        inquire(file=TRIM(dir_arquivo)//"_"//TRIM(num_arquivo)//".txt", exist=arquivo_existe)
      END DO
      dir_arquivo = TRIM(dir_arquivo)//"_"//TRIM(num_arquivo)//".txt"
    endif

    WRITE(*,'(a)') '  > arquivo: ' // trim(dir_arquivo)

    ! Abre um arquivo
    OPEN(newunit=u,file=dir_arquivo)

    WRITE(u,'(*(g0,1x))') "! Configs"
    WRITE(u,'(*(g0,1x))') "modo vi"
    WRITE(u,'(*(g0,1x))') "nome ", nome
    WRITE(u,'(*(g0,1x))') "integrador ", metodo
    WRITE(u,'(*(g0,1x))') "timestep ", timestep
    WRITE(u,'(*(g0,1x))') "passos ", floor((tf-t0)/timestep)
    WRITE(u,'(*(g0,1x))') "t0 ", t0
    WRITE(u,'(*(g0,1x))') "tf ", tf
    WRITE(u,'(*(g0,1x))') "corretor ", corretor
    WRITE(u,'(*(g0,1x))') "colisoes ", colisoes

    WRITE(u,*)

    WRITE(u,'(*(g0,1x))') "! Valores do problema"
    WRITE(u,'(*(g0,1x))') "N ", size(massas)
    WRITE(u,'(*(g0,1x))') "G ", G
    
    WRITE(u,*) 

    WRITE(u,'(*(g0,1x))') "! Massas"
    do i = 1, size(massas)
      WRITE(u,'(*(g0,1x))') massas(i)
    end do

    WRITE(u,*)

    WRITE(u,'(*(g0,1x))') "! Posicoes"
    do i = 1, size(massas)
      WRITE(u,'(*(g0,1x,","))') R(i,:) 
    end do

    WRITE(u,*)

    WRITE(u,'(*(g0,1x))') "! Momentos"
    do i = 1, size(massas)
      WRITE(u,'(*(g0,1x,","))') P(i,:) 
    end do

    CLOSE(u)

    WRITE(*,'(a)') '  > arquivo salvo!'
    WRITE(*,*)

  end subroutine salvar_sorteio

   ! para remover espacos vazios
  function espacosVazios (valor)

    implicit none
    integer, intent(in)           :: valor
    character(7)                  :: valor_str
    character(:), allocatable     :: valor_str_parcial, espacosVazios
    integer                       :: i = 1

    ! transforma o valor em string
    write(valor_str, '(I7)') valor

    ! alinha a esquerda para facilitar
    valor_str = adjustl(valor_str)

    ! onde ficara salvo
    valor_str_parcial = ""
    
    ! elimina os caracteres vazios
    do while (.true.)
      if (valor_str(i:i).eq." ") then
        i = 1
        exit
      else
        valor_str_parcial = valor_str_parcial // valor_str(i:i)
        i = i + 1
      end if     
    end do

    ! aloca a string para poder salvar
    allocate( character(len_trim(valor_str_parcial)) :: espacosVazios)
    ! enfim, salva
    espacosVazios = trim(valor_str_parcial)

  end function espacosVazios

  !*****************************************************************************
  !
  !! Retorna uma unidade FORTRAN que esteja livre
  !
  !  Objetivos:
  !
  !    Uma unidade de FORTRAN "livre" eh um inteiro entre 1 e 99 que nao esta
  !    associado a nenhum dispositivo I/O, e eh utilizado para abrir arquivos.
  !    Se a unidade eh nula, entao nao ha nenhuma unidade FORTRAN livre.
  !    
  !    Os numeros 5, 6 e 9 sao reservados, entao nunca sao retornados.
  !    
  !    O codigo foi baseado na biblioteca GNUFOR de John Burkardt.
  !
  !  Modificado:
  !
  !    02 de fevereiro de 2024
  !
  !  Autoria:
  !
  !    oap
  !
  !  Parametros:
  !
  !    Output, integer ( kind = 4 ) IUNIT, o numero de unidade livre.
  !
  subroutine capturar_unidade ( iunit )
    IMPLICIT NONE
    INTEGER ( kind = 4 ) i
    INTEGER ( kind = 4 ) ios
    INTEGER ( kind = 4 ) iunit
    LOGICAL lopen
    iunit = 0
    DO i = 1, 99
      IF (i /= 5 .and. i /= 6 .and. i /= 9) THEN
        INQUIRE ( unit = i, opened = lopen, iostat = ios )
        IF ( ios == 0 ) THEN
          IF ( .not. lopen ) THEN
            iunit = i
            RETURN
          END IF
        END IF
      END IF
    END DO
    RETURN
  end
end module arquivos