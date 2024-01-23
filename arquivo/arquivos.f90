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
  ! pacote de strings
  use strings

  implicit none
  private
  public arquivo, ler_csv, criar_dir

  ! classe de arquivo
  type :: arquivo

  ! id do arquivo
  integer :: idarq, qntdCorpos_int, dimensao_int
  ! nome do arquivo, qntd de corpos, formato e dimensao
  character(:), allocatable :: nomearq, formato, formatoMassas, qntdCorpos, dimensao
  ! extensão
  character(4) :: extensao = '.csv'
  ! diretório padrão (fora da pasta build)
  character(8) :: dir = "../data/"
  
  contains
    procedure :: criar, escrever, fechar, nomeArquivo, criarFormato, escrever_massas, escrever_cabecalho

  end type

contains

  ! para criacao do nome do arquivo
  subroutine nomeArquivo (self)

    implicit none
    class(arquivo), intent(inout) :: self

    ! para iterar e nao repetir arquivo
    integer :: i = 1
    character(3) :: numero
    logical :: existe

    ! para capturar a data
    character(9) :: datahoje

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

    ! cria formatacao
    call self % criarFormato(qntdCorpos, dimensao)
    WRITE (*,*) 'Formato : ', self % formato

    ! criacao do nome do arquivo
    call self % nomeArquivo()
    WRITE (*,*) 'Arquivo de saída: ', self % nomearq

    ! agora cria o arquivo
    self % idarq = idarq
    open(idarq, file = self % nomearq, status='new')

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
    
    ! captura as massas
    rewind(iu)
    allocate(massas(qntdCorpos))
    read(iu, *) massas

    ! aloca os tamanhos
    allocate(R(qntdLinhas,qntdCorpos,3))
    allocate(P(qntdLinhas,qntdCorpos,3))

    ! captura as posicoes e momentos
    do i = 1, qntdLinhas-1
      read(iu,*) R(i,:,:),P(i,:,:)
      ! R(i,:,:) = transpose(R(i,:,:))
      ! P(i,:,:) = transpose(P(i,:,:))
    end do
    
    close(iu)

  end subroutine ler_csv

  ! Criacao de diretorio
  subroutine criar_dir (dir)

    IMPLICIT NONE
    CHARACTER(LEN=*) :: dir

    call SYSTEM("mkdir " // trim(dir))

  end subroutine criar_dir
end module arquivos