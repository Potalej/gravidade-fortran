! Arquivo
! 
! Para o manejo das questões voltadas para arquivos .CSV, como a
! abertura, a leitura e a escrita destes.
! 
! O objeto `arquivo` tem as propriedades de criação, escrita, fechamento,
! criação de nome e criação de formato.
! 
! = subroutine criar (self, idarq, qntdCorpos, dimensao)
! cria um arquivo .CSV que se configura para fazer a formatação de uma
! determinada quantidade de partículas e uma determinada quantidade de 
! dimensões. O `idarq` é utilizado como identificar único do arquivo.
! 
! = subroutine escrever (self, array)
! dado um real :: array(2,3,3), cuja ideia é ser (/Rk, Pk/), os vetores
! são salvos no CSV.
! 
! = subroutine fechar (self)
! fecha o arquivo quando termina a simulação.
! 
! = subroutine criarFormato (self, qntdCorpos, dimensao)
! para uma determinada quantidade de partículas e dimensões, é criada a
! formatação para transformar o array em uma string corretamente. Sendo
! N := qntdCorpos e D := dimensao, por exemplo, deve ser gerada a seguinte
! formatação: '(2 (N ( D( F15.7, :, "," ) )))'
! 
! = subroutine nomeArquivo (self)
! cria o nome do arquivo baseado no dia corrente e contando a partir do 10,
! ou seja, não havendo arquivos CSV do mesmo dia é criado um AAAAMMDD_10.csv,
! e caso haja é criado AAAAMMDD_11.csv, AAAAMMDD_12.csv, etc.
! 

module arquivos

  ! pacote de strings
  use strings

  implicit none
  private
  public arquivo

  ! classe de arquivo
  type :: arquivo

  ! id do arquivo
  integer :: idarq
  ! nome do arquivo, qntd de corpos, formato e dimensão
  character(:), allocatable :: nomearq, formato, qntdCorpos, dimensao
  ! extensão
  character(4) :: extensao = '.csv'
  ! diretório padrão (fora da pasta build)
  character(3) :: dir = "../"
  
  contains
    procedure :: criar, escrever, fechar, nomeArquivo, criarFormato

  end type

contains

  ! para criação do nome do arquivo
  subroutine nomeArquivo (self)

    implicit none
    class(arquivo), intent(inout) :: self

    ! para iterar e não repetir arquivo
    integer :: i = 10
    character(2) :: numero
    logical :: existe = .true.

    ! para capturar a data
    character(9) :: datahoje
 
    ! em string
    call date_and_time(datahoje)

    do while (existe)
      write(numero, '(I2)') i
      i = i + 1

      ! cria nome 
      self % nomearq = self % dir//trim(datahoje)//"_"//trim(numero)//self % extensao

      ! verifica se existe
      inquire(file=trim(self % nomearq), exist=existe)
    end do

  end subroutine nomeArquivo

  ! formatação do arquivo
  subroutine criarFormato (self, qntdCorpos, dimensao)

    implicit none
    class(arquivo), intent(inout) :: self
    integer, intent(in)           :: qntdCorpos, dimensao
    character                     :: formato
    integer                       :: i = 1

    ! salva a quantidade de corpos e dimensão
    self % qntdCorpos = espacosVazios(qntdCorpos)
    self % dimensao = espacosVazios(dimensao)

    self % formato = '(2(' // self % qntdCorpos // '(' // self % dimensao // '(F15.7, :, ","))))'

  end subroutine criarFormato

  ! criação do arquivo
  subroutine criar (self, idarq, qntdCorpos, dimensao)

    implicit none
    class(arquivo), intent(inout) :: self
    integer, intent(in)           :: idarq, qntdCorpos, dimensao
    character(len=18)             :: nomearq
    logical                       :: existe ! para verificar se já existe ou não

    ! cria formatação
    call self % criarFormato(qntdCorpos, dimensao)
    print *, 'Formato : ', self % formato

    ! criação do nome do arquivo
    call self % nomeArquivo()
    print *, 'Arquivo de saída: ', self % nomearq

    ! agora cria o arquivo
    self % idarq = idarq
    open(idarq, file = self % nomearq, status='new')

  end subroutine criar

  ! escrever no arquivo
  subroutine escrever (self, array)

    implicit none
    class(arquivo), intent(in) :: self
    real, intent(in)           :: array(2,3,3)

    ! salva 
    write (self % idarq, self % formato) array

  end subroutine escrever

  ! fechar arquivo
  subroutine fechar (self)

    implicit none
    class(arquivo), intent(in) :: self
    
    close(self % idarq)

  end subroutine fechar

end module arquivos