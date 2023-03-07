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
  public arquivo, ler_csv

  ! classe de arquivo
  type :: arquivo

  ! id do arquivo
  integer :: idarq, qntdCorpos_int, dimensao_int
  ! nome do arquivo, qntd de corpos, formato e dimensão
  character(:), allocatable :: nomearq, formato, formatoMassas, qntdCorpos, dimensao
  ! extensão
  character(4) :: extensao = '.csv'
  ! diretório padrão (fora da pasta build)
  character(3) :: dir = "../"
  
  contains
    procedure :: criar, escrever, fechar, nomeArquivo, criarFormato, escrever_massas

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

    self % qntdCorpos_int = qntdCorpos
    self % dimensao_int = dimensao

    self % formato = '(2(' // self % qntdCorpos // '(' // self % dimensao // '(F15.7, :, ","))))'
    self % formatoMassas = '(' // self % qntdCorpos // '(F15.7, :, ","))'


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

  ! escrever massas
  subroutine escrever_massas (self, massas)

      implicit none
      class(arquivo), intent(in) :: self
      real, intent(in)           :: massas(:)

      ! salva
      write (self % idarq, self % formatoMassas) massas

  end subroutine escrever_massas

  ! escrever no arquivo
  subroutine escrever (self, array)

    implicit none
    class(arquivo), intent(in) :: self
    real, intent(in)           :: array(2,self % qntdCorpos_int,self % dimensao_int)

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
  subroutine ler_csv (nome, massas, R, P)

    character(len=*), intent(in) :: nome
    real, allocatable, intent(inout) :: R(:,:,:), P(:,:,:), massas(:)
    character(len=10000) :: massas_string
    integer :: iu, i, qntdLinhas = 0, io, qntdCorpos = 0

    open(newunit=iu,file=nome,status='old',action='read')

    ! captura a string de massas
    read(iu,'(A)') massas_string    
    
    ! captura a quantidade de corpos a partir da quantidade de vírgulas
    do i = 1, len(massas_string)
      if (massas_string(i:i) == ',') then
        qntdCorpos = qntdCorpos + 1
      end if
    end do
    qntdCorpos = qntdCorpos + 1 ! número de vírgulas = N - 1

    ! captura o número de linhas do CSV
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

    ! captura as posições e momentos
    do i = 1, qntdLinhas
      read(iu,*) R(i,:,:),P(i,:,:)
      R(i,:,:) = transpose(R(i,:,:))
      P(i,:,:) = transpose(P(i,:,:))
    end do
    
    close(iu)

  end subroutine ler_csv

end module arquivos