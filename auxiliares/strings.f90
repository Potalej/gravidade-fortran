! strings
! 
! métodos para facilitar o trabalho com strings
! 
! = function espacosVazios (valor)
! dado um valor integer(7), cria uma string sem espaços
! vazios com esse inteiro.
! 

module strings

  implicit none
  private
  public espacosVazios
  
contains
  ! para remover espaços vazios
  function espacosVazios (valor)

    implicit none
    integer, intent(in)           :: valor
    character(7)                  :: valor_str
    character(:), allocatable     :: valor_str_parcial, espacosVazios
    integer                       :: i = 1

    ! transforma o valor em string
    write(valor_str, '(I7)') valor

    ! alinha à esquerda para facilitar
    valor_str = adjustl(valor_str)

    ! onde ficará salvo
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

end module strings