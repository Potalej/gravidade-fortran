MODULE string_utils
  IMPLICIT NONE
CONTAINS
FUNCTION int_para_string (inteiro) RESULT(string)
  INTEGER, INTENT(IN)           :: inteiro
  CHARACTER(32)                 :: buffer
  CHARACTER(:), ALLOCATABLE     :: string_parcial, string
  INTEGER                       :: i = 1

  ! transforma o valor em string
  WRITE(buffer, '(I32)') inteiro

  ! alinha a esquerda para facilitar
  buffer = ADJUSTL(buffer)

  ! onde ficara salvo
  string_parcial = ""
  
  ! elimina os caracteres vazios
  DO WHILE (.TRUE.)
    IF (buffer(i:i).eq." ") THEN
      i = 1
      exit
    ELSE
      string_parcial = string_parcial // buffer(i:i)
      i = i + 1
    ENDIF     
  END DO

  ! aloca a string para poder salvar
  ALLOCATE(CHARACTER(LEN_TRIM(string_parcial)) :: string)
  ! enfim, salva
  string = TRIM(string_parcial)
END FUNCTION int_para_string

! ************************************************************
!! Captura a data e a hora no formato 11:58 24/07/2025
!
! Modificado:
!   24 de julho de 2025
!
! Autoria:
!   oap
!
SUBROUTINE data_hora_string (data_hora_str)
  CHARACTER(LEN=20), INTENT(OUT) :: data_hora_str
  INTEGER :: v(8)

  ! v = [ano, mês, dia, fuso, hora, min, seg, milisseg]
  call date_and_time(values=v)

  WRITE(data_hora_str, '(I2.2,":",I2.2," ",I2.2,"/",I2.2,"/",I4)') v(5), v(6), v(3), v(2), v(1)
END SUBROUTINE data_hora_string

END MODULE string_utils