MODULE diretorio
  IMPLICIT NONE
  PUBLIC
CONTAINS

! ************************************************************
!! Criacao de diretorio
!
! Objetivos:
!   Cria um diretorio em algum lugar.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
SUBROUTINE criar_dir (dir, onde)

  IMPLICIT NONE
  CHARACTER(LEN=*) :: dir
  CHARACTER(LEN=*),OPTIONAL :: onde
  CHARACTER(LEN=LEN(dir)) :: res
  CHARACTER(:), ALLOCATABLE :: comando
  INTEGER :: i
  res = dir
  ! Remove o "./" se tiver
  DO i = 1, LEN(dir)
    IF (dir(i:i) == "/" .OR. dir(i:i) == ".") THEN
      res(i:i) = " "
    ENDIF
  END DO

  IF (PRESENT(onde)) THEN
    ALLOCATE(CHARACTER(3+LEN(onde)+10+LEN(res)) :: comando)
    comando = "cd "//onde//" && mkdir "// TRIM(res)
  ELSE
    comando = "mkdir "//TRIM(res)
  ENDIF

  CALL SYSTEM(comando)

  DEALLOCATE(comando)

END SUBROUTINE criar_dir

! ************************************************************
!! Verifica se o diretorio de saida existe, e caso nao, cria
!
! Modificado:
!   10 de novembro de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE diretorio_out (dir_param)
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dir_param
  CHARACTER(:), ALLOCATABLE :: dir
  LOGICAL :: existe
  
  IF (PRESENT(dir_param)) THEN
    ALLOCATE(CHARACTER(LEN(TRIM(dir_param))) :: dir)
    dir = dir_param
  ELSE
    ALLOCATE(CHARACTER(5) :: dir)
    dir = "./out"
  ENDIF

  INQUIRE(file=dir, exist=existe)
  IF (.NOT. existe) CALL criar_dir(dir)
END SUBROUTINE diretorio_out

! ************************************************************
!! Verifica se o diretorio de "data" existe, e caso nao, cria
!
! Modificado:
!   10 de novembro de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE diretorio_data (dir_param)
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dir_param
  CHARACTER(:), ALLOCATABLE :: dir
  LOGICAL :: existe
  
  IF (PRESENT(dir_param)) THEN
    ALLOCATE(CHARACTER(LEN(TRIM(dir_param))) :: dir)
    dir = dir_param
  ELSE
    ALLOCATE(CHARACTER(5) :: dir)
    dir = "./out"
  ENDIF
  
  ! Verifica se existe o diretorio de saida
  CALL diretorio_out(dir)

  ! Verifica se existe o subdiretorio "data"
  INQUIRE(file=dir//"/data", exist=existe)
  IF (.NOT. existe) CALL criar_dir('data', dir)
END SUBROUTINE diretorio_data

! ************************************************************
!! Verifica se o diretorio "valores_iniciais" existe, e 
!! caso nao, cria
!
! Modificado:
!   10 de novembro de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE diretorio_vi (dir_param)
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dir_param
  CHARACTER(:), ALLOCATABLE :: dir
  LOGICAL :: existe
  
  IF (PRESENT(dir_param)) THEN
    ALLOCATE(CHARACTER(LEN(TRIM(dir_param))) :: dir)
    dir = dir_param
  ELSE
    ALLOCATE(CHARACTER(5) :: dir)
    dir = "./out"
  ENDIF
  
  ! Verifica se existe o diretorio de saida
  CALL diretorio_out(dir)

  ! Verifica se existe o subdiretorio "data"
  INQUIRE(file=dir//"/valores_iniciais", exist=existe)
  IF (.NOT. existe) CALL criar_dir('valores_iniciais', dir)
END SUBROUTINE diretorio_vi

END MODULE