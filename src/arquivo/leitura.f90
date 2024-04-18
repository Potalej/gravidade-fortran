! *****************************************************************
!! LEITURA DE ARQUIVOS
!
! Objetivos:
!   Para facilitar a leitura e analise de arquivos, esta classe contem
!   alguns metodos "ajudadores" e outros prontos para uso.
! 
! Modificado:
!   15 de marco de 2024
! 
! Autoria:
!   oap
!
MODULE leitura

  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64

  IMPLICIT NONE
  PRIVATE
  PUBLIC preset_config

  ! Classe para leitura
  TYPE :: preset_config
    
    ! Integrais primeiras
    REAL(pf) :: Etot
    REAL(pf), DIMENSION(3) :: Jtot, Ptot, Rcm
    ! Constantes
    REAL(pf) :: G
    ! Intervalos para geracao
    REAL(pf), DIMENSION(2) :: int_massas, int_posicoes, int_momentos
    ! Integrador numerico
    CHARACTER(10) :: integrador
    ! Timestep
    REAL(pf) :: timestep
    ! Quantidade de passos antes de salvar
    INTEGER :: passos_antes_salvar = 1
    ! Quantidade de corpos
    INTEGER :: N
    ! Intervalo de integracao
    REAL(pf) :: t0, tf
    ! Uso de corretor e colisoes
    LOGICAL  :: corretor, colisoes

    ! Massas, posicoes e momentos
    REAL(pf), allocatable :: R(:,:), P(:,:), massas(:)
    ! Nome do problema de valores iniciais
    CHARACTER(50)         :: nome 

    CONTAINS
      PROCEDURE :: config, tratamento, valores_iniciais
  END TYPE

CONTAINS

! ************************************************************
!! Leitura de preset
!
! Objetivos:
!   Le um preset de condicoes iniciais para geracao (sorteio).
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE config (self, arquivo)
  CLASS(preset_config),  INTENT(INOUT) :: self
  CHARACTER(256), INTENT(INOUT) :: arquivo
  CHARACTER(len=48)             :: atributo, valor

  WRITE (*,'(a)') "LEITURA DE PRESET"
  WRITE (*,'(a)') '  > lendo o arquivo "' // TRIM(arquivo) // '"'

  OPEN(2,file=arquivo)
  READ(2,*) ! Comentario
  READ(2,*) ! Modo
  
  ! Integrais primeiras
  READ(2,*) atributo, self % Etot
  READ(2,*) atributo, self % Jtot
  READ(2,*) atributo, self % Ptot

  READ(2,*) ! Espaco
  READ(2,*) ! Comentario

  ! Constantes
  READ(2,*) atributo, self % G

  READ(2,*) ! Espaco
  READ(2,*) ! Comentario

  READ(2,*) atributo, self % N

  READ(2,*) ! Espaco
  READ(2,*) ! Comentario

  ! Intervalos para geracao
  READ(2,*) atributo, self % int_massas
  READ(2,*) atributo, self % int_posicoes
  READ(2,*) atributo, self % int_momentos

  READ(2,*) ! Espaco
  READ(2,*) ! Comentario  

  ! Integracao
  READ(2,*) atributo, self % integrador
  READ(2,*) atributo, self % timestep
  READ(2,*) atributo, self % passos_antes_salvar

  READ(2,*) ! Espaco
  READ(2,*) ! Comentario  

  READ(2,*) atributo, self % t0 ! tempo inicial
  READ(2,*) atributo, self % tf ! tempo final

  READ(2,*) ! Espaco
  READ(2,*) ! Comentario  

  READ(2,*) atributo, self % corretor ! uso do corretor
  READ(2,*) atributo, self % colisoes ! uso de colisoes

  CLOSE(2)

  WRITE (*,'(a)') '  > arquivo lido com sucesso!'
  WRITE (*,*)

END SUBROUTINE config

! ************************************************************
!! Tratamento dos dados.
!
! Objetivos:
!   Tratamento dos dados. Por enquanto, so verifica se o zero
!   esta contido, para padronizacao.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE tratamento (self)
  CLASS(preset_config),  INTENT(INOUT) :: self

  ! Se o tempo inicial for negativo e o final tambem, entao tem erro
  IF (self%t0 < 0 .AND. self%tf < 0) THEN
    ERROR STOP "O intervalo de tempo deve conter o zero!"
  ! Idem se o tempo inicial for positivo
  ELSE IF (self%t0 > 0) THEN
    ERROR STOP "O intervalo de tempo deve conter o zero!"
  ENDIF

END SUBROUTINE tratamento

! ************************************************************
!! Leitura de arquivo de valores iniciais
!
! Objetivos:
!   Faz a leitura de um arquivo com preset de valores iniciais,
!   ou seja, com valores de posicao e velocidades prontos.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE valores_iniciais (self, arquivo)
  CLASS(preset_config),  INTENT(INOUT) :: self
  CHARACTER(256), INTENT(INOUT) :: arquivo
  CHARACTER(len=48)             :: atributo, valor
  REAL(pf), DIMENSION(3)        :: R, P
  INTEGER                       :: a ! iterador

  WRITE (*,'(a)') "LEITURA DE VALORES INICIAIS"
  WRITE (*,'(a)') '  > lendo o arquivo "' // TRIM(arquivo) // '"'

  OPEN(2,file=arquivo)
  READ(2,*) ! Comentario ("Configs")
  READ(2,*) ! Modo

  READ(2,*) atributo, self % nome
  READ(2,*) atributo, self % integrador
  READ(2,*) atributo, self % timestep
  READ(2,*) atributo, self % passos_antes_salvar
  READ(2,*) atributo, self % t0 ! tempo inicial
  READ(2,*) atributo, self % tf ! tempo final
  READ(2,*) atributo, self % corretor
  READ(2,*) atributo, self % colisoes 

  READ(2,*) ! Espaco
  READ(2,*) ! Comentario

  ! Constantes
  READ(2,*) atributo, self % N
  READ(2,*) atributo, self % G

  READ(2,*) ! Espaco
  READ(2,*) ! Comentario

  ! Aloca os espacos na memoria
  ALLOCATE(self%massas(self%N))
  ALLOCATE(self%R(self%N,3))
  ALLOCATE(self%P(self%N,3))

  ! Leitura das massas
  DO a = 1, self%N
    READ(2,*) self%massas(a)
  END DO

  READ(2,*) ! Espaco
  READ(2,*) ! Comentario

  ! Leitura das posicoes
  DO a = 1, self%N
    READ(2,*) self%R(a,:)
  END DO

  READ(2,*) ! Espaco
  READ(2,*) ! Comentario

  ! Leitura dos momentos lineares
  DO a = 1, self%N
    READ(2,*) self%P(a,:)
  END DO

  CLOSE(2)

  WRITE (*,'(a)') '  > arquivo lido com sucesso!'
  WRITE (*,*)

END SUBROUTINE valores_iniciais
END MODULE leitura