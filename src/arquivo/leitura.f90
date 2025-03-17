! *****************************************************************
!! LEITURA DE ARQUIVOS
!
! Objetivos:
!   Para facilitar a leitura e analise de arquivos, esta classe contem
!   alguns metodos "ajudadores" e outros prontos para uso.
! 
! Modificado:
!   16 de marco de 2025
! 
! Autoria:
!   oap
!
MODULE leitura

  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  USE arquivos
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
    REAL(pf), DIMENSION(2) :: int_massas, int_momentos
    REAL(pf), DIMENSION(3) :: int_posicoes
    ! Integrador numerico
    CHARACTER(10) :: integrador
    ! Modo
    CHARACTER(20) :: modo
    ! Timestep
    REAL(pf) :: timestep
    ! Quantidade de passos antes de salvar
    INTEGER :: passos_antes_salvar = 1
    ! Quantidade de corpos
    INTEGER :: N
    ! Intervalo de integracao
    INTEGER :: t0, tf
    ! Softening ("amolecimento") do potencial 
    REAL(pf) :: potsoft
    ! Uso de corretor
    LOGICAL  :: corretor
    REAL(pf) :: corretor_margem_erro
    INTEGER  :: corretor_max_num_tentativas

    ! Uso de colisoes
    LOGICAL        :: colisoes
    CHARACTER(10)  :: colisoes_modo
    REAL(pf) :: colisoes_max_distancia

    ! Uso de paralelisacao
    LOGICAL  :: paralelo

    ! Geracao dos valores iniciais
    CHARACTER(20) :: vi_dist
    CHARACTER(20) :: vi_regiao
    REAL(pf)      :: vi_raio

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
!   03 de fevereiro de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE config (self, arquivo)
  CLASS(preset_config),  INTENT(INOUT) :: self
  CHARACTER(256), INTENT(INOUT) :: arquivo
  CHARACTER(len=48)             :: atributo, valor
  INTEGER(kind=4)               :: unidade

  call capturar_unidade(unidade)

  WRITE (*,'(a)') "LEITURA DE PRESET"
  WRITE (*,'(a)') '  > lendo o arquivo "' // TRIM(arquivo) // '"'

  OPEN(unidade,file=arquivo)
  READ(unidade,*) atributo, self % modo
  READ(unidade,*) ! Espaco
  READ(unidade,*) ! Comentario
  
  ! Integrais primeiras
  READ(unidade,*) atributo, self % Etot
  READ(unidade,*) atributo, self % Jtot
  READ(unidade,*) atributo, self % Ptot

  READ(unidade,*) ! Espaco
  READ(unidade,*) ! Comentario

  ! Constantes
  READ(unidade,*) atributo, self % G

  READ(unidade,*) ! Espaco
  READ(unidade,*) ! Comentario

  READ(unidade,*) atributo, self % N

  READ(unidade,*) ! Espaco
  READ(unidade,*) ! Comentario

  ! Intervalos para geracao
  READ(unidade,*) atributo, self % int_massas
  READ(unidade,*) atributo, self % int_posicoes
  READ(unidade,*) atributo, self % int_momentos

  READ(unidade,*) ! Espaco
  READ(unidade,*) ! Comentario  

  ! Integracao
  READ(unidade,*) atributo, self % integrador
  READ(unidade,*) atributo, self % timestep
  READ(unidade,*) atributo, self % potsoft
  READ(unidade,*) atributo, self % passos_antes_salvar

  READ(unidade,*) ! Espaco
  READ(unidade,*) ! Comentario  

  READ(unidade,*) atributo, self % t0 ! tempo inicial
  READ(unidade,*) atributo, self % tf ! tempo final

  READ(unidade,*) ! Espaco
  READ(unidade,*) ! Comentario  

  ! Uso de paralelisacao
  READ(unidade,*) atributo, self % paralelo

  READ(unidade,*) ! Espaco
  READ(unidade,*) ! Comentario  

  ! Correcao numerica
  READ(unidade,*) atributo, self % corretor
  READ(unidade,*) atributo, self % corretor_margem_erro
  READ(unidade,*) atributo, self % corretor_max_num_tentativas

  READ(unidade,*) ! Espaco
  READ(unidade,*) ! Comentario  

  ! Colisoes
  READ(unidade,*) atributo, self % colisoes_modo
  IF (TRIM(self % colisoes_modo) == 'F') THEN
    self % colisoes = .FALSE.
  ELSE
    self % colisoes = .TRUE.
  ENDIF
  READ(unidade,*) atributo, self % colisoes_max_distancia

  READ(unidade,*) ! Espaco
  READ(unidade,*) ! Comentario  

  ! Valores iniciais
  READ(unidade,*) atributo, self % vi_dist
  READ(unidade,*) atributo, self % vi_regiao
  READ(unidade,*) atributo, self % vi_raio

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
  INTEGER(kind=4)               :: unidade

  call capturar_unidade(unidade)

  WRITE (*,'(a)') "LEITURA DE VALORES INICIAIS"
  WRITE (*,'(a)') '  > lendo o arquivo "' // TRIM(arquivo) // '"'

  OPEN(unidade,file=arquivo)
  READ(unidade,*) ! Comentario ("Configs")
  READ(unidade,*) atributo, self % modo

  READ(unidade,*) atributo, self % nome
  READ(unidade,*) atributo, self % integrador
  READ(unidade,*) atributo, self % timestep
  READ(unidade,*) atributo, self % potsoft
  READ(unidade,*) atributo, self % passos_antes_salvar
  READ(unidade,*) atributo, self % t0 ! tempo inicial
  READ(unidade,*) atributo, self % tf ! tempo final
  
  READ(unidade,*) ! Espaco
  READ(unidade,*) ! Comentario

  ! Uso de paralelisacao
  READ(unidade,*) atributo, self % paralelo

  READ(unidade,*) ! Espaco
  READ(unidade,*) ! Comentario

  ! Correcao
  READ(unidade,*) atributo, self % corretor
  READ(unidade,*) atributo, self % corretor_margem_erro
  READ(unidade,*) atributo, self % corretor_max_num_tentativas

  READ(unidade,*) ! Espaco
  READ(unidade,*) ! Comentario

  ! Colisoes
  READ(unidade,*) atributo, self % colisoes_modo
  IF (TRIM(self % colisoes_modo) == 'F') THEN
    self % colisoes = .FALSE.
  ELSE
    self % colisoes = .TRUE.
  ENDIF
  READ(unidade,*) atributo, self % colisoes_max_distancia

  READ(unidade,*) ! Espaco
  READ(unidade,*) ! Comentario

  ! Constantes
  READ(unidade,*) atributo, self % N
  READ(unidade,*) atributo, self % G

  READ(unidade,*) ! Espaco
  READ(unidade,*) ! Comentario

  ! Aloca os espacos na memoria
  ALLOCATE(self%massas(self%N))
  ALLOCATE(self%R(self%N,3))
  ALLOCATE(self%P(self%N,3))

  ! Leitura das massas
  DO a = 1, self%N
    READ(unidade,*) self%massas(a)
  END DO

  READ(unidade,*) ! Espaco
  READ(unidade,*) ! Comentario

  ! Leitura das posicoes
  DO a = 1, self%N
    READ(unidade,*) self%R(a,:)
  END DO

  READ(unidade,*) ! Espaco
  READ(unidade,*) ! Comentario

  ! Leitura dos momentos lineares
  DO a = 1, self%N
    READ(unidade,*) self%P(a,:)
  END DO

  CLOSE(unidade)

  WRITE (*,'(a)') '  > arquivo lido com sucesso!'
  WRITE (*,*)

END SUBROUTINE valores_iniciais
END MODULE leitura