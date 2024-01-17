! Leitura de arquivos
!
! Para facilitar a leitura e analise de arquivos, esta classe contem alguns
! metodos "ajudadores" e outros prontos para uso.
!
! = function config (self, arquivo)
! Faz a leitura de um arquivo de presets para o sorteio de condicoes iniciais.
module leitura

  use, intrinsic :: iso_fortran_env, only: pf=>real64

  implicit none
  private
  public preset_config

  ! Classe para leitura
  type :: preset_config
    
    ! Integrais primeiras
    REAL(pf) :: Etot
    REAL(pf), dimension(3) :: Jtot, Ptot, Rcm
    ! Constantes
    REAL(pf) :: G
    ! Intervalos para geracao
    REAL(pf), dimension(2) :: int_massas, int_posicoes, int_momentos
    ! Integrador numerico
    CHARACTER(10) :: integrador
    ! Timestep
    REAL(pf) :: timestep
    ! Quantidade de passos
    INTEGER :: passos
    ! Quantidade de corpos
    INTEGER :: N
    ! Intervalo de integracao
    REAL(pf) :: t0, tf

    contains
      procedure :: config, tratamento
  end type

contains

  subroutine config (self, arquivo)
    class(preset_config),  intent(inout) :: self
    character(256), intent(inout) :: arquivo
    character(len=48)             :: atributo, valor

    WRITE (*,*) "= Leitura do arquivo de presets"
    WRITE (*,*) '> Lendo o arquivo "', TRIM(arquivo), '"...'

    OPEN(2,file=arquivo)
    READ(2,*) ! Comentario
    
    ! Integrais primeiras
    READ(2,*) atributo, self % Etot
    READ(2,*) atributo, self % Jtot
    READ(2,*) atributo, self % Ptot
    READ(2,*) atributo, self % Rcm

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
    READ(2,*) atributo, self % passos

    READ(2,*) ! Espaco
    READ(2,*) ! Comentario  

    READ(2,*) atributo, self % t0 ! tempo inicial
    READ(2,*) atributo, self % tf ! tempo final

    CLOSE(2)

    WRITE (*,*) '> Arquivo lido com sucesso!'

  end subroutine config

  subroutine tratamento (self)
    class(preset_config),  intent(inout) :: self

    ! Se o tempo inicial for negativo e o final tambem, entao tem erro
    if (self%t0 < 0 .AND. self%tf < 0) then
      ERROR STOP "O intervalo de tempo deve conter o zero!"
    ! Idem se o tempo inicial for positivo
    else if (self%t0 > 0) then
      ERROR STOP "O intervalo de tempo deve conter o zero!"
    end if

  end subroutine tratamento

end module