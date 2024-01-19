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
    ! Uso de corretor e colisoes
    LOGICAL  :: corretor, colisoes

    ! Massas, posicoes e momentos
    REAL(pf), allocatable :: R(:,:), P(:,:), massas(:)
    ! Nome do problema de valores iniciais
    CHARACTER(50)         :: nome 


    contains
      procedure :: config, tratamento, valores_iniciais
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

  subroutine valores_iniciais (self, arquivo)
    class(preset_config),  intent(inout) :: self
    character(256), intent(inout) :: arquivo
    character(len=48)             :: atributo, valor
    real(pf), dimension(3)        :: R, P
    integer                       :: a ! iterador

    WRITE (*,*) "= Leitura do arquivo de valores iniciais"
    WRITE (*,*) '> Lendo o arquivo "', TRIM(arquivo), '"...'

    OPEN(2,file=arquivo)
    READ(2,*) ! Comentario ("Configs")
    
    READ(2,*) atributo, self % nome
    READ(2,*) atributo, self % integrador
    READ(2,*) atributo, self % timestep
    READ(2,*) atributo, self % passos
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
    allocate(self%massas(self%N))
    allocate(self%R(self%N,3))
    allocate(self%P(self%N,3))

    ! Leitura das massas
    do a = 1, self%N
      READ(2,*) self%massas(a)
    end do

    READ(2,*) ! Espaco
    READ(2,*) ! Comentario

    ! Leitura das posicoes
    do a = 1, self%N
      READ(2,*) self%R(a,:)
    end do

    READ(2,*) ! Espaco
    READ(2,*) ! Comentario

    ! Leitura dos momentos lineares
    do a = 1, self%N
      READ(2,*) self%P(a,:)
    end do

    CLOSE(2)

    WRITE (*,*) '> Arquivo lido com sucesso!'

  end subroutine valores_iniciais
end module