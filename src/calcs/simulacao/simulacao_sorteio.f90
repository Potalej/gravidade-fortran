module simulacao_sorteio

  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use OMP_LIB
  use simulacao
  use condicoesArtigo
  use leitura
  use arquivos
  implicit none
  private
  public simular_sorteio, sorteio_salvar

  ! Instanciamento da classe
  type(simular) :: Sim_rk4, Sim_verlet, Sim_corrigir
  type(preset_config) :: configs

contains

  ! Metodo principal da classe
  subroutine simular_sorteio (arquivo)
    character(256), intent(inout) :: arquivo

    ! Tempo de execucao
    real :: t0, tf
    ! Vetores
    real(pf), allocatable :: massas(:), posicoes(:,:), momentos(:,:)
    ! Quantidade total de passos
    integer :: qntd_total_passos
    CHARACTER(9) :: datahoje
    CHARACTER(3) :: numero
    CHARACTER(16) :: nome_arq
    CHARACTER(13) :: nome_sorteio
    INTEGER :: i = 1
    LOGICAL :: arquivo_existe = .TRUE.

    ! Le o arquivo de configuracoes
    call configs % config(arquivo)

    call gerar_condicionado(configs%G, &
       configs%N, &
       massas,    &
       posicoes,  &
       momentos,  &
       configs%int_posicoes, & ! Intervalo de posicoes
       configs%int_momentos, & ! Intervalo de momentos
       configs%int_massas)     ! Intervalo de massas

    ! Gera o nome
    call date_and_time(datahoje)

    do while (arquivo_existe)
      write(numero, '(I3.3)') i
      i = i + 1

      ! cria nome 
      nome_sorteio = trim(datahoje)//"_"//trim(numero)
      nome_arq = nome_sorteio // ".txt"

      ! verifica se existe
      inquire(file='out/auto_vi/'//nome_arq, exist=arquivo_existe)
    end do

    ! Salva o preset gerado
    CALL salvar_sorteio('out/', 'auto_vi/', nome_arq,  &
      "Sorteio_"//nome_sorteio, &
      configs % G,          &
      massas,               &
      posicoes,             &
      momentos,             &
      configs % t0,         &
      configs % tf,         &
      configs % timestep,   &
      configs % integrador, &
      configs % corretor,   &
      configs % colisoes)

    ! Se o instante inicial for negativo, entao vai rodar ao contrario
    if (configs%t0 < 0) then
      if (configs%tf == 0) then
        ! Roda apenas o passado
        WRITE (*,*) " > Intervalo [", configs%t0, ",", configs%tf, "]"
        call rodar(-configs%timestep,massas,posicoes,momentos)
      else if (configs%tf > 0) then
        ! Roda o passado e o futuro
        WRITE (*,*) " > Intervalo [", configs%t0, ",", 0, "]"
        call rodar(-configs%timestep,massas,posicoes,momentos)
        WRITE (*,*) " > Intervalo [", 0, ",", configs%tf, "]"
        call rodar(configs%timestep,massas,posicoes,momentos)
      end if
    ! Se for positivo, apenas roda normal
    else
      ! Roda apenas o futuro
      WRITE (*,*) " > Intervalo [", 0, ",", configs%tf, "]"
      call rodar(configs%timestep,massas,posicoes,momentos)
    end if

  end subroutine simular_sorteio

  subroutine rodar (timestep, massas, posicoes, momentos)
    REAL(pf), allocatable :: massas(:), posicoes(:,:), momentos(:,:)
    REAL(pf), intent(in)  :: timestep
    REAL(pf)              :: t0, tf
    INTEGER               :: qntd_total_passos

    qntd_total_passos = (configs%tf - configs%t0) / configs%timestep
    WRITE (*,*)

    ! timer
    t0 = omp_get_wtime()

    SELECT CASE (configs%integrador)
      CASE ("verlet")
        Sim_verlet % corrigir = .FALSE.
        Sim_verlet % colidir  = .FALSE.
        call Sim_verlet%Iniciar(configs%G, massas, posicoes, momentos, timestep, configs%passos)
        call Sim_verlet%rodar_verlet(qntd_total_passos)
      ! Adicionar outros casos posteriormente
    END SELECT

    tf = omp_get_wtime()
    WRITE (*,*) ' * tempo ', configs%integrador, ': ', tf - t0
    WRITE (*,*) ' * tempor por passo: ', (tf-t0) / qntd_total_passos
  end subroutine rodar

  subroutine sorteio_salvar (dir)
    IMPLICIT NONE
    CHARACTER(LEN=*) :: dir
    REAL(pf), ALLOCATABLE :: massas(:), posicoes(:,:), momentos(:,:)
    CHARACTER(8) :: datahoje
    CHARACTER(3) :: numero
    CHARACTER(16) :: nome_arq
    CHARACTER(12) :: nome_sorteio
    INTEGER :: i = 1
    LOGICAL :: arquivo_existe = .TRUE.

    ! em string
    call date_and_time(datahoje)

    ! Le o arquivo de configuracoes
    call configs % config(dir)

    ! Gera os valores
    call gerar_condicionado(configs%G, &
       configs%N, &
       massas,    &
       posicoes,  &
       momentos,  &
       configs%int_posicoes, & ! Intervalo de posicoes
       configs%int_momentos, & ! Intervalo de momentos
       configs%int_massas)     ! Intervalo de massas

    ! Gera o nome
    do while (arquivo_existe)
      write(numero, '(I3.3)') i
      i = i + 1

      ! cria nome 
      nome_sorteio = trim(datahoje)//"_"//trim(numero)
      nome_arq = nome_sorteio // ".txt"

      ! verifica se existe
      inquire(file="./out/auto_vi/"//nome_arq, exist=arquivo_existe)
    end do

    ! Salva o preset gerado
    CALL salvar_sorteio('out/','auto_vi/', nome_arq,  &
      "Sorteio_"//nome_sorteio, &
      configs % G,          &
      massas,               &
      posicoes,             &
      momentos,             &
      configs % t0,         &
      configs % tf,         &
      configs % timestep,   &
      configs % integrador, &
      configs % corretor,   &
      configs % colisoes)
  end subroutine sorteio_salvar

end module