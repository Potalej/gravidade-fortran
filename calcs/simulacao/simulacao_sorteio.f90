module simulacao_sorteio

  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use OMP_LIB
  use simulacao
  use condicoesArtigo
  use leitura
  implicit none
  private
  public simular_sorteio

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

    ! Le o arquivo de configuracoes
    call configs % config(arquivo)

    WRITE (*,*)

    call gerar_condicionado(configs%G, &
       configs%N, &
       massas,    &
       posicoes,  &
       momentos,  &
       configs%int_posicoes, & ! Intervalo de posicoes
       configs%int_momentos, & ! Intervalo de momentos
       configs%int_massas)     ! Intervalo de massas

    WRITE (*,*)

    ! Se o instante inicial for negativo, entao vai rodar ao contrario
    if (configs%t0 < 0) then
      if (configs%tf == 0) then
        ! Roda apenas o passado
        WRITE (*,*) "Intervalo [", configs%t0, ",", configs%tf, "]"
        call rodar(-configs%timestep,massas,posicoes,momentos)
      else if (configs%tf > 0) then
        ! Roda o passado e o futuro
        WRITE (*,*) "Intervalo [", configs%t0, ",", 0, "]"
        call rodar(-configs%timestep,massas,posicoes,momentos)
        WRITE (*,*) "Intervalo [", 0, ",", configs%tf, "]"
        call rodar(configs%timestep,massas,posicoes,momentos)
      end if
    end if

  end subroutine simular_sorteio

  subroutine rodar (timestep, massas, posicoes, momentos)
    REAL(pf), allocatable :: massas(:), posicoes(:,:), momentos(:,:)
    REAL(pf), intent(in)  :: timestep
    REAL(pf)              :: t0, tf
    INTEGER               :: qntd_total_passos

    qntd_total_passos = configs%passos / configs%timestep
    ! timer
    t0 = omp_get_wtime()

    SELECT CASE (configs%integrador)
      CASE ("verlet")
        Sim_verlet % corrigir = .FALSE.
        Sim_verlet % colidir  = .FALSE.
        call Sim_verlet%Iniciar(configs%G, massas, posicoes, momentos, timestep, configs%passos)
        call Sim_verlet%rodar_verlet(qntd_total_passos)
        STOP
      ! Adicionar outros casos posteriormente
    END SELECT

    tf = omp_get_wtime()
    WRITE (*,*) '= Tempo ', configs%integrador, ': ', tf - t0
    WRITE (*,*) '= Tempor por passo: ', (tf-t0) / qntd_total_passos
  end subroutine rodar

end module