module simulacao_vi

  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use OMP_LIB
  use simulacao
  use leitura
  implicit none
  private
  public simular_vi

  ! Instanciamento da classe
  type(simular) :: Sim_rk4, Sim_verlet, Sim_corrigir
  type(preset_config) :: configs

contains

  ! Metodo principal da classe
  subroutine simular_vi (arquivo)
    character(LEN=*), intent(inout) :: arquivo

    ! Tempo de execucao
    real :: t0, tf
    ! Quantidade total de passos
    integer :: qntd_total_passos

    ! Le o arquivo de configuracoes
    call configs % valores_iniciais(arquivo)

    WRITE (*,*)

    ! Se o instante inicial for negativo, entao vai rodar ao contrario
    if (configs%t0 < 0) then
      if (configs%tf == 0) then
        ! Roda apenas o passado
        WRITE (*,*) "Intervalo [", configs%t0, ",", configs%tf, "]"
        call rodar(-configs%timestep,configs%massas,configs%R,configs%P)
      else if (configs%tf > 0) then
        ! Roda o passado e o futuro
        WRITE (*,*) "Intervalo [", configs%t0, ",", 0, "]"
        call rodar(-configs%timestep,configs%massas,configs%R,configs%P)
        WRITE (*,*) "Intervalo [", 0, ",", configs%tf, "]"
        call rodar(configs%timestep,configs%massas,configs%R,configs%P)
      end if
    ! Se for positivo, apenas roda normal
    else
      ! Roda apenas o futuro
      WRITE (*,*) "Intervalo [", 0, ",", configs%tf, "]"
      call rodar(configs%timestep,configs%massas,configs%R,configs%P)
    end if

  end subroutine simular_vi

  subroutine rodar (timestep, massas, posicoes, momentos)
    REAL(pf), allocatable :: massas(:), posicoes(:,:), momentos(:,:)
    REAL(pf), intent(in)  :: timestep
    REAL(pf)              :: t0, tf
    INTEGER               :: qntd_total_passos

    qntd_total_passos = (configs%tf - configs%t0) / configs%timestep
    ! timer
    t0 = omp_get_wtime()

    WRITE (*,*) 'Rodando com ', qntd_total_passos, ' passos'

    SELECT CASE (configs%integrador)
      CASE ("verlet")
        Sim_verlet % corrigir = .FALSE.
        Sim_verlet % colidir  = .FALSE.
        call Sim_verlet%Iniciar(configs%G, massas, posicoes, momentos, timestep, configs%passos)
        call Sim_verlet%rodar_verlet(qntd_total_passos)
      ! Adicionar outros casos posteriormente
    END SELECT

    tf = omp_get_wtime()
    WRITE (*,*) '= Tempo ', configs%integrador, ': ', tf - t0
    WRITE (*,*) '= Tempor por passo: ', (tf-t0) / qntd_total_passos
  end subroutine rodar

end module simulacao_vi