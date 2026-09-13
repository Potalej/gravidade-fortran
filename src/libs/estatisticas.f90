! *****************************************************************
!! ESTATISTICAS 
!
! Objetivos:
!   Calcula estatisticas sobre uma determinada simulacao e exporta
!   em um arquivo.
! 
! Modificado:
!   13 de setembro de 2026 (criado)
!   13 de setembro de 2026 (modificado)
! 
! Autoria:
!   oap
! 
MODULE estatisticas_mod
  USE tipos
  USE utilidades
  USE arquivos_mod
CONTAINS    

SUBROUTINE estatisticas_trajetoria (arquivo_in, arquivo_out)
  CHARACTER(LEN=*), INTENT(IN) :: arquivo_in, arquivo_out
  
  ! variaveis gerais
  REAL(pf) :: dt, G, eps
  REAL(pf), ALLOCATABLE :: massas(:), Rs(:,:,:), Ps(:,:,:)
  INTEGER :: i, qntd
  INTEGER :: cron_0, cron_1, cron_rate

  ! variaveis de estatisticas
  REAL(pf), ALLOCATABLE :: E(:), E_err_rel(:)
  REAL(pf) :: E_final, E_medio
  REAL(pf) :: J(3), Qcm(3)
  REAL(pf), ALLOCATABLE :: Jx(:), Jy(:), Jz(:)
  REAL(pf), ALLOCATABLE :: Px(:), Py(:), Pz(:)
  REAL(pf), ALLOCATABLE :: Qcmx(:), Qcmy(:), Qcmz(:)
  REAL(pf), ALLOCATABLE :: rmm(:) ! raio de meia massa

  ! namelists
  NAMELIST /energia_nml/ E, E_err_rel, E_final, E_medio
  NAMELIST /integrais_nml/ Jx, Jy, Jz, Px, Py, Pz, Qcmx, Qcmy, Qcmz
  NAMELIST /dinamica_nml/ rmm

  ! le o arquivo com a trajetoria
  CALL ler_data(arquivo_in, dt, G, eps, massas, Rs, Ps)

  CALL SYSTEM_CLOCK(count_rate=cron_rate)
  CALL SYSTEM_CLOCK(cron_0)

  ! alocando
  qntd = SIZE(Rs,1)
  ALLOCATE(E(qntd),  E_err_rel(qntd-1))
  ALLOCATE(Jx(qntd), Jy(qntd), Jz(qntd))
  ALLOCATE(Px(qntd), Py(qntd), Pz(qntd))
  ALLOCATE(Qcmx(qntd), Qcmy(qntd), Qcmz(qntd))
  ALLOCATE(rmm(qntd))

  !> calculo das estatisticas iniciais
  E(1) = energia_total_par(G, massas, Rs(1,:,:), Ps(1,:,:), eps)
  J = momento_angular_total(Rs(1,:,:), Ps(1,:,:))
  Jx(1) = J(1); Jy(1) = J(2); Jz(1) = J(3)
  Px(1) = SUM(Ps(1,:,1))
  Py(1) = SUM(Ps(1,:,2))
  Pz(1) = SUM(Ps(1,:,3))
  Qcm = centro_massas(massas, Rs(1,:,:))
  Qcmx(1) = Qcm(1); Qcmy(1) = Qcm(2); Qcmz(1) = Qcm(3)
  rmm(1) = raio_meia_massa(massas, Rs(1,:,:))

  ! erro relativo
  DO i = 1, qntd - 1
    ! energia total
    E(i+1) = energia_total_par(G, massas, Rs(i+1,:,:), Ps(i+1,:,:), eps)
    E_err_rel(i) = ABS(E(i+1) - E(1))/ABS(E(1))

    ! angular
    J = momento_angular_total(Rs(i+1,:,:), Ps(i+1,:,:))
    Jx(i+1) = J(1)
    Jy(i+1) = J(2)
    Jz(i+1) = J(3)

    ! linear
    Px(i+1) = SUM(Ps(i+1,:,1))
    Py(i+1) = SUM(Ps(i+1,:,2))
    Pz(i+1) = SUM(Ps(i+1,:,3))

    ! centro de massas
    Qcm = centro_massas(massas, Rs(i+1,:,:))
    Qcmx(i+1) = Qcm(1)
    Qcmy(i+1) = Qcm(2)
    Qcmz(i+1) = Qcm(3)

    ! raio de meia massa
    rmm(i+1) = raio_meia_massa(massas, Rs(i+1,:,:))
  END DO

  E_final = E(qntd)
  E_medio = SUM(E) / qntd

  CALL SYSTEM_CLOCK(cron_1)
  WRITE (*,'(a,F10.4,a)') "  > tempo estatisticas: ", REAL(cron_1-cron_0)/REAL(cron_rate), "s" 

  ! armazena em um arquivo
  OPEN(unit=10, file=arquivo_out, status="replace", action="write")
  WRITE(unit=10, nml=energia_nml)
  WRITE(unit=10, nml=integrais_nml)
  WRITE(unit=10, nml=dinamica_nml)
  CLOSE(unit=10)
END SUBROUTINE

END MODULE