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
  REAL(pf) :: V, T, V_vir

  ! variaveis de estatisticas
  REAL(pf), ALLOCATABLE :: E(:), E_err_rel(:), E_cin(:), E_pot(:)
  REAL(pf) :: E_final, E_medio
  REAL(pf) :: J(3), Qcm(3), J0(3), P0(3)
  REAL(pf), ALLOCATABLE :: Jx(:), Jy(:), Jz(:)
  REAL(pf), ALLOCATABLE :: Px(:), Py(:), Pz(:)
  REAL(pf), ALLOCATABLE :: Qcmx(:), Qcmy(:), Qcmz(:)
  REAL(pf), ALLOCATABLE :: r10(:), r25(:), r50(:), r70(:), r90(:) ! raios lagrangianos
  REAL(pf), ALLOCATABLE :: momine(:), momdil(:), virial(:)
  REAL(pf) :: f_prod_q
  INTEGER, ALLOCATABLE :: escapes(:)

  ! namelists
  NAMELIST /energia_nml/ E, E_err_rel, E_final, E_medio, E_cin, E_pot
  NAMELIST /integrais_nml/ Jx, Jy, Jz, Px, Py, Pz, Qcmx, Qcmy, Qcmz
  NAMELIST /dinamica_nml/ r10, r25, r50, r70, r90, momine, momdil, virial, escapes

  ! le o arquivo com a trajetoria
  CALL ler_data(arquivo_in, dt, G, eps, massas, Rs, Ps)

  CALL SYSTEM_CLOCK(count_rate=cron_rate)
  CALL SYSTEM_CLOCK(cron_0)

  ! alocando
  qntd = SIZE(Rs,1)
  ALLOCATE(E(qntd),  E_err_rel(qntd-1))
  ALLOCATE(E_cin(qntd), E_pot(qntd))
  ALLOCATE(Jx(qntd-1), Jy(qntd-1), Jz(qntd-1))
  ALLOCATE(Px(qntd-1), Py(qntd-1), Pz(qntd-1))
  ALLOCATE(Qcmx(qntd), Qcmy(qntd), Qcmz(qntd))
  ALLOCATE(r10(qntd), r25(qntd), r50(qntd), r70(qntd), r90(qntd))
  ALLOCATE(momine(qntd), momdil(qntd), virial(qntd))
  ALLOCATE(escapes(qntd))

  !> calculo das estatisticas iniciais
  CALL interacoes_mutuas(G, massas, Rs(1,:,:), Ps(1,:,:), eps, &
                        E_pot(1), f_prod_q, E_cin(1), escapes(1))

  E(1) = E_pot(1) + E_cin(1)
  J0 = momento_angular_total(Rs(1,:,:), Ps(1,:,:))
  P0(1) = SUM(Ps(1,:,1))
  P0(2) = SUM(Ps(1,:,2))
  P0(3) = SUM(Ps(1,:,3))
  ! Jx(1) = J(1); Jy(1) = J(2); Jz(1) = J(3)
  ! Jx(1) = 0.0_pf; Jy(1) = 0.0_pf; Jz(1) = 0.0_pf
  ! Px(1) = SUM(Ps(1,:,1))
  ! Py(1) = SUM(Ps(1,:,2))
  ! Pz(1) = SUM(Ps(1,:,3))
  Qcm = centro_massas(massas, Rs(1,:,:))
  Qcmx(1) = Qcm(1); Qcmy(1) = Qcm(2); Qcmz(1) = Qcm(3)
  
  r10(1) = raio_meia_massa(massas, Rs(1,:,:), 0.1_pf)
  r25(1) = raio_meia_massa(massas, Rs(1,:,:), 0.25_pf)
  r50(1) = raio_meia_massa(massas, Rs(1,:,:))
  r70(1) = raio_meia_massa(massas, Rs(1,:,:), 0.70_pf)
  r90(1) = raio_meia_massa(massas, Rs(1,:,:), 0.9_pf)

  IF (eps > 0) THEN
    virial(1) = 2.0_pf * E_cin(1) + f_prod_q
  ELSE
    virial(1) = 2.0_pf * E_cin(1) + E_pot(1)
  ENDIF

  momine(1) = momento_inercia_vec(massas, Rs(1,:,:))
  momdil(1) = momento_dilatacao(Rs(1,:,:), Ps(1,:,:))

  ! erro relativo
  DO i = 1, qntd - 1
    CALL interacoes_mutuas(G, massas, Rs(i+1,:,:), Ps(i+1,:,:), eps, &
                        E_pot(i+1), f_prod_q, E_cin(i+1), escapes(i+1))
    
    ! energia total
    E(i+1) = E_cin(i+1) + E_pot(i+1)
    E_err_rel(i) = ABS(E(i+1) - E(1))/ABS(E(1))

    ! angular
    J = momento_angular_total(Rs(i+1,:,:), Ps(i+1,:,:))
    Jx(i) = J(1) - J0(1)
    Jy(i) = J(2) - J0(2)
    Jz(i) = J(3) - J0(3)

    ! linear
    Px(i) = SUM(Ps(i+1,:,1)) - P0(1)
    Py(i) = SUM(Ps(i+1,:,2)) - P0(2)
    Pz(i) = SUM(Ps(i+1,:,3)) - P0(3)

    ! centro de massas
    Qcm = centro_massas(massas, Rs(i+1,:,:))
    Qcmx(i+1) = Qcm(1)
    Qcmy(i+1) = Qcm(2)
    Qcmz(i+1) = Qcm(3)

    ! raios lagrangianos
    r10(i+1) = raio_meia_massa(massas, Rs(i+1,:,:), 0.1_pf)
    r25(i+1) = raio_meia_massa(massas, Rs(i+1,:,:), 0.25_pf)
    r50(i+1) = raio_meia_massa(massas, Rs(i+1,:,:))
    r70(i+1) = raio_meia_massa(massas, Rs(i+1,:,:), 0.70_pf)
    r90(i+1) = raio_meia_massa(massas, Rs(i+1,:,:), 0.90_pf)
    
    IF (eps > 0) THEN
      virial(i+1) = 2.0_pf * E_cin(i+1) + f_prod_q
    ELSE
      virial(i+1) = 2.0_pf * E_cin(i+1) + E_pot(i+1)
    ENDIF

    momine(i+1) = momento_inercia_vec(massas, Rs(i+1,:,:))
    momdil(i+1) = momento_dilatacao(Rs(i+1,:,:), Ps(i+1,:,:))
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

SUBROUTINE interacoes_mutuas (G, m, R, P, eps, pot_geral, f_prod_q, E_cin, escapes)
  REAL(pf), INTENT(IN) :: G, m(:), R(:,:), P(:,:), eps
  REAL(pf), INTENT(OUT) :: pot_geral, f_prod_q, E_cin
  INTEGER,  INTENT(OUT) :: escapes

  REAL(pf) :: Fab(3), dist2, denominador, ep, ec, pot_ind(SIZE(m))
  INTEGER :: a, b
  
  f_prod_q = 0.0_pf
  pot_geral = 0.0_pf
  pot_ind = 0.0_pf

  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP SHARED(m, R, P, G, eps) &
  !$OMP PRIVATE(a, b, Fab, dist2, denominador, ep) &
  !$OMP SCHEDULE(DYNAMIC) REDUCTION(+:pot_geral, pot_ind, f_prod_q)
  DO a = 2, SIZE(m)
    DO b = 1, a - 1
      Fab = R(b,:) - R(a,:)
      dist2 = DOT_PRODUCT(Fab, Fab)
      denominador = SQRT(dist2 + eps*eps)

      ! potencial
      ep = - G * m(a) * m(b) / denominador
      pot_geral  = pot_geral  + ep
      pot_ind(a) = pot_ind(a) + ep
      pot_ind(b) = pot_ind(b) + ep

      ! forca
      Fab = G * m(a) * m(b) * Fab / (denominador*denominador*denominador)

      ! virial
      f_prod_q = f_prod_q + DOT_PRODUCT(R(a,:), Fab)
      f_prod_q = f_prod_q - DOT_PRODUCT(R(b,:), Fab)
    END DO
  END DO
  !$OMP END PARALLEL DO

  escapes = 0
  E_cin = 0.0_pf
  DO a = 1, SIZE(m)
    ec = 0.5_pf * DOT_PRODUCT(P(a,:), P(a,:)) / m(a)
    E_cin = E_cin + ec
    IF (2.0_pf * ec > - pot_ind(a)) escapes = escapes + 1
  END DO
END SUBROUTINE


END MODULE