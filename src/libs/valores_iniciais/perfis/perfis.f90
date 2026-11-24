! ************************************************************
!! Initial values
!
!> Objectives
!  This module provides routines to generate random initial
!  values. The available density profiles are:
!   - uniform / uniform_cubic (default)
!   - uniform_spherical
!   - homogeneous_sphere
!   - plummer   (Plummer, 1911)
!   - hernquist (Hernquist, 1990)
!
!  The velocities are sampled using p(v|r) \propto v^2 f, where
!  f can be an ergodic distribution function (if the systems
!  are isotropic) or other distribution function. The available
!  options for this are. For the Plummer and the Hernquist we
!  sample isotropic velocities, and the velocities are radially
!  uniform for the others.
!
!> Modified
!  2026.11.24
!
!> Created
!  2026.07.28
!
!> Author
!  oap
!
MODULE initial_values_mod
    USE iv_plummer_mod
    USE iv_homogeneous_sphere_mod
    USE iv_hernquist_mod
    USE tipos
    IMPLICIT NONE
    PUBLIC generate_initial_values, api
    PRIVATE

    REAL(pf), PARAMETER :: PI = 4.0_pf * ATAN(1.0_pf)
CONTAINS

SUBROUTINE api (size_qs, size_ps, nm, td, pr, pp, m, qs, ps)
    REAL(pf), INTENT(IN) :: size_qs, size_ps
    LOGICAL, INTENT(IN) :: nm ! normalized masses
    LOGICAL, INTENT(IN) :: td ! two dimensional
    CHARACTER(LEN=*), INTENT(IN) :: pr ! mass profile
    REAL(pf),          INTENT(IN) :: pp(:) ! profile parameters
    REAL(pf), INTENT(INOUT) :: m(:), qs(:,:), ps(:,:)

    CALL generate_initial_values(SIZE(m), m, qs, ps, size_qs, size_ps, nm, td, pr, pp)
END SUBROUTINE

SUBROUTINE generate_initial_values (N, m, qs, ps, size_qs, size_ps, nm, td, pr, pp, seed)
    INTEGER,  INTENT(IN) :: N
    REAL(pf), INTENT(IN) :: size_qs, size_ps
    LOGICAL,  INTENT(IN) :: nm ! normalized masses
    LOGICAL,  INTENT(IN) :: td ! two dimensional
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: pr ! mass profile
    REAL(pf), INTENT(IN), OPTIONAL         :: pp(:) ! profile parameters
    INTEGER, INTENT(IN), OPTIONAL         :: seed
    REAL(pf), INTENT(INOUT) :: m(N), qs(3,N), ps(3,N)
    REAL(pf) :: r(N), vr(N) ! radius and radial velocities

    REAL(pf) :: tlm(3), com(3), big_ps, max_ps(3)
    INTEGER :: p

    INTEGER :: seed_size
    INTEGER, ALLOCATABLE :: seed_vec(:)

    IF (PRESENT(seed)) THEN
        IF (seed .NE. -1) THEN
            CALL RANDOM_SEED(size=seed_size)
            ALLOCATE(seed_vec(seed_size))
            seed_vec = seed
            CALL RANDOM_SEED(put=seed_vec)
        ENDIF
    ENDIF

    ! the massas can be equal or no
    IF (nm) THEN 
        m = 1.0_pf / N
    ELSE
        CALL random_number(m)
        DO WHILE (MINVAL(m) == 0.0_pf)
            CALL random_number(m)
        END DO
        m = m / SUM(m)
    ENDIF

    ! if profile not present, generate uniform (cubic)
    IF (.NOT. PRESENT(pr)) THEN
        CALL generate_initial_values_uniform(N, qs, ps, size_qs, size_ps)
        RETURN
    ENDIF

    SELECT CASE (TRIM(pr))
        ! uniform (cubic), the velocities are cubic too
        CASE ("uniform")
            CALL generate_initial_values_uniform(N, qs, ps, size_qs, size_ps)

        ! Uniform distribution of the radius, generating a uniform sphere.
        ! The velocities are isotropic.
        CASE ("uniform_sphere")
            CALL generate_initial_values_uniform_spheric(N, r, vr, size_qs, size_ps)
            CALL generate_cartesian_from_radius(N, r,  qs, td)
            CALL generate_cartesian_from_radius(N, vr, ps, td)
        
        ! Plummer sphere (Plummer 1911). The velocities are isotropic.
        CASE ("plummer")
            IF (.NOT. PRESENT(pp) .OR. SIZE(pp) .NE. 2) THEN
                STOP "For the Plummer profile, it is necessary to set the parameters G and b."
            ENDIF
            CALL generate_initial_values_plummer(N, r, vr, size_qs, pp)
            CALL generate_cartesian_from_radius(N, r,  qs, td)
            CALL generate_cartesian_from_radius(N, vr, ps, td)

        ! Homogeneous sphere, with isotropic velocities.
        CASE ("homogeneous")
            IF (.NOT. PRESENT(pp) .OR. SIZE(pp) .NE. 2) THEN
                STOP "For the homogeneous sphere profile, it is necessary to set the parameters G and R."
            ENDIF
            CALL generate_initial_values_homogeneous_sphere(N, r, vr, size_qs, size_ps, pp)
            CALL generate_cartesian_from_radius(N, r,  qs, td)
            CALL generate_cartesian_from_radius(N, vr, ps, td)

        ! Hernquist model (Hernquist 1990), with isotropic velocities
        CASE ("hernquist")
            IF (.NOT. PRESENT(pp) .OR. SIZE(pp) .NE. 2) THEN
                STOP "For the Hernquist sphere profile, it is necessary to set the parameters G and a"
            ENDIF
            CALL generate_initial_values_hernquist_isotropic(N, r, vr, size_qs, pp)
            CALL generate_cartesian_from_radius(N, r,  qs, td)
            CALL generate_cartesian_from_radius(N, vr, ps, td)

        CASE DEFAULT
            STOP 'Unknow density profile: "'//TRIM(pr)//'"'
    END SELECT

    IF (size_ps >= 0) THEN
        max_ps(1) = MAXVAL(ABS(ps(1,:)))
        max_ps(2) = MAXVAL(ABS(ps(2,:)))
        max_ps(3) = MAXVAL(ABS(ps(3,:)))
        big_ps = MAXVAL(max_ps)

        IF (big_ps > size_ps) ps = ps * size_ps / big_ps
    ENDIF

    !> 2d
    IF (td) THEN
        ps(3,:) = 0.0_pf
        qs(3,:) = 0.0_pf
    ENDIF

    ! velocities to momentum
    DO p = 1, N
        ps(:,p) = ps(:,p) * m(p)
    END DO

    ! total linear momentum -> 0
    tlm(1) = sum(ps(1,:))/N
    tlm(2) = sum(ps(2,:))/N
    tlm(3) = sum(ps(3,:))/N
    DO p = 1, N
        ps(:,p) = ps(:,p) - tlm
    END DO

    ! center of mass to the origin
    com = 0.0_pf
    DO p = 1, N
        com = com + m(p) * qs(:,p)
    END DO
    com = com / SUM(m)

    DO p = 1, N
        qs(:,p) = qs(:,p) - com
    END DO
END SUBROUTINE

SUBROUTINE generate_initial_values_uniform (N, qs, ps, size_qs, size_ps)
! Generates random initial values with uniform distribution (positions
! and momenta). The masses can be 1/N if nm == .TRUE. or uniformly
! distributed with M = 1.
    INTEGER, INTENT(IN)    :: N
    REAL(pf), INTENT(IN)    :: size_qs, size_ps
    REAL(pf), INTENT(INOUT) :: qs(3,N), ps(3,N)

    CALL random_number(qs)
    qs = size_qs * (2.0_pf * qs - 1.0_pf)

    CALL random_number(ps)
    ps = size_ps * (2.0_pf * ps - 1.0_pf)
END SUBROUTINE

SUBROUTINE generate_initial_values_uniform_spheric (N, r, vr, size_qs, size_ps)
! Generates random radius with uniform distribution. The masses can 
! be 1/N if nm == .TRUE. or uniformly distributed with M = 1.
    INTEGER, INTENT(IN)    :: N
    REAL(pf), INTENT(IN)    :: size_qs, size_ps
    REAL(pf), INTENT(INOUT) :: r(N), vr(N)

    CALL random_number(r)
    r = size_qs * (2.0_pf * r - 1.0_pf)

    CALL random_number(vr)
    vr = size_ps * (2.0_pf * vr - 1.0_pf)
END SUBROUTINE

SUBROUTINE generate_initial_values_plummer (N, r, vr, size_qs, pars)
! Generates random initial values with Plummer density and isotropic velocities. 
! The masses can be 1/N if nm == .TRUE. or uniformly distributed with M = 1.
    INTEGER, INTENT(IN)    :: N
    REAL(pf), INTENT(IN)    :: pars(2) ! G and plummer parameter
    REAL(pf), INTENT(IN)    :: size_qs
    REAL(pf), INTENT(INOUT) :: r(N), vr(N)

    REAL(pf) :: X
    
    INTEGER  :: i
    REAL(pf) :: potential, v_escape
    REAL(pf) :: G, plupar
    REAL(pf) :: qvn, yvn, edf, maximum

    ! parameters
    G = pars(1)
    plupar = pars(2)

    !> POSITIONS
    i = 1
    DO WHILE (i <= N)
        ! uniformly distributed values
        CALL random_number(X)

        ! apply P^-1(u) = r
        ! r(i) = X**(1.0_pf/3.0_pf) * plupar / (SQRT(1.0_pf - X**(2.0_pf/3.0_pf)))
        r(i) = plummer_inverse_P(X, pars)
        
        IF (size_qs >= 0 .AND. r(i) > size_qs) CYCLE
        i = i + 1
    END DO

    !> MOMENTA
    DO i = 1, N
        ! potential = - G * 1.0_pf / SQRT(r(i)**2 + plupar**2)
        potential = plummer_potential(r(i), pars)
        v_escape = SQRT(-2.0_pf * potential)
        ! maximum = 0.1_pf
        maximum = plummer_maximum_pdf()

        von_neumann: DO WHILE (.TRUE.)
            CALL random_number(qvn)
            CALL random_number(yvn)
            yvn = maximum * yvn

            edf = plummer_edf_norm(qvn)

            ! IF (yvn < qvn*qvn*(1.0_pf - qvn*qvn)**(3.5_pf)) EXIT von_neumann
            IF (yvn < qvn*qvn*edf) EXIT von_neumann
        END DO von_neumann

        vr(i) = qvn * v_escape
    END DO
END SUBROUTINE

SUBROUTINE generate_initial_values_homogeneous_sphere (N, r, vr, size_qs, size_ps, pars)
! Generates random initial values with homogeneous density and uniform velocities.
! The masses can be 1/N if nm == .TRUE. or uniformly distributed with M = 1.
    INTEGER, INTENT(IN)    :: N
    REAL(pf), INTENT(IN)    :: pars(2) ! G, R
    REAL(pf), INTENT(IN)    :: size_qs, size_ps
    REAL(pf), INTENT(INOUT) :: r(N), vr(N)

    REAL(pf) :: X
    INTEGER  :: i

    !> POSITIONS
    i = 1
    DO WHILE (i <= N)
        ! uniformly distributed values
        CALL random_number(X)

        ! apply P^-1(u) = r
        r(i) = homogeneous_sphere_inverse_P(x, pars)
        IF (size_qs >= 0 .AND. r(i) > size_qs) CYCLE

        i = i + 1
    END DO

    !> MOMENTA
    CALL random_number(vr)
    vr = size_ps * (2.0_pf * vr - 1.0_pf)
END SUBROUTINE

SUBROUTINE generate_initial_values_hernquist_isotropic (N, r, vr, size_qs, pars)
! Generates random initial values with Hernquist density and isotropic velocities.
! The masses can be 1/N if nm == .TRUE. or uniformly distributed with M = 1.
    INTEGER, INTENT(IN)    :: N
    REAL(pf), INTENT(IN)    :: pars(3) ! G, R and rho0
    REAL(pf), INTENT(IN)    :: size_qs
    REAL(pf), INTENT(INOUT) :: r(N), vr(N)

    REAL(pf) :: X
    
    INTEGER :: i
    REAL(pf) :: potential, v_escape
    REAL(pf) :: G, apar ! parameters of the profile

    REAL(pf) :: qvn, yvn
    REAL(pf) :: edf

    INTEGER :: counter
    REAL(pf) :: max_qvn

    ! parameters
    G = pars(1)
    apar = pars(2)

    !> POSITIONS
    i = 1
    DO WHILE (i <= N)
        ! uniformly distributed values
        CALL random_number(X)

        ! apply P^-1(u) = r
        r(i) = hernquist_inverse_P(x, pars)
        IF (size_qs >= 0 .AND. r(i) > size_qs) CYCLE

        i = i + 1
    END DO

    !> MOMENTA
    DO i = 1, N
        potential = hernquist_potential(r(i), pars)
        v_escape = SQRT(-2.0_pf * potential)

        ! get the maximum of the distribution v^2 f(E)
        max_qvn = hernquist_maximum_pdf(r(i), pars)

        ! TODO
        counter = 0
        von_neumann: DO WHILE (.TRUE.)
            CALL random_number(yvn)
            IF (counter == 100) THEN
                STOP "Too many attempts to sample with von Neumann, try again"
            ENDIF
            counter = counter + 1
            CALL random_number(qvn)
            yvn = max_qvn * yvn

            edf = hernquist_edf_norm(r(i), qvn, pars)
            ! print *, 'aqui: ', yvn,SQRT(-potential*(1- qvn**2)),  qvn*qvn*fE

            IF (yvn < qvn*qvn*edf) EXIT von_neumann
        END DO von_neumann

        vr(i) = v_escape * qvn
    END DO
END SUBROUTINE

SUBROUTINE generate_cartesian_from_radius (N, r, qs, td)
    INTEGER, INTENT(IN)  :: N
    REAL(pf), INTENT(IN)  :: r(N)
    REAL(pf), INTENT(OUT) :: qs(3,N)
    LOGICAL, INTENT(IN)  :: td ! 2d
    REAL(pf) :: X(N), phi(N), cos_theta(N), sin_theta(N)

    ! phi angle variables
    CALL random_number(X)
    phi = 2.0_pf * PI * X

    ! theta angle variables
    ! if two dimensional, theta = PI
    IF (td) THEN
        cos_theta = 0.0_pf
        sin_theta = 1.0_pf
    ELSE
        CALL random_number(X)
        cos_theta = 2.0_pf * X - 1.0_pf
        sin_theta = SQRT(1.0_pf - cos_theta*cos_theta)
    ENDIF

    ! now the cartesian coordinates
    qs(1,:) = r * sin_theta * COS(phi)
    qs(2,:) = r * sin_theta * SIN(phi)
    qs(3,:) = r * cos_theta
END SUBROUTINE

END MODULE