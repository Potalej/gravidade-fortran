! *****************************************************************
!! MECANICA
!
! Objetivos:
!   Este arquivo contem helpers de mecanica. O modulo contem
!   funcionalidades como o calculo de momento angular, energia,
!   momento de inercia, etc.
!
!   O modulo "auxiliares" eh utilizado em algumas funcoes/rotinas,
!   como no calculo do momento angular, que usa o produto vetorial.
! 
! Modificado:
!   02 de fevereiro de 2024
! 
! Autoria:
!   oap
! 
MODULE mecanica
  USE tipos
  USE auxiliares
  IMPLICIT NONE

  INTERFACE energia_cinetica
    MODULE PROCEDURE energia_cinetica_vec
    MODULE PROCEDURE energia_cinetica_esc
  END INTERFACE

  INTERFACE energia_potencial
    MODULE PROCEDURE energia_potencial_vec
    MODULE PROCEDURE energia_potencial_esc
  END INTERFACE
  
  INTERFACE energia_total
    MODULE PROCEDURE energia_total_vec
    MODULE PROCEDURE energia_total_esc
  END INTERFACE
  
  INTERFACE momento_inercia
    MODULE PROCEDURE momento_inercia_vec
    MODULE PROCEDURE momento_inercia_esc
  END INTERFACE
CONTAINS

! ************************************************************
!! Momento angular individual
!
! Objetivos:
!   Calcula o momento angular de um corpo dado sua posicao e
!   seu momento linear.
!
! Modificado:
!   03 de fevereiro de 2024
!
! Autoria:
!   oap
! 
FUNCTION momento_angular_individual (R, P)
  IMPLICIT NONE
  REAL(pf), DIMENSION(3), INTENT(IN) :: R, P
  REAL(pf), DIMENSION(3)             :: momento_angular_individual
  momento_angular_individual = produto_vetorial(R,P)
END FUNCTION momento_angular_individual

! ************************************************************
!! Momento angular total
!
! Objetivos:
!   Calcula o momento angular de uma lista de corpos dadas as
!   suas posicoes e seus momentos lineares.
!
! Modificado:
!   03 de fevereiro de 2024
!
! Autoria:
!   oap
!
FUNCTION momento_angular_total (Rs, Ps)
  IMPLICIT NONE
  REAL(pf), DIMENSION(:,:), INTENT(IN) :: Rs, Ps
  REAL(pf), DIMENSION(3)               :: momento_angular_total
  INTEGER                              :: i
  momento_angular_total = (/0.0_pf,0.0_pf,0.0_pf/)
  DO i=1, SIZE(Rs,1)
    momento_angular_total = momento_angular_total &
                          + momento_angular_individual(Rs(i,:),Ps(i,:))
  END DO
END FUNCTION momento_angular_total

! ************************************************************
!! Energia cinetica
!
! Objetivos:
!   Calcula a energia cinetica de um conjunto de corpos dadas
!   as suas massas e momentos linares.
!
! Modificado:
!   03 de fevereiro de 2024
!
! Autoria:
!   oap
!
FUNCTION energia_cinetica_vec (m, P) RESULT(ec)
  IMPLICIT NONE
  REAL(pf) :: m(:), P(:,:), ec
  INTEGER  :: i
  ec=0.0_pf
  DO i=1, SIZE(m)
    ec = ec + DOT_PRODUCT(P(i,:), P(i,:))/m(i)
  END DO
  ec = 0.5_pf * ec
END FUNCTION energia_cinetica_vec

FUNCTION energia_cinetica_esc (m, P) RESULT(ec)
  IMPLICIT NONE
  REAL(pf) :: m ! nesse caso, m = m_inv
  REAL(pf) :: P(:,:), ec
  INTEGER  :: i
  ec = 0.5_pf * sum(sum(P**2, dim=2)) * m
END FUNCTION energia_cinetica_esc

! ************************************************************
!! Energia potencial
!
! Objetivos:
!   Calcula a energia potencial newtoniana dada a constante
!   de gravitacao universal, as massas e as posicoes.
!
!   A formula utilizada eh a seguinte:
!   V = - G sum (i < j) m_i m_j / r_ij
!
! Modificado:
!   03 de fevereiro de 2024
!
! Autoria:
!   oap
!
FUNCTION energia_potencial_vec (G,m,R) RESULT(ep)
  IMPLICIT NONE
  REAL(pf) :: m(:), R(:,:)
  REAL(pf) :: distancia, G, ep, distancia_inv
  INTEGER  :: i,j
  ep=0.0_pf
  DO i=2, SIZE(m)
    DO j=1,i-1
      distancia = NORM2(R(i,:)-R(j,:))
      distancia_inv = 1.0_pf/distancia
      ep = ep + m(i)*m(j)*distancia_inv
    END DO
  END DO
  ep = -G*ep
END FUNCTION energia_potencial_vec

FUNCTION energia_potencial_esc (G,m,R) RESULT(ep)
  IMPLICIT NONE
  REAL(pf) :: m ! aqui m = m*m
  REAL(pf) :: R(:,:)
  REAL(pf) :: distancia, G, ep, distancia_inv
  INTEGER  :: i,j
  ep=0.0_pf
  DO i=2, SIZE(R,1)
    DO j=1,i-1
      distancia = NORM2(R(i,:)-R(j,:))
      distancia_inv = 1.0_pf/distancia
      ep = ep + distancia_inv
    END DO
  END DO
  ep = -G*ep*m
END FUNCTION energia_potencial_esc

! ************************************************************
!! Energia total do problema de N corpos
!
! Objetivos:
!   Calcula a energia total do problema de N corpos dados o
!   o valor da constante de gravitacao universal, as massas,
!   as posicoes e os momentos lineares.
!
!   A energia total eh dada pela soma da energia cinetica com
!   a energia potencial newtoniana
!
! Modificado:
!   03 de fevereiro de 2024
!
! Autoria:
!   oap
!
FUNCTION energia_total_vec (G, m, R, P) RESULT(e_tot)
  IMPLICIT NONE
  REAL(pf) :: m(:), R(:,:), P(:,:), G
  REAL(pf) :: e_tot
  e_tot = energia_cinetica(m,P) + energia_potencial(G,m,R)
END FUNCTION energia_total_vec

FUNCTION energia_total_esc (G, m, R, P) RESULT(e_tot)
  IMPLICIT NONE
  REAL(pf) :: m
  REAL(pf) :: R(:,:), P(:,:), G
  REAL(pf) :: e_tot
  e_tot = energia_cinetica(1/m,P) + energia_potencial(G,m*m,R)
END FUNCTION energia_total_esc

! ************************************************************
!! Momento de dilatacao
!
! Objetivos:
!   Calcula o momento de dilatacao.
!
! Modificado:
!   05 de janeiro de 2025
!
! Autoria:
!   oap
!
FUNCTION momento_dilatacao (R, P)
  IMPLICIT NONE
  REAL(pf) :: R(:,:), P(:,:)
  REAL(pf) :: momento_dilatacao
  INTEGER  :: i
  momento_dilatacao=0.0_pf
  DO i=1, SIZE(R,1)
    momento_dilatacao = momento_dilatacao + DOT_PRODUCT(R(i,:),P(i,:))
  END DO
END FUNCTION momento_dilatacao

! ************************************************************
!! Momento de inercia
!
! Objetivos:
!   Calcula o momento de inercia.
!
! Modificado:
!   07 de janeiro de 2025
!
! Autoria:
!   oap
!
FUNCTION momento_inercia_vec (m, R) RESULT(momine)
  IMPLICIT NONE
  REAL(pf) :: m(:), R(:,:)
  REAL(pf) :: momine
  INTEGER  :: i
  momine=0.0_pf
  DO i=1, SIZE(R,1)
    momine = momine + m(i) * DOT_PRODUCT(R(i,:),R(i,:))
  END DO
END FUNCTION momento_inercia_vec

FUNCTION momento_inercia_esc (m, R) RESULT(momine)
  IMPLICIT NONE
  REAL(pf) :: m, R(:,:)
  REAL(pf) :: momine
  INTEGER  :: i
  momine=m*sum(R(:,1)**2 + R(:,2)**2 + R(:,3)**2)
END FUNCTION momento_inercia_esc

! ************************************************************
!! Raio de meia massa (half-mass ratio)
!
! Objetivos:
!   Calcula o raio de meia massa.
!
! Modificado:
!   30 de abril de 2025
!
! Autoria:
!   oap
!
FUNCTION raio_meia_massa (m, R)
  IMPLICIT NONE
  REAL(pf) :: m(:), R(:,:)
  REAL(pf) :: raio_meia_massa, M_met, M_soma, qcm(3)
  INTEGER :: i, j
  INTEGER, ALLOCATABLE :: idx_ord(:)
  REAL(pf), ALLOCATABLE :: raios(:)

  ALLOCATE(idx_ord(SIZE(m)))
  ALLOCATE(raios(SIZE(m)))
  raios = 0.0_pf
  
  ! Calcula os raios
  qcm = centro_massas(m, R)
  DO i=1, SIZE(m)
    raios(i) = NORM2(R(i,:) - qcm)
  END DO

  ! Ordena a lista em ordem crescente
  idx_ord = ordenar_lista_crescente(raios)

  ! Agora vai pegando as massas ate chegar na metade
  M_met = 0.5_pf * SUM(m)
  M_soma = 0.0_pf
  DO i=1, SIZE(m)
    j = idx_ord(i)
    M_soma = M_soma + m(j)
    IF (M_soma .GE. M_met) THEN
      raio_meia_massa = raios(j)
      EXIT
    END IF
  END DO

END FUNCTION raio_meia_massa

! ************************************************************
!! Tempo de relaxacao do raio de meia massa t_rh
!
! Objetivos:
!   Calcula o tempo de relaxacao do raio de meia massa t_rh
!
! Modificado:
!   30 de abril de 2025
!
! Autoria:
!   oap
!
FUNCTION tempo_relaxacao_rh (m, R)
  IMPLICIT NONE
  REAL(pf) :: G=1.0_pf, m(:), R(:,:), tempo_relaxacao_rh
  REAL(pf) :: rh, m_media, log_coulomb
  INTEGER :: N
  
  N = SIZE(m)
  ! Calcula o raio de meia massa
  rh = raio_meia_massa(m, R)
  WRITE(*,*) 'rh:', rh

  ! Media das massas
  m_media = SUM(m) / N

  ! Logaritmo de coulomb
  log_coulomb = LOG(0.4_pf * N)

  ! Aproximacao do valor de t_rh
  tempo_relaxacao_rh = (0.138_pf * N / log_coulomb) * SQRT(rh*rh*rh / G*m_media)

END FUNCTION tempo_relaxacao_rh

END MODULE mecanica