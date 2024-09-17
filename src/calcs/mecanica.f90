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
  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  USE auxiliares
  IMPLICIT NONE
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
FUNCTION energia_cinetica (m, P)
  IMPLICIT NONE
  REAL(pf) :: m(:), P(:,:), energia_cinetica
  INTEGER  :: i
  energia_cinetica=0.0_pf
  DO i=1, SIZE(m)
    energia_cinetica = energia_cinetica + NORM2(P(i,:))**2/m(i)
  END DO
  energia_cinetica = 0.5_pf * energia_cinetica
END FUNCTION energia_cinetica

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
FUNCTION energia_potencial (G,m,R)
  IMPLICIT NONE
  REAL(pf) :: m(:), R(:,:)
  REAL(pf) :: distancia, G, energia_potencial, distancia_inv
  INTEGER  :: i,j
  energia_potencial=0.0_pf
  DO i=2, SIZE(m)
    DO j=1,i-1
      distancia = NORM2(R(i,:)-R(j,:))
      distancia_inv = 1.0_pf/distancia
      energia_potencial = energia_potencial + m(i)*m(j)*distancia_inv
    END DO
  END DO
  energia_potencial = -G*energia_potencial
END FUNCTION energia_potencial

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
FUNCTION energia_total (G, m, R, P)
  IMPLICIT NONE
  REAL(pf) :: m(:), R(:,:), P(:,:), G
  REAL(pf) :: energia_total
  energia_total = energia_cinetica(m,P) + energia_potencial(G,m,R)
END FUNCTION energia_total

END MODULE mecanica