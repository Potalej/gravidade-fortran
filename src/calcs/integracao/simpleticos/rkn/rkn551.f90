! ************************************************************
!! METODO NUMERICO: RKN55 (1)
!
! Objetivos:
!   Aplicacao do metodo simpletico RKN55 (1).
!   Runge-Kutta-Nystrom de Ordem 5 e 5 Estagios. Eh o primeiro
!   metodo na tabela 1.
!   Referencia: (Okunbor & Skeel, 1992, p.378)
!
! Modificado:
!   03 de junho de 2025
!
! Autoria:
!   oap
!
MODULE rkn551
  USE tipos
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_rkn551

  REAL(pf), DIMENSION(5)   :: g, c, b
  REAL(pf), DIMENSION(4,4) :: aij

  REAL(pf128), DIMENSION(5)   :: g_128, c_128, b_128
  REAL(pf128), DIMENSION(4,4) :: aij_128

  TYPE, EXTENDS(integracao) :: integracao_rkn551

    CONTAINS
      PROCEDURE :: metodo, metodo_mi, atualizar_constantes

  END TYPE

CONTAINS

! ************************************************************
!! Atualizacao das constantes
!
! Modificado:
!   03 de junho de 2025
!
! Autoria:
!   oap
!
SUBROUTINE atualizar_constantes (self)
  IMPLICIT NONE
  class(integracao_rkn551), INTENT(IN) :: self

  ! Atualiza a constante c
  c_128 = (/ 0.69491389107017931259_pf128, 0.63707199676998338411_pf128, &
            -0.02055756998211598005_pf128, 0.79586189634575355001_pf128, &
             0.30116624272377778837_pf128 /)

  ! Atualiza a constante aij
  aij_128(1, 1) =  9.66427576646961680621728658676147461E-0002_pf128
  aij_128(2, 1) =  1.19541605539808219588717185072600831_pf128
  aij_128(2, 2) = -0.803254424095807837234189875982701731_pf128
  aij_128(3, 1) = -0.168664878862877104515577328205108636_pf128      
  aij_128(3, 2) =  0.193952231299022844326977914571762075_pf128      
  aij_128(3, 3) =  7.22491700286848878265806124545633774E-0002_pf128
  aij_128(4, 1) =  0.657877049937500205642537677288055438_pf128      
  aij_128(4, 2) = -0.410288390837650406220850950479507422_pf128      
  aij_128(4, 3) =  2.84709990640056776782580302096903318E-0002_pf128
  aij_128(4, 4) = -0.474893431444463157194499740600585919_pf128

  ! Atualiza a constante b
  b_128(1) = -0.509740638168924449099219894409179678_pf128
  b_128(2) =  0.443294487048759456108835303783416735_pf128
  b_128(3) =  9.03144034820916384258222416974604182E-0002_pf128
  b_128(4) =  0.195966631667991276791308071613311759_pf128
  b_128(5) =  0.280165106288905703470427826046943676_pf128

  ! Atualiza a constante g
  g_128 = (/ -1.67080892327314312060_pf128, 1.22143909230997538270_pf128, &
              0.08849515813253908125_pf128, 0.95997088013770159876_pf128, &
              0.40090379269297793385_pf128 /)
  
  ! Se for massas iguais
  IF (self % mi) THEN
    c = c_128 * self % m_inv_128
    aij = aij_128 * self % m_esc_128
    g = g_128 * self % m2_128

  ! Se nao
  ELSE
    c = c_128
    aij = aij_128
    g = g_128
  ENDIF
END SUBROUTINE atualizar_constantes

! ************************************************************
!! Metodo numerico
!
! Objetivos:
!   Aplicacao do metodo em si.
!
! Modificado:
!   17 de setembro de 2024
!
! Autoria:
!   oap
!
FUNCTION metodo (self, R, P, FSomas_ant)

  IMPLICIT NONE
  class(integracao_rkn551), INTENT(IN) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, FSomas_prox
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo
  
  INTEGER :: i, j

  REAL(pf), DIMENSION(SIZE(c), self%N, self%dim) :: y, fi

  R1 = R + self % h * P * self % massasInvertidas
  P1 = P

  DO i = 1, SIZE(c)

    ! Calcula a base do yi
    y(i,:,:) = R + c(i) * self % h * P * self % massasInvertidas

    ! Se i > 1, calcula os termos a mais
    IF (i .GT. 1) THEN
      ! Somatorio
      DO j = 1, i-1
        y(i,:,:) = y(i,:,:) + self % h * self % h * aij(i-1,j) * fi(j,:,:)
      END DO
    ENDIF

    ! Calcula fi
    FSomas_prox = self%forcas(y(i,:,:))
    fi(i,:,:) = FSomas_prox * self % massasInvertidas

    ! Adiciona nas posicoes e momentos
    R1 = R1 + self % h * self % h * b(i) * fi(i,:,:)
    P1 = P1 + self % h * g(i) * FSomas_prox
  END DO

  metodo(1,:,:) = R1
  metodo(2,:,:) = P1

END FUNCTION metodo

! ************************************************************
!! Metodo numerico (massas iguais)
!
! Objetivos:
!   Aplicacao do metodo em si.
!
! Modificado:
!   03 de junho de 2025
!
! Autoria:
!   oap
!
FUNCTION metodo_mi (self, R, P, FSomas_ant)
  IMPLICIT NONE
  class(integracao_rkn551), INTENT(IN) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, FSomas_prox
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo_mi
  INTEGER :: i, j
  REAL(pf), DIMENSION(SIZE(c), self%N, self%dim) :: fi
  REAL(pf), DIMENSION(self%N, self%dim) :: yi

  R1 = R + self % h * (self % m_inv * P)
  P1 = P

  DO i = 1, SIZE(c)

    ! Calcula a base do yi
    yi = R + self%h * (c(i) * P)

    ! Se i > 1, calcula os termos a mais
    IF (i .GT. 1) THEN
      ! Somatorio
      DO j = 1, i-1
        yi = yi + (self%h * self% h) * (aij(i-1,j) * fi(j,:,:))
      END DO
    ENDIF

    ! Calcula fi
    FSomas_prox = self%forcas(yi)
    fi(i,:,:) = FSomas_prox

    ! Adiciona nas posicoes e momentos
    R1 = R1 + (self%h * self% h) * (b(i) * FSomas_prox)
    P1 = P1 + self%h * (g(i) * FSomas_prox)
  END DO

  metodo_mi(1,:,:) = R1
  metodo_mi(2,:,:) = P1

END FUNCTION metodo_mi

END MODULE rkn551