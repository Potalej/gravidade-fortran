! ************************************************************
!! METODO NUMERICO: svcp10s35
!
! Objetivos:
!   Aplicacao do metodo simpletico svcp10s35.
!   Stormer-Verlet Composto de Ordem 10 e 35 Estagios.
!   Referencia: (Hairer et al, 2006, p.158)
!
! Modificado:
!   03 de junho de 2025
!
! Autoria:
!   oap
!
MODULE svcp10s35
  USE tipos
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_svcp10s35

  REAL(pf),    DIMENSION(35) :: s, a, b
  REAL(pf128), DIMENSION(35) :: s_128

  REAL(pf),    DIMENSION(34) :: a_consec
  REAL(pf128), DIMENSION(34) :: s_consec

  TYPE, EXTENDS(integracao) :: integracao_svcp10s35

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
  class(integracao_svcp10s35), INTENT(IN) :: self
  INTEGER :: i

  s_128 = (/ 0.07879572252168641926390768_pf128, &
         0.31309610341510852776481247_pf128, &
         0.02791838323507806610952027_pf128, &
        -0.22959284159390709415121340_pf128, &
         0.13096206107716486317465686_pf128, &
        -0.26973340565451071434460973_pf128, &
         0.07497334315589143566613711_pf128, &
         0.11199342399981020488957508_pf128, &
         0.36613344954622675119314812_pf128, &
        -0.39910563013603589787862981_pf128, &
         0.10308739852747107731580277_pf128, &
         0.41143087395589023782070412_pf128, &
        -0.00486636058313526176219566_pf128, &
        -0.39203335370863990644808194_pf128, &
         0.05194250296244964703718290_pf128, &
         0.05066509075992449633587434_pf128, &
         0.04967437063972987905456880_pf128, &
         0.04931773575959453791768001_pf128, &
         0.04967437063972987905456880_pf128, &
         0.05066509075992449633587434_pf128, &
         0.05194250296244964703718290_pf128, &
        -0.39203335370863990644808194_pf128, &
        -0.00486636058313526176219566_pf128, &
         0.41143087395589023782070412_pf128, &
         0.10308739852747107731580277_pf128, &
        -0.39910563013603589787862981_pf128, &
         0.36613344954622675119314812_pf128, &
         0.11199342399981020488957508_pf128, &
         0.07497334315589143566613711_pf128, &
        -0.26973340565451071434460973_pf128, &
         0.13096206107716486317465686_pf128, &
        -0.22959284159390709415121340_pf128, &
         0.02791838323507806610952027_pf128, &
         0.31309610341510852776481247_pf128, &
         0.07879572252168641926390768_pf128 /)

  IF (self % mi) THEN
    a = self % m2_128 * (0.5_pf128 * s_128)
    b = self % m_inv_128 * s_128

    DO i = 1, SIZE(s_128)-1
      s_consec(i) = s_128(i) + s_128(i+1)
      a_consec(i) = self % m2_128 * (0.5_pf128 * s_consec(i))
    END DO
  ELSE
    s = s_128
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
  class(integracao_svcp10s35), INTENT(IN) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: P_meio, R1, P1, FSomas_prox
  REAL(pf), DIMENSION(self%N, self%dim) :: FSomas
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo
  INTEGER :: i
  REAL(pf) :: h ! tamanho de passo dinamico

  R1 = R
  P1 = P
  FSomas = self%forcas(R1)

  DO i=size(s),1,-1
    h = s(i) * self % h

    P_meio = P1 + 0.5_pf * h * FSomas
    R1 = R1 + h * P_meio * self % massasInvertidas
    FSomas = self%forcas(R1)
    P1 = P_meio + 0.5_pf * h * FSomas
  END DO

  metodo(1,:,:) = R1
  metodo(2,:,:) = P1
  metodo(3,:,:) = 0.0_pf

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
  class(integracao_svcp10s35), INTENT(IN) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: P_meio, R1, P1, FSomas_prox
  REAL(pf), DIMENSION(self%N, self%dim) :: FSomas
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo_mi
  INTEGER :: i
  REAL(pf) :: h ! tamanho de passo dinamico

  R1 = R
  P1 = P
  FSomas = self%m2 * self%forcas(R1)
  
  i = size(s_128)
  P_meio = P1 + self%h * (a(i) * FSomas)
  R1 = R1 + self%h * b(i) * P_meio
  FSomas = self%forcas(R1)
  P1 = P_meio + self%h * (a(i) * FSomas)

  DO i=size(s_128)-1,1,-1
    P_meio = P_meio + self%h * (a_consec(i) * FSomas)
    R1 = R1 + self%h * b(i) * P_meio
    FSomas = self%forcas(R1)
    P1 = P_meio + self%h * (a(i) * FSomas)
  END DO

  metodo_mi(1,:,:) = R1
  metodo_mi(2,:,:) = P1
  metodo_mi(3,:,:) = 0.0_pf

END FUNCTION metodo_mi

END MODULE svcp10s35