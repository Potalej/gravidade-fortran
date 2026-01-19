! ************************************************************
!! METODO NUMERICO: svcp10s35
!
! Objetivos:
!   Aplicacao do metodo simpletico svcp10s35.
!   Stormer-Verlet Composto de Ordem 10 e 35 Estagios.
!   Referencia: (Hairer et al, 2006, p.158)
!
! Modificado:
!   17 de setembro de 2024 (criado)
!   18 de janeiro de 2026 (modificado)
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
!   18 de janeiro de 2026
!
! Autoria:
!   oap
!
SUBROUTINE metodo (self, R, P, FSomas)

  IMPLICIT NONE
  class(integracao_svcp10s35), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  INTEGER :: i
  REAL(pf) :: h ! tamanho de passo dinamico

  DO i=size(s),1,-1
    h = s(i) * self % h

    P = P + 0.5_pf * h * FSomas
    R = R + h * P * self % massasInvertidas
    FSomas = self%forcas(R)
    P = P + 0.5_pf * h * FSomas
  END DO

END SUBROUTINE metodo

! ************************************************************
!! Metodo numerico (massas iguais)
!
! Objetivos:
!   Aplicacao do metodo em si.
!
! Modificado:
!   18 de janeiro de 2026
!
! Autoria:
!   oap
!
SUBROUTINE metodo_mi (self, R, P, FSomas)
  IMPLICIT NONE
  class(integracao_svcp10s35), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  INTEGER :: i

  ! i = size(s_128)
  ! P = P + self%h * (a(i) * FSomas)
  ! R = R + self%h * b(i) * P
  ! FSomas = self%forcas(R)
  ! P = P + self%h * (a(i) * FSomas)

  ! DO i=size(s_128)-1,1,-1
  !   P = P + self%h * (a_consec(i) * FSomas)
  !   R = R + self%h * b(i) * P
  !   FSomas = self%forcas(R)
  !   P = P + self%h * (a(i) * FSomas)
  ! END DO

  i = size(s_128)
  P = P + self%h * (a(i) * FSomas)
  R = R + self%h * b(i) * P
  FSomas = self%forcas(R)

  DO i=size(s_128)-1,1,-1
    P = P + self%h * (a_consec(i) * FSomas)
    R = R + self%h * b(i) * P
    FSomas = self%forcas(R)
  END DO

  P = P + self%h * (a(size(s_128)) * FSomas)

END SUBROUTINE metodo_mi

END MODULE svcp10s35