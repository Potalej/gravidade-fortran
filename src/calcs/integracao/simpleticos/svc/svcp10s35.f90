! ************************************************************
!! METODO NUMERICO: svcp10s35
!
! Objetivos:
!   Aplicacao do metodo simpletico svcp10s35.
!   Stormer-Verlet Composto de Ordem 10 e 35 Estagios.
!   Referencia: (Hairer et al, 2006, p.158)
!
! Modificado:
!   17 de setembro de 2024
!
! Autoria:
!   oap
!
MODULE svcp10s35
  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_svcp10s35

  REAL(pf), DIMENSION(35) :: s = (/  &
         0.07879572252168641926390768_pf, &
         0.31309610341510852776481247_pf, &
         0.02791838323507806610952027_pf, &
        -0.22959284159390709415121340_pf, &
         0.13096206107716486317465686_pf, &
        -0.26973340565451071434460973_pf, &
         0.07497334315589143566613711_pf, &
         0.11199342399981020488957508_pf, &
         0.36613344954622675119314812_pf, &
        -0.39910563013603589787862981_pf, &
         0.10308739852747107731580277_pf, &
         0.41143087395589023782070412_pf, &
        -0.00486636058313526176219566_pf, &
        -0.39203335370863990644808194_pf, &
         0.05194250296244964703718290_pf, &
         0.05066509075992449633587434_pf, &
         0.04967437063972987905456880_pf, &
         0.04931773575959453791768001_pf, &
         0.04967437063972987905456880_pf, &
         0.05066509075992449633587434_pf, &
         0.05194250296244964703718290_pf, &
        -0.39203335370863990644808194_pf, &
        -0.00486636058313526176219566_pf, &
         0.41143087395589023782070412_pf, &
         0.10308739852747107731580277_pf, &
        -0.39910563013603589787862981_pf, &
         0.36613344954622675119314812_pf, &
         0.11199342399981020488957508_pf, &
         0.07497334315589143566613711_pf, &
        -0.26973340565451071434460973_pf, &
         0.13096206107716486317465686_pf, &
        -0.22959284159390709415121340_pf, &
         0.02791838323507806610952027_pf, &
         0.31309610341510852776481247_pf, &
         0.07879572252168641926390768_pf /)


  TYPE, EXTENDS(integracao) :: integracao_svcp10s35

    CONTAINS
        PROCEDURE :: metodo

  END TYPE

CONTAINS

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
  REAL(pf), DIMENSION(self%N, self%dim) :: R_meio, R1, P1, FSomas_prox
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo
  INTEGER :: i
  REAL(pf) :: h ! tamanho de passo dinamico

  R1 = R
  P1 = P

  DO i=size(s),1,-1
    h = s(i) * self % h

    R_meio = R1 + 0.5_pf * h * P1 * self % massasInvertidas
    P1 = P1 + h * self%forcas(R_meio)
    R1 = R_meio + 0.5_pf * h * P1 * self % massasInvertidas
  END DO

  metodo(1,:,:) = R1
  metodo(2,:,:) = P1
  metodo(3,:,:) = 0.0_pf

END FUNCTION metodo

END MODULE svcp10s35