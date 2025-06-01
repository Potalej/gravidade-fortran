! ************************************************************
!! METODO NUMERICO: svcp8s15
!
! Objetivos:
!   Aplicacao do metodo simpletico svcp8s15.
!   Stormer-Verlet Composto de Ordem 8 e 15 Estagios.
!   Referencia: (Hairer et al, 2006, p.157)
!
! Modificado:
!   17 de setembro de 2024
!
! Autoria:
!   oap
!
MODULE svcp8s15
  USE tipos
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_svcp8s15

  REAL(pf), DIMENSION(15) :: s = (/  &
        0.74167036435061295344822780_pf, &
       -0.40910082580003159399730010_pf, &
        0.19075471029623837995387626_pf, &
       -0.57386247111608226665638773_pf, &
        0.29906418130365592384446354_pf, &
        0.33462491824529818378495798_pf, &
        0.31529309239676659663205666_pf, &
       -0.79688793935291635401978884_pf, &
        0.31529309239676659663205666_pf, &
        0.33462491824529818378495798_pf, &
        0.29906418130365592384446354_pf, &
       -0.57386247111608226665638773_pf, &
        0.19075471029623837995387626_pf, &
       -0.40910082580003159399730010_pf, &
        0.74167036435061295344822780_pf /)

  TYPE, EXTENDS(integracao) :: integracao_svcp8s15

    CONTAINS
      PROCEDURE :: metodo, metodo_mi

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
  class(integracao_svcp8s15), INTENT(IN) :: self
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
!   01 de junho de 2025
!
! Autoria:
!   oap
!
FUNCTION metodo_mi (self, R, P, FSomas_ant)
  IMPLICIT NONE
  class(integracao_svcp8s15), INTENT(IN) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: P_meio, R1, P1, FSomas_prox
  REAL(pf), DIMENSION(self%N, self%dim) :: FSomas
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo_mi
  INTEGER :: i
  REAL(pf) :: h ! tamanho de passo dinamico

  R1 = R
  P1 = P
  FSomas = self%m2 * self%forcas(R1)

  DO i=size(s),1,-1
    h = s(i) * self % h

    P_meio = P1 + 0.5_pf * h * FSomas
    R1 = R1 + h * P_meio * self % m_inv
    FSomas = self%m2 * self%forcas(R1)
    P1 = P_meio + 0.5_pf * h * FSomas
  END DO

  metodo_mi(1,:,:) = R1
  metodo_mi(2,:,:) = P1
  metodo_mi(3,:,:) = 0.0_pf

END FUNCTION metodo_mi

END MODULE svcp8s15