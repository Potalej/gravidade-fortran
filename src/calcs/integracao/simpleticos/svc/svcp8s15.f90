! ************************************************************
!! METODO NUMERICO: svcp8s15
!
! Objetivos:
!   Aplicacao do metodo simpletico svcp8s15.
!   Stormer-Verlet Composto de Ordem 8 e 15 Estagios.
!   Referencia: (Hairer et al, 2006, p.157)
!
! Modificado:
!   12 de julho de 2025
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

  REAL(pf),    DIMENSION(15) :: s, a, b
  REAL(pf128), DIMENSION(15) :: s_128

  REAL(pf),    DIMENSION(14) :: a_consec
  REAL(pf128), DIMENSION(14) :: s_consec

  TYPE, EXTENDS(integracao) :: integracao_svcp8s15

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
  class(integracao_svcp8s15), INTENT(IN) :: self
  INTEGER :: i

  s_128 = (/0.74167036435061295344822780_pf128, &
       -0.40910082580003159399730010_pf128, &
        0.19075471029623837995387626_pf128, &
       -0.57386247111608226665638773_pf128, &
        0.29906418130365592384446354_pf128, &
        0.33462491824529818378495798_pf128, &
        0.31529309239676659663205666_pf128, &
       -0.79688793935291635401978884_pf128, &
        0.31529309239676659663205666_pf128, &
        0.33462491824529818378495798_pf128, &
        0.29906418130365592384446354_pf128, &
       -0.57386247111608226665638773_pf128, &
        0.19075471029623837995387626_pf128, &
       -0.40910082580003159399730010_pf128, &
        0.74167036435061295344822780_pf128 /)

  IF (self % mi) THEN
    a = self % m2_128 * (0.5_pf128 * s_128)
    b = self % m_inv_128 * s_128

    DO i = 1, SIZE(s)-1
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
!   12 de julho de 2025
!
! Autoria:
!   oap
!
FUNCTION metodo (self, R, P, FSomas_ant)

  IMPLICIT NONE
  class(integracao_svcp8s15), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: P_meio, R1, P1
  REAL(pf), DIMENSION(self%N, self%dim) :: FSomas
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo
  INTEGER :: i
  REAL(pf) :: h ! tamanho de passo dinamico

  R1 = R
  P1 = P
  FSomas = FSomas_ant

  DO i=size(s),1,-1
    h = s(i) * self % h

    P_meio = P1 + 0.5_pf * h * FSomas
    R1 = R1 + h * P_meio * self % massasInvertidas
    FSomas = self%forcas(R1)
    P1 = P_meio + 0.5_pf * h * FSomas
  END DO

  metodo(1,:,:) = R1
  metodo(2,:,:) = P1
  metodo(3,:,:) = FSomas

END FUNCTION metodo

! ************************************************************
!! Metodo numerico (massas iguais)
!
! Objetivos:
!   Aplicacao do metodo em si.
!
! Modificado:
!   12 de julho de 2025
!
! Autoria:
!   oap
!
FUNCTION metodo_mi (self, R, P, FSomas_ant)
  IMPLICIT NONE
  class(integracao_svcp8s15), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: P_meio, R1, P1
  REAL(pf), DIMENSION(self%N, self%dim) :: FSomas
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo_mi
  INTEGER :: i

  R1 = R
  P1 = P
  FSomas = FSomas_ant

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
  metodo_mi(3,:,:) = FSomas

END FUNCTION metodo_mi

END MODULE svcp8s15