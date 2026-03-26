! ************************************************************
!! METODO NUMERICO: svcp8s17
!
! Objetivos:
!   Aplicacao do metodo simpletico svcp8s17.
!   Stormer-Verlet Composto de Ordem 8 e 17 Estagios.
!   Referencia: (Hairer et al, 2006, p.157)
!
! Modificado:
!   18 de janeiro de 2026
!
! Autoria:
!   oap
!
MODULE svcp8s17
  USE tipos
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_svcp8s17

  REAL(pf),    DIMENSION(17) :: s, a, b
  REAL(pf128), DIMENSION(17) :: s_128

  REAL(pf),    DIMENSION(16) :: a_consec
  REAL(pf128), DIMENSION(16) :: s_consec

  TYPE, EXTENDS(integracao) :: integracao_svcp8s17

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
  class(integracao_svcp8s17), INTENT(IN) :: self
  INTEGER :: i

  s_128 = (/ 0.13020248308889008087881763_pf128, &
             0.56116298177510838456196441_pf128, &
            -0.38947496264484728640807860_pf128, &
             0.15884190655515560089621075_pf128, &
            -0.39590389413323757733623154_pf128, &
             0.18453964097831570709183254_pf128, &
             0.25837438768632204729397911_pf128, &
             0.29501172360931029887096624_pf128, &
            -0.60550853383003451169892108_pf128, &
             0.29501172360931029887096624_pf128, &
             0.25837438768632204729397911_pf128, &
             0.18453964097831570709183254_pf128, &
            -0.39590389413323757733623154_pf128, &
             0.15884190655515560089621075_pf128, &
            -0.38947496264484728640807860_pf128, &
             0.56116298177510838456196441_pf128, &
             0.13020248308889008087881763_pf128 /)

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
!   18 de janeiro de 2026
!
! Autoria:
!   oap
!
SUBROUTINE metodo (self, R, P, FSomas)

  IMPLICIT NONE
  class(integracao_svcp8s17), INTENT(INOUT) :: self
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
  class(integracao_svcp8s17), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  INTEGER :: i

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

END MODULE svcp8s17