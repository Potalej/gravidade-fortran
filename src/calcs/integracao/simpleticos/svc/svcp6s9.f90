! ************************************************************
!! METODO NUMERICO: svcp6s9 (Kahan-Li)
!
! Objetivos:
!   Aplicacao do metodo de composicao de Kahan e Li (1997), com
!   ordem 6 e 9 estagios. Metodo base: Velocity-Verlet.
!   Referencia: (Hairer et al, 2006, p.157)
!
! Modificado:
!   04 de marco de 2026
!
! Autoria:
!   oap
!
MODULE svcp6s9
  USE tipos
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_svcp6s9

  REAL(pf),    DIMENSION(9) :: s, a, b
  REAL(pf128), DIMENSION(9) :: s_128

  REAL(pf),    DIMENSION(9) :: a_consec
  REAL(pf128), DIMENSION(9) :: s_consec

  TYPE, EXTENDS(integracao) :: integracao_svcp6s9

    CONTAINS
      PROCEDURE :: metodo, metodo_mi, atualizar_constantes

  END TYPE

CONTAINS

! ************************************************************
!! Atualizacao das constantes
!
! Modificado:
!   04 de marco de 2026
!
! Autoria:
!   oap
!
SUBROUTINE atualizar_constantes (self)
  IMPLICIT NONE
  class(integracao_svcp6s9), INTENT(IN) :: self
  INTEGER :: i

  s_128 = (/0.39216144400731413927925056_pf128, &
        0.33259913678935943859974864_pf128, &
       -0.70624617255763935980996482_pf128, &
        0.08221359629355080023149045_pf128, &
        
        0.79854399093482996339895035_pf128, &

        0.08221359629355080023149045_pf128, &
       -0.70624617255763935980996482_pf128, &
        0.33259913678935943859974864_pf128, &
        0.39216144400731413927925056_pf128 /)

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
!   04 de marco de 2026
!
! Autoria:
!   oap
!
SUBROUTINE metodo (self, R, P, FSomas)

  IMPLICIT NONE
  class(integracao_svcp6s9), INTENT(INOUT) :: self
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
!   04 de marco de 2026
!
! Autoria:
!   oap
!
SUBROUTINE metodo_mi (self, R, P, FSomas)
  IMPLICIT NONE
  class(integracao_svcp6s9), INTENT(INOUT) :: self
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

END MODULE svcp6s9