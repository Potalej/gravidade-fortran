! ************************************************************
!! METODO NUMERICO: ecp4s6 (Blanes-Moan ordem 4)
!
! Objetivos:
!   Aplicacao do metodo de composicao de Blanes & Moan (2002), 
!   com ordem 4 e 6 estagios, usando como base o Euler 
!   Simpletico.
!   Referencia: (Hairer et al, 2006, p.153)
!
! Modificado:
!   04 de marco de 2026
!
! Autoria:
!   oap
!
MODULE ecp4s6
  USE tipos
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_ecp4s6

  REAL(pf) :: coef_t(7), coef_v(6)
  
  TYPE, EXTENDS(integracao) :: integracao_ecp4s6

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
  class(integracao_ecp4s6), INTENT(IN) :: self
  INTEGER :: i
  REAL(pf) :: a(6), b(6)

  b(1) = 0.082984406417405_pf
  b(2) = 0.23399525073150_pf
  b(3) = -0.40993371990193_pf
  b(4) = 0.059762097006575_pf
  b(5) = 0.37087741497958_pf
  b(6) = 0.16231455076687_pf

  a(1) = b(6)
  a(2) = b(5)
  a(3) = b(4)
  a(4) = b(3)
  a(5) = b(2)
  a(6) = b(1)

  coef_t(1) = b(1)
  coef_v(1) = a(1) + b(1)
 
  coef_t(2) = b(2) + a(1)
  coef_v(2) = a(2) + b(2)

  coef_t(3) = b(3) + a(2)
  coef_v(3) = a(3) + b(3)
 
  coef_t(4) = b(4)+ a(3)
  coef_v(4) = a(4)+ b(4)
  
  coef_t(5) = b(5) + a(4)
  coef_v(5) = a(5) + b(5)

  coef_t(6) = b(6) + a(5)
  coef_v(6) = a(6) + b(6)

  coef_t(7) = a(6)
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
  class(integracao_ecp4s6), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas

  R = R + coef_t(1) * self % h * P * self % massasInvertidas
  FSomas = self%forcas(R)
  P = P + coef_v(1) * self % h * FSomas

  R = R + coef_t(2) * self % h * P * self % massasInvertidas
  FSomas = self%forcas(R)
  P = P + coef_v(2) * self % h * FSomas

  R = R + coef_t(3) * self % h * P * self % massasInvertidas
  FSomas = self%forcas(R)
  P = P + coef_v(3) * self % h * FSomas

  R = R + coef_t(4) * self % h * P * self % massasInvertidas
  FSomas = self%forcas(R)
  P = P + coef_v(4) * self % h * FSomas

  R = R + coef_t(5) * self % h * P * self % massasInvertidas
  FSomas = self%forcas(R)
  P = P + coef_v(5) * self % h * FSomas

  R = R + coef_t(6) * self % h * P * self % massasInvertidas
  FSomas = self%forcas(R)
  P = P + coef_v(6) * self % h * FSomas

  R = R + coef_t(7) * self % h * P * self % massasInvertidas
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
  class(integracao_ecp4s6), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas

  R = R + coef_t(1) * self % h * self % m_inv * P
  FSomas = self%forcas(R)
  P = P + coef_v(1) * self % h * (self % m2 * FSomas)

  R = R + coef_t(2) * self % h * self % m_inv * P
  FSomas = self%forcas(R)
  P = P + coef_v(2) * self % h * (self % m2 * FSomas)

  R = R + coef_t(3) * self % h * self % m_inv * P
  FSomas = self%forcas(R)
  P = P + coef_v(3) * self % h * (self % m2 * FSomas)

  R = R + coef_t(4) * self % h * self % m_inv * P
  FSomas = self%forcas(R)
  P = P + coef_v(4) * self % h * (self % m2 * FSomas)

  R = R + coef_t(5) * self % h * self % m_inv * P
  FSomas = self%forcas(R)
  P = P + coef_v(5) * self % h * (self % m2 * FSomas)

  R = R + coef_t(6) * self % h * self % m_inv * P
  FSomas = self%forcas(R)
  P = P + coef_v(6) * self % h * (self % m2 * FSomas)

  R = R + coef_t(7) * self % h * self % m_inv * P
END SUBROUTINE metodo_mi

END MODULE ecp4s6