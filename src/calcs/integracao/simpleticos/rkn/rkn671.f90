! ************************************************************
!! METODO NUMERICO: RKN67 (1)
!
! Objetivos:
!   Aplicacao do metodo simpletico RKN67 (1).
!   Runge-Kutta-Nystrom de Ordem 6 e 7 Estagios. Eh o primeiro
!   metodo na tabela 2.
!   Referencia: (Okunbor & Skeel, 1992, p.380)
!
! Modificado:
!   26 de marco de 2026
!
! Autoria:
!   oap
!  
MODULE rkn671
  USE tipos
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_rkn671

  REAL(pf), DIMENSION(7)   :: g, c, b
  REAL(pf), DIMENSION(6,6) :: aij

  TYPE, EXTENDS(integracao) :: integracao_rkn671

    CONTAINS
      PROCEDURE :: metodo, metodo_mi, atualizar_constantes

  END TYPE

CONTAINS

SUBROUTINE atualizar_constantes (self)
  IMPLICIT NONE
  class(integracao_rkn671), INTENT(IN) :: self
  
  ! Atualiza a constante g
  g = (/ -0.68774007118557290171_pf, 0.13118241020105280626_pf, &
          0.92161977504885189358_pf, 0.26987577187133640373_pf, &
          0.92161977504885189358_pf, 0.13118241020105280626_pf, &
          -0.68774007118557290171_pf /)

  ! Atualiza a constante c
  c = (/  0.94413392188212619249_pf, 0.34626230516255218639_pf, &
          0.93479137012319657440_pf, 0.50000000000000000000_pf, &
          0.06520862987680341172_pf, 0.65373769483744781361_pf, &
          0.05586607811787376587_pf /)

  ! Atualiza a constante b
  b = (/ -0.03842134054164531021_pf, 0.08575888644805676475_pf, &
          0.06009756279830340969_pf, 0.13493788593566818923_pf, &
          0.86152221225054848031_pf, 0.04542352375299604783_pf, &
          -0.64931873064392764405_pf /)

  ! Atualiza a constante aij
  aij = reshape( &
  (/ 0.41118026824255338170_pf,  0.0_pf,                     0.0_pf,                     0.0_pf, 0.0_pf, 0.0_pf, &
     0.00642524721174115524_pf,  0.07720466121490930644_pf,  0.0_pf,                     0.0_pf, 0.0_pf, 0.0_pf, &
     0.30544869505114113917_pf,  0.02016768134753035846_pf, -0.40071232472612250408_pf,  0.0_pf, 0.0_pf, 0.0_pf, &
     0.60447214289054107539_pf, -0.03686929851984859646_pf, -0.80142464945224500816_pf, -0.11733965661499358435_pf, &
     0.0_pf,                     0.0_pf,                     0.19971712185972889664_pf,  0.04033536269506071692_pf, &
    -0.25902462499350481506_pf,  0.04149007905997619677_pf,  0.54240002445874024861_pf,  0.0_pf, &
     0.61089739010228227833_pf, -0.03809487697701307435_pf, -0.81003492990269199137_pf, -0.11986098498218263064_pf, &
     -0.00861028045044699535_pf, -0.07843023967207378433_pf /), &
     shape(aij), order=(/2,1/))

  ! Se for massas iguais
  IF (self % mi) THEN
    c = c * self % m_inv_128
    aij = aij * self % m_esc_128
    g = g * self % m2_128
    b = b * self % m_esc_128
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
SUBROUTINE metodo (self, R, P, FSomas)

  IMPLICIT NONE
  class(integracao_rkn671), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1
  REAL(pf), DIMENSION(self%N, self%dim) :: P_prod
  
  INTEGER :: i, j

  REAL(pf), DIMENSION(SIZE(c), self%N, self%dim) :: y, fi

  P_prod = self % h * P * self % massasInvertidas
  R1 = R + P_prod
  P1 = P

  DO i = 1, SIZE(c)

    ! Calcula a base do yi
    y(i,:,:) = R + c(i) * P_prod

    ! Se i > 1, calcula os termos a mais
    IF (i .GT. 1) THEN
      ! Somatorio
      DO j = 1, i-1
        y(i,:,:) = y(i,:,:) + self % h * self % h * aij(i-1,j) * fi(j,:,:)
      END DO
    ENDIF

    ! Calcula fi
    FSomas = self%forcas(y(i,:,:))
    fi(i,:,:) = FSomas * self % massasInvertidas

    ! Adiciona nas posicoes e momentos
    R1 = R1 + self % h * self % h * b(i) * fi(i,:,:)
    P1 = P1 + self % h * g(i) * FSomas
  END DO

  R = R1
  P = P1

END SUBROUTINE metodo

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
SUBROUTINE metodo_mi (self, R, P, FSomas)
  IMPLICIT NONE
  class(integracao_rkn671), INTENT(INOUT) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(INOUT) :: R, P, FSomas
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1
  REAL(pf), DIMENSION(SIZE(c), self%N, self%dim) :: fi
  REAL(pf), DIMENSION(self%N, self%dim) :: yi
  INTEGER :: i, j

  R1 = R + self % h * (self % m_inv * P)
  P1 = P

  DO i = 1, SIZE(c)

    ! Calcula a base do yi
    yi = R + self%h * (c(i) * P)

    ! Se i > 1, calcula os termos a mais
    IF (i .GT. 1) THEN
      ! Somatorio
      DO j = 1, i-1
        yi = yi + (self%h * self%h) * (aij(i-1,j) * fi(j,:,:))
      END DO
    ENDIF

    ! Calcula fi
    FSomas = self%forcas(yi)
    fi(i,:,:) = FSomas

    ! Adiciona nas posicoes e momentos
    R1 = R1 + (self%h * self%h) * (b(i) * FSomas)
    P1 = P1 + self % h * (g(i) * FSomas)
  END DO

  R = R1
  P = P1

END SUBROUTINE metodo_mi

END MODULE rkn671