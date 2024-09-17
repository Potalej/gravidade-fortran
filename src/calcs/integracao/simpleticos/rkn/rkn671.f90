! ************************************************************
!! METODO NUMERICO: RKN67 (1)
!
! Objetivos:
!   Aplicacao do metodo simpletico RKN67 (1).
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
MODULE rkn671
  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_rkn671

  REAL(pf), DIMENSION(7) :: g = (/ -0.68774007118557290171_pf, 0.13118241020105280626_pf, &
                                    0.92161977504885189358_pf, 0.26987577187133640373_pf, &
                                    0.92161977504885189358_pf, 0.13118241020105280626_pf, &
                                   -0.68774007118557290171_pf /)

  REAL(pf), DIMENSION(7) :: c = (/  0.94413392188212619249_pf, 0.34626230516255218639_pf, &
                                    0.93479137012319657440_pf, 0.50000000000000000000_pf, &
                                    0.06520862987680341172_pf, 0.65373769483744781361_pf, &
                                    0.05586607811787376587_pf /)

  REAL(pf), DIMENSION(7) :: b = (/ -0.03842134054164531021_pf, 0.08575888644805676475_pf, &
                                    0.06009756279830340969_pf, 0.13493788593566818923_pf, &
                                    0.86152221225054848031_pf, 0.04542352375299604783_pf, &
                                    -0.64931873064392764405_pf /)
  
  REAL(pf), DIMENSION(6,6) :: aij = reshape( &
  (/ 0.41118026824255338170_pf,  0.0_pf,                     0.0_pf,                     0.0_pf, 0.0_pf, 0.0_pf, &
     0.00642524721174115524_pf,  0.07720466121490930644_pf,  0.0_pf,                     0.0_pf, 0.0_pf, 0.0_pf, &
     0.30544869505114113917_pf,  0.02016768134753035846_pf, -0.40071232472612250408_pf,  0.0_pf, 0.0_pf, 0.0_pf, &
     0.60447214289054107539_pf, -0.03686929851984859646_pf, -0.80142464945224500816_pf, -0.11733965661499358435_pf, &
     0.0_pf,                     0.0_pf,                     0.19971712185972889664_pf,  0.04033536269506071692_pf, &
    -0.25902462499350481506_pf,  0.04149007905997619677_pf,  0.54240002445874024861_pf,  0.0_pf, &
     0.61089739010228227833_pf, -0.03809487697701307435_pf, -0.81003492990269199137_pf, -0.11986098498218263064_pf, &
     -0.00861028045044699535_pf, -0.07843023967207378433_pf /), &
     shape(aij), order=(/2,1/))

  
  TYPE, EXTENDS(integracao) :: integracao_rkn671

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
!   15 de marco de 2024
!
! Autoria:
!   oap
!
FUNCTION metodo (self, R, P, FSomas_ant)

  IMPLICIT NONE
  class(integracao_rkn671), INTENT(IN) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, FSomas_prox
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo
  
  INTEGER :: a, i, j

  REAL(pf), DIMENSION(7, self%N, self%dim) :: y, fi

  R1 = R + self % h * P * self % massasInvertidas
  P1 = P

  DO i = 1, 7

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

  ! Calcula as novas forcas
  FSomas_prox = self%forcas(R1)

  metodo(1,:,:) = R1
  metodo(2,:,:) = P1
  metodo(3,:,:) = FSomas_prox

END FUNCTION metodo

END MODULE rkn671