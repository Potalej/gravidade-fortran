! ************************************************************
!! METODO NUMERICO: RKN55 (1)
!
! Objetivos:
!   Aplicacao do metodo simpletico RKN55 (1).
!   Runge-Kutta-Nystrom de Ordem 5 e 5 Estagios. Eh o primeiro
!   metodo na tabela 1.
!   Referencia: (Okunbor & Skeel, 1992, p.378)
!
! Modificado:
!   17 de setembro de 2024
!
! Autoria:
!   oap
!
MODULE rkn551
  USE tipos
  USE OMP_LIB
  USE integrador

  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao_rkn551

  REAL(pf), DIMENSION(5) :: g = (/ -1.67080892327314312060_pf, 1.22143909230997538270_pf, &
                                    0.08849515813253908125_pf, 0.95997088013770159876_pf, &
                                    0.40090379269297793385_pf /)

  REAL(pf), DIMENSION(5) :: c = (/  0.69491389107017931259_pf, 0.63707199676998338411_pf, &
                                   -0.02055756998211598005_pf, 0.79586189634575355001_pf, &
                                    0.30116624272377778837_pf /)

  REAL(pf), DIMENSION(5) :: b = (/ -0.5097405931666266_pf, 0.4432944508391433_pf, &
                                    0.09031440353892717_pf, 0.19596663503460832_pf, &
                                    0.28016510375392145_pf /)
  
  REAL(pf), DIMENSION(4,4) :: aij = reshape( &
    (/  0.09664275313578931_pf,  0.0_pf,                 0.0_pf,                   0.0_pf,                  &
        1.1954161014734481_pf,  -0.8032544610898866_pf,  0.0_pf,                   0.0_pf,                  &
       -0.16866482800105376_pf,  0.19395219080582085_pf, 0.07224916977516799_pf,   0.0_pf,                  &
        0.6578770843749833_pf,  -0.4102884193238952_pf,  0.028470999680411457_pf, -0.4748934220077836_pf/), &
       shape(aij), order=(/2,1/) )

  TYPE, EXTENDS(integracao) :: integracao_rkn551

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
  class(integracao_rkn551), INTENT(IN) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, FSomas_prox
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo
  
  INTEGER :: i, j

  REAL(pf), DIMENSION(SIZE(c), self%N, self%dim) :: y, fi

  R1 = R + self % h * P * self % massasInvertidas
  P1 = P

  DO i = 1, SIZE(c)

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
  class(integracao_rkn551), INTENT(IN) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, FSomas_prox
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo_mi
  
  INTEGER :: i, j

  REAL(pf), DIMENSION(SIZE(c), self%N, self%dim) :: y, fi

  R1 = R + self % h * P * self % m_inv
  P1 = P

  DO i = 1, SIZE(c)

    ! Calcula a base do yi
    y(i,:,:) = R + c(i) * self % h * P * self % m_inv

    ! Se i > 1, calcula os termos a mais
    IF (i .GT. 1) THEN
      ! Somatorio
      DO j = 1, i-1
        y(i,:,:) = y(i,:,:) + self % h * self % h * aij(i-1,j) * fi(j,:,:)
      END DO
    ENDIF

    ! Calcula fi
    FSomas_prox = self%forcas(y(i,:,:))
    fi(i,:,:) = FSomas_prox * self % m_esc

    ! Adiciona nas posicoes e momentos
    R1 = R1 + self % h * self % h * b(i) * fi(i,:,:)
    P1 = P1 + self % h * g(i) * (self % m2 * FSomas_prox)
  END DO

  metodo_mi(1,:,:) = R1
  metodo_mi(2,:,:) = P1
  metodo_mi(3,:,:) = 0.0_pf
END FUNCTION metodo_mi

END MODULE rkn551