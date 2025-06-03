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
!   17 de setembro de 2024
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

  REAL(pf128), DIMENSION(7)   :: g_128, c_128, b_128
  REAL(pf128), DIMENSION(6,6) :: aij_128
  
  TYPE, EXTENDS(integracao) :: integracao_rkn671

    CONTAINS
      PROCEDURE :: metodo, metodo_mi, atualizar_constantes

  END TYPE

CONTAINS

SUBROUTINE atualizar_constantes (self)
  IMPLICIT NONE
  class(integracao_rkn671), INTENT(IN) :: self

  ! Atualiza a constante c    
  c_128(4) = 0.5_pf128
  c_128(5) = 0.06520862987680341024_pf128
  c_128(3) = 1.0_pf128 - 0.06520862987680341024_pf128
  c_128(6) = 0.65373769483744778901_pf128
  c_128(2) = 1.0_pf128 - 0.65373769483744778901_pf128
  c_128(7) = 0.05586607811787376572_pf128
  c_128(1) = 1.0_pf128 - 0.05586607811787376572_pf128

  ! Atualiza a constante aij
  aij_128(1,1) =  0.411180268242553396631920556447064423_pf128      
  aij_128(2,1) =  6.42524721174117304546999945731517070E-0003_pf128
  aij_128(2,2) =  7.72046612149093048359128245096028702E-0002_pf128
  aij_128(3,1) =  0.305448695051141170626198463759854879_pf128      
  aij_128(3,2) =  2.01676813475303541865157066166928855E-0002_pf128
  aij_128(3,3) = -0.400712324726122545078139295292374502_pf128      
  aij_128(4,1) =  0.604472142890541168206926928062394563_pf128      
  aij_128(4,2) = -3.68692985198485964628814112762170932E-0002_pf128
  aij_128(4,3) = -0.801424649452245090156278590584749003_pf128      
  aij_128(4,4) = -0.117339656614993593462036188533278538_pf128      
  aij_128(5,1) =  0.199717121859728944620476371072645335_pf128      
  aij_128(5,2) =  4.03353626950607083730314132333857710E-0002_pf128
  aij_128(5,3) = -0.259024624993504874689062107754851698_pf128      
  aij_128(5,4) =  4.14900790599761918687210201176132095E-0002_pf128
  aij_128(5,5) =  0.542400024458740215467216482829897305_pf128      
  aij_128(6,1) =  0.610897390102282341252396927519709759_pf128      
  aij_128(6,2) = -3.80948769770130725002204085208673273E-0002_pf128
  aij_128(6,3) = -0.810034929902692084723923846393547217_pf128      
  aij_128(6,4) = -0.119860984982182642862914255823422648_pf128      
  aij_128(6,5) = -8.61028045044699456764525580879826849E-0003_pf128
  aij_128(6,6) = -7.84302396720737808732518217542530983E-0002_pf128

  ! Atualiza a constante b
  b_128(1) = -3.84213405416452802288015362401451317E-0002_pf128
  b_128(2) =  8.57588864480567573165157066166928894E-0002_pf128
  b_128(3) =  6.00975627983034017118607047076255281E-0002_pf128
  b_128(4) =  0.134937885935668201865000000000000000_pf128
  b_128(5) =  0.861522212250548491868139295292374513_pf128
  b_128(6) =  4.54235237529960489434842933833071184E-0002_pf128
  b_128(7) = -0.649318730643927621481198463759854890_pf128
  b = b_128 * self % m_esc_128

  ! Atualiza a constante g
  g_128(1) = -0.68774007118557290171_pf128
  g_128(2) =  0.13118241020105280626_pf128
  g_128(3) =  0.92161977504885189358_pf128
  g_128(4) =  0.26987577187133640373_pf128
  g_128(5) =  0.92161977504885189358_pf128
  g_128(6) =  0.13118241020105280626_pf128
  g_128(7) = -0.68774007118557290171_pf128

  ! Se for massas iguais
  IF (self % mi) THEN
    c = c_128 * self % m_inv_128
    aij = aij_128 * self % m_esc_128
    g = g_128 * self % m2_128

  ! Se nao
  ELSE
    c = c_128
    aij = aij_128
    g = g_128
  ENDIF
END SUBROUTINE atualizar_constantes

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
  class(integracao_rkn671), INTENT(IN) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, FSomas_prox
  REAL(pf), DIMENSION(self%N, self%dim) :: P_prod
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo
  
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
  class(integracao_rkn671), INTENT(IN) :: self
  REAL(pf), DIMENSION(self%N, self%dim), INTENT(IN) :: R, P, FSomas_ant
  REAL(pf), DIMENSION(self%N, self%dim) :: R1, P1, FSomas_prox
  REAL(pf), DIMENSION(3, self%N, self%dim) :: metodo_mi
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
    FSomas_prox = self%forcas(yi)
    fi(i,:,:) = FSomas_prox

    ! Adiciona nas posicoes e momentos
    R1 = R1 + (self%h * self%h) * (b(i) * FSomas_prox)
    P1 = P1 + self % h * (g(i) * FSomas_prox)
  END DO

  metodo_mi(1,:,:) = R1
  metodo_mi(2,:,:) = P1
  metodo_mi(3,:,:) = 0.0_pf

END FUNCTION metodo_mi

END MODULE rkn671