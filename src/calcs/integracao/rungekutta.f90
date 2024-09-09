! ************************************************************
!! METODO NUMERICO: RUNGE-KUTTA (BASE)
!
! Objetivos:
!   Base para os metodos de Runge-Kutta. Os metodos basicos
!   aqui servem apenas para contruir outras classes que usem
!   de fato o metodo, como RK4 e RKF45.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
MODULE rungekutta
  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64

  IMPLICIT NONE
  PRIVATE
  PUBLIC RK

  TYPE :: RK

    ! m: Massas
    ! massasInvertidas : Matriz como inverso das massas para facilitar a integracao dos momentos
    REAL(pf), ALLOCATABLE :: m(:), massasInvertidas(:,:)

    ! h: Passo de integracao
    ! G: Constante de gravitacao
    ! potsoft: Softening do potencial
    REAL(pf) :: h, G, potsoft

    ! dim: dimensao do problema
    ! N: quantidade de particulas
    INTEGER :: dim = 3, N

    CONTAINS
      PROCEDURE :: Iniciar, forcas

  END TYPE

CONTAINS

! ************************************************************
!! Construtor da classe
!
! Objetivos:
!   Define o principal, salvando os valores.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE Iniciar (self, N, massas, G, h, potsoft)
  IMPLICIT NONE
  CLASS(RK), INTENT(INOUT) :: self
  REAL(pf), ALLOCATABLE :: massas(:)
  INTEGER, INTENT(INOUT) :: N
  REAL(pf)               :: G, h, potsoft
  INTEGER :: a, i

  ! quantidade de particulas
  self % N = N
  ! Softening do potencial
  self % potsoft = potsoft
  ! massas
  ALLOCATE(self % m (self % N))
  self % m = massas
  ! vetor de massas invertidas
  ALLOCATE(self % massasInvertidas (self % N, self % dim))
  DO a = 1, self % N
    DO i = 1, self % dim
      self % massasInvertidas(a,i) = 1/(massas(a))
    END DO
  END DO

  ! gravidade
  self % G = G
  ! passo
  self % h = h
END SUBROUTINE Iniciar  

! ************************************************************
!! Matriz de forcas
!
! Objetivos:
!   Calcula a matriz de forcas a partir das posicoes. Eh 
!   igual em todo metodo de Runge-Kutta.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
FUNCTION forcas (self, R)
  IMPLICIT NONE
  CLASS(RK), INTENT(IN) :: self
  REAL(pf), DIMENSION(self % N, self % dim), INTENT(IN) :: R
  REAL(pf), DIMENSION(self % dim) :: Fab
  INTEGER :: a, b
  REAL(pf) :: distancia
  REAL(pf), DIMENSION(self % N, self % dim) :: forcas
  
  forcas(:,:) = 0

  DO a = 2, self % N
    DO b = 1, a - 1
      ! distancia entre os corpos
      distancia = (norm2(R(b,:) - R(a,:))**2 + self % potsoft**2)**(3/2)
      ! forca entre os corpos a e b
      Fab = - self % G * self % m(a) * self % m(b) * (R(b,:) - R(a,:))/distancia
      ! Adiciona na matriz
      forcas(a,:) = forcas(a,:) - Fab
      forcas(b,:) = forcas(b,:) + Fab
    END DO
  END DO

END FUNCTION forcas

END MODULE rungekutta