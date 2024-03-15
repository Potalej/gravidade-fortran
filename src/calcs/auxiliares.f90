! *****************************************************************
!! FUNCOES AUXILIARES GERAIS
!
! Objetivos:
!   Este arquivo contem helpers de varias coisas, como vetores e
!   tal.
!
! Modificado:
!   15 de marco de 2024
! 
! Autoria:
!   oap
! 
MODULE auxiliares
  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  IMPLICIT NONE
  EXTERNAL dgesv
CONTAINS

! ************************************************************
!! Produto vetorial
!
! Objetivos:
!   Calcula o produto vetorial entre dois vetores do R3.
!
! Modificado:
!   03 de fevereiro de 2024
!
! Autoria:
!   oap
! 
FUNCTION produto_vetorial (u, v)

  IMPLICIT NONE
  REAL(pf), DIMENSION(:), INTENT(IN) :: u, v
  REAL(pf), DIMENSION(3)             :: produto_vetorial

  produto_vetorial(:) = 0.0_pf

  produto_vetorial(1) =  u(2)*v(3)-v(2)*u(3)
  produto_vetorial(2) = -u(1)*v(3)+v(1)*u(3)
  produto_vetorial(3) =  u(1)*v(2)-v(1)*u(2)

END FUNCTION produto_vetorial

! ************************************************************
!! Determinante de matriz 3x3
!
! Objetivos:
!   Calcula o determinante de uma matriz 3x3.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
REAL FUNCTION determinante (M) RESULT (det)

  REAL(pf), DIMENSION(3,3), INTENT(IN) :: M
  det = M(1,1)*(-M(3,2)*M(2,3)) - M(1,2)*(M(2,1)*M(3,3)-M(2,3)*M(3,1)) + M(1,3)*(M(2,1)*M(3,2)-M(2,2)*M(3,1))
  det = det + M(1,1)*M(2,2)*M(3,3)

END FUNCTION determinante

! ************************************************************
!! Tensor de inercia
!
! Objetivos:
!   Dada a massa e a posicao de um corpo, calcula o tensor de
!   inercia respectivo.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
FUNCTION tensor_inercia (m, R)

  IMPLICIT NONE
  REAL(pf), DIMENSION(3), INTENT(IN) :: R
  REAL(pf), INTENT(IN)               :: m
  REAL(pf), DIMENSION(3,3) :: tensor_inercia
  INTEGER              :: a, b

  DO a = 1, 3
    DO b = 1, 3
      IF (a /= b) THEN
        tensor_inercia(a,b) = m * R(a) * R(b)
      ENDIF
    END DO
  END DO

  tensor_inercia(1,1) = - m * (R(2)**2 + R(3)**2)
  tensor_inercia(2,2) = - m * (R(1)**2 + R(3)**2)
  tensor_inercia(3,3) = - m * (R(1)**2 + R(2)**2)

END FUNCTION tensor_inercia

! ************************************************************
!! Tensor de inercia geral
!
! Objetivos:
!   Dadas as massas e posicoes dos corpos, calcula o tensor de
!   inercia geral, que eh basicamente a soma dos tensores de
!   inercia individuais.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
FUNCTION tensor_inercia_geral (massas, posicoes)

  IMPLICIT NONE
  REAL(pf), DIMENSION(:), INTENT(IN) :: massas
  INTEGER :: a
  REAL(pf), DIMENSION(SIZE(massas),3) :: posicoes
  REAL(pf), DIMENSION(3,3) :: tensor_inercia_geral

  tensor_inercia_geral(:,:) = 0.0_pf
  
  DO a = 1, SIZE(massas)
    tensor_inercia_geral = tensor_inercia_geral + tensor_inercia(massas(a), posicoes(a,:))
  END DO   

END FUNCTION tensor_inercia_geral

! ************************************************************
!! Sistema linear de 3 equacoes (dgesv)
!
! Objetivos:
!   Resolve um sistema de 3 equacoes lineares usando o dgesv
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
FUNCTION sistema_linear3 (A, b)

  IMPLICIT NONE
  REAL(pf) :: A(3,3), b(3)
  REAL(pf) :: sistema_linear3(3)
  INTEGER  :: PIVOS(3), INFO

  CALL dgesv(3,1,A,3,PIVOS,b,3,INFO)

  IF (INFO == 0) THEN
    sistema_linear3 = b
  ELSE IF (INFO < 0) THEN
    WRITE (*, '(a)') 'O ', -INFO, '-ÉSIMO PARAMETRO TEM UM VALOR ILEGAL'
  ELSE
    WRITE (*, '(a)') 'MATRIZ SINGULAR! SEM SOLUCAO'
  ENDIF
END FUNCTION sistema_linear3

! ************************************************************
!! Centro de massas
!
! Objetivos:
!   Dadas as massas e posicoes dos corpos, calcula o centro de
!   massas, dado pela media dos corpos ponderada pelas massas.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
FUNCTION centro_massas (massas, posicoes)

  IMPLICIT NONE
  REAL(pf), INTENT(IN)   :: massas(:), posicoes(:,:)
  REAL(pf), DIMENSION(3) :: centro_massas
  INTEGER                :: a

  centro_massas(:) = 0.0_pf

  DO a = 1, SIZE(massas)
    centro_massas = centro_massas + massas(a) * posicoes(a,:)
  END DO
  centro_massas = centro_massas / SUM(massas)

END FUNCTION centro_massas

! ************************************************************
!! Momento linear total
!
! Objetivos:
!   Soma os momentos lineares e retorna tal vetor.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
FUNCTION momentoLinear_total (momentos)

  IMPLICIT NONE
  REAL(pf), INTENT(IN)   :: momentos(:,:)
  REAL(pf), DIMENSION(3) :: momentoLinear_total
  INTEGER                :: a
  momentoLinear_total = 0.0_pf
  DO a = 1, SIZE(momentos,1)
    momentoLinear_total = momentoLinear_total + momentos(a,:)
  END DO

END FUNCTION momentoLinear_total

END MODULE auxiliares