! ************************************************************
!! INTEGRADOR NUMERICO
!
! Objetivos:
!   Estrutura basica de um integrador numerico. As configuracoes
!   gerais de um metodo devem vir aqui, como tamanho de passo,
!   dimensao, massas, se corrige ou nao, se colide ou nao, etc.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
MODULE integrador

  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64
  IMPLICIT NONE
  PRIVATE
  PUBLIC integracao

  TYPE :: integracao
    ! m: Massas
    REAL(pf), ALLOCATABLE :: m(:)

    ! h: Passo de integracao
    ! G: Constante de gravitacao
    REAL(pf) :: h, G

    ! dim: Dimensao do problema
    ! N: Quantidade de partículas
    INTEGER :: dim = 3, N

    ! Se vai ou nao corrigir
    LOGICAL :: corrigir = .FALSE.

    ! Se vai ou nao colidir
    LOGICAL :: colidir = .FALSE.

    ! vetores para aplicar a correcao
    REAL(pf), ALLOCATABLE :: grads(:,:), gradsT(:,:), vetorCorrecao(:)
  END type integracao

END MODULE integrador