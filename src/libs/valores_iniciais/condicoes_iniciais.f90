! *****************************************************************
!! METODOS DE CONDICIONAMENTO DE VALORES INICIAIS
!
! Objetivos:
!   Funcoes para o condicionamento de valores iniciais.
! 
! Modificado:
!   17 de dezembro de 2025
! 
! Autoria:
!   oap
! 
MODULE condicoes_iniciais
  use iso_fortran_env, only: output_unit
  USE tipos
  USE utilidades
  USE condicionamento ! Metodos gerais de restricao
  USE sorteio_mod
  IMPLICIT NONE

  PUBLIC gerar_condicionar, gerar, condicionar
CONTAINS

! ************************************************************
!! Geracao de valores iniciais e condicionamento
!
! Objetivos:
!   Subrotina que sorteia valores iniciais e os condiciona.
!
! Modificado:
!   05 de janeiro de 2026
!
! Autoria:
!   oap
!
SUBROUTINE gerar_condicionar (sorteio_infos, massas, posicoes, momentos)
  TYPE(sorteio), POINTER, INTENT(INOUT) :: sorteio_infos
  REAL(pf), INTENT(INOUT), ALLOCATABLE  :: massas(:), posicoes(:,:), momentos(:,:)

  WRITE (output_unit,'(a)') "GERACAO DOS VALORES INICIAIS ("//sorteio_infos % modo//")"

  ! Gera
  CALL gerar(sorteio_infos, massas, posicoes, momentos)

  ! Condiciona
  CALL condicionar(sorteio_infos, massas, posicoes, momentos)

END SUBROUTINE gerar_condicionar

! ************************************************************
!! Geracao de valores iniciais
!
! Objetivos:
!   Subrotina que sorteia valores iniciais. Esta nao faz o
!   condicionamento.
!
! Modificado:
!   05 de janeiro de 2026
!
! Autoria:
!   oap
!
SUBROUTINE gerar (sorteio_infos, massas, posicoes, momentos)
  TYPE(sorteio), POINTER, INTENT(INOUT) :: sorteio_infos
  REAL(pf), INTENT(INOUT), ALLOCATABLE  :: massas(:), posicoes(:,:), momentos(:,:)
  INTEGER                               :: N
  REAL(pf)                              :: G, eps, ed

  N   = sorteio_infos % N
  G   = sorteio_infos % G
  eps = sorteio_infos % amortecedor

  WRITE (output_unit,'(a)') "  > gerando valores..."
  
  ALLOCATE(massas(N))
  ALLOCATE(posicoes(N,3))
  ALLOCATE(momentos(N,3))

  ! Sorteio
  CALL gerar_valores(N, sorteio_infos, massas, posicoes, momentos)

END SUBROUTINE

! ************************************************************
!! Condicionamento geral
!
! Objetivos:
!   Subrotina geral para condicionamento. Todas as outras 
!   rotinas de condicionamento sao chamadas a partir desta.
!
! Modificado:
!   05 de janeiro de 2026
!
! Autoria:
!   oap
!
SUBROUTINE condicionar (sorteio_infos, massas, posicoes, momentos)
  TYPE(sorteio), POINTER, INTENT(INOUT) :: sorteio_infos
  REAL(pf), INTENT(INOUT), ALLOCATABLE  :: massas(:), posicoes(:,:), momentos(:,:) ! out
  REAL(pf)                              :: G, eps, ed
  REAL(pf), DIMENSION(:), ALLOCATABLE   :: pd, jd

  G   = sorteio_infos % G
  eps = sorteio_infos % amortecedor

  ! Captura as integrais primeiras
  ed = sorteio_infos % ed
  pd = sorteio_infos % jd
  jd = sorteio_infos % pd

  ! Agora faz o condicionamento de acordo com o metodo desejado
  WRITE (output_unit,*) ""
  WRITE (output_unit,'(a)') "  > condicionando..."
  
  SELECT CASE (TRIM(sorteio_infos % modo))
    ! Sem nenhum condicionamento
    CASE("sorteio")
      ! Apenas nao faz nada

    ! Condicionamento por integrais primeiras iterativamente
    CASE("sorteio_ip_iterativo")
      CALL condicionar_ip_iterativo(G, massas, posicoes, momentos, eps, ed, jd, pd, 50)

    ! Condicionamento por integrais primeiras diretamente
    CASE("sorteio_ip_direto")
      CALL condicionar_ip_direto(G, massas, posicoes, momentos, eps, ed, jd, pd)

    ! Condicionamento de Aarseth
    CASE("sorteio_aarseth")
      CALL condicionar_aarseth(G, massas, posicoes, momentos, eps)

    ! Condicionamento de Aarseth Modificado
    CASE("sorteio_aarseth_modificado")
      CALL condicionar_aarseth_modificado(G, massas, posicoes, momentos, eps, jd)

    ! Outro caso: erro
    CASE DEFAULT
      ERROR STOP "!! Modo de condicionamento desconhecido. !!"
  END SELECT

  ! Outputs com informacoes sobre os dados iniciais sorteados
  CALL condicionamento_outputs(G, massas, posicoes, momentos, eps)
END SUBROUTINE

SUBROUTINE condicionamento_outputs (G, massas, posicoes, momentos, eps)
  REAL(pf), INTENT(IN) :: G, massas(:), posicoes(:,:), momentos(:,:), eps
  
  ! informacoes para exibir
  REAL(pf) :: potencial, cinetica
  REAL(pf) :: lintot(3), angtot(3)
  REAL(pf) :: inercia, dilatacao
  REAL(pf) :: anitenine ! anisotropia do tensor de inercia
  REAL(pf) :: f_prod_q, virial
  INTEGER  :: a
  REAL(pf) :: dist, mais_distante

  cinetica  = energia_cinetica(massas, momentos)
  
  virial = cinetica + cinetica
  IF (eps == 0) THEN
    potencial = energia_potencial(G, massas, posicoes, eps)
    virial = virial + potencial
  ELSE
    f_prod_q = virial_potencial_amortecido(G, massas, posicoes, eps, potencial)
    virial = virial + f_prod_q
  ENDIF

  lintot    = momento_linear_total(momentos)
  angtot    = momento_angular_total(posicoes, momentos)
  
  inercia   = momento_inercia(massas, posicoes)
  dilatacao = momento_dilatacao(posicoes, momentos)

  anitenine = anisotropia_tensor_inercia(massas, posicoes)

  mais_distante = 0
  DO a = 1, SIZE(massas)
    dist = NORM2(posicoes(a,:))
    IF (dist > mais_distante) mais_distante = dist
  END DO

  WRITE (output_unit,*) ""
  WRITE (output_unit,'(a)') "  > valores iniciais condicionados!"
  WRITE (output_unit,*) "     * V   = ", potencial
  WRITE (output_unit,*) "     * T   = ", cinetica
  WRITE (output_unit,*) "     * E   = ", potencial + cinetica
  WRITE (output_unit,*) "     * Vir = ", virial
  WRITE (output_unit,*) "     * J   = ", angtot
  WRITE (output_unit,*) "     * P   = ", lintot
  WRITE (output_unit,*) "     * I   = ", inercia
  WRITE (output_unit,*) "     * D   = ", dilatacao
  WRITE (output_unit,*) "     * A_I = ", anitenine
  WRITE (output_unit,*) "     * MD  = ", mais_distante

  CALL FLUSH(OUTPUT_UNIT)
END SUBROUTINE

END MODULE condicoes_iniciais