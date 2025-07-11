! *****************************************************************
!! METODOS DE CONDICIONAMENTO DE VALORES INICIAIS
!
! Objetivos:
!   Funcoes para o condicionamento de valores iniciais.
! 
! Modificado:
!   02 de maio de 2025
! 
! Autoria:
!   oap
! 
MODULE condicoesIniciais
  USE tipos
  USE mecanica
  USE auxiliares
  USE condicionamento ! Metodos gerais de restricao
  USE json_utils_mod
  IMPLICIT NONE
CONTAINS

! ************************************************************
!! Condicionamento geral
!
! Objetivos:
!   Subrotina geral para condicionamento. Todas as outras sao
!   chamadas a partir desta.
!
! Modificado:
!   02 de maio de 2025
!
! Autoria:
!   oap
!
SUBROUTINE condicionar (dados, massas, posicoes, momentos, metodo)
  TYPE(json_value), POINTER, INTENT(IN) :: dados ! in
  TYPE(json_value), POINTER             :: sorteio
  CHARACTER(LEN=*), INTENT(IN)          :: metodo ! in
  REAL(pf), INTENT(INOUT), ALLOCATABLE  :: massas(:), posicoes(:,:), momentos(:,:) ! out
  REAL(pf) :: ed
  REAL(pf), DIMENSION(:), ALLOCATABLE :: pd, jd
  INTEGER  :: N
  REAL(pf) :: G
  LOGICAL  :: mi, encontrado ! massas iguais

  CALL json % get(dados, "sorteio", sorteio)
  CALL json % get(dados, "N", N)
  G = json_get_float(dados, "G")
  CALL json % get(dados, "massas_iguais", mi, encontrado)
  IF (.NOT. encontrado) mi = .FALSE.

  WRITE(*,'(a)') "GERACAO DAS VALORES INICIAIS ("//metodo//")"

  ! Sorteio dos valores na regiao e com distribuicao desejadas
  WRITE(*,'(a)') "  > gerando valores..."
  ALLOCATE(massas(N))
  ALLOCATE(posicoes(N,3))
  ALLOCATE(momentos(N,3))
  CALL gerar_valores(N, sorteio, massas, posicoes, momentos, mi)

  ! Captura as integrais primeiras
  ed = json_get_float(sorteio, "integrais.energia_total")
  pd = json_get_float_vec(sorteio, "integrais.linear_total")
  jd = json_get_float_vec(sorteio, "integrais.angular_total")

  ! Agora faz o condicionamento de acordo com o metodo desejado
  WRITE(*,*) ""
  WRITE(*,'(a)') "  > condicionando..."
  
  SELECT CASE (TRIM(metodo))
    ! Condicionamento por integrais primeiras iterativamente
    CASE("sorteio_ip_iterativo")
      CALL condicionar_ip_iterativo(G, massas, posicoes, momentos, ed, pd, jd)

    ! Condicionamento por integrais primeiras diretamente
    CASE("sorteio_ip_direto")
      CALL condicionar_ip_direto(G, massas, posicoes, momentos, ed, pd, jd)

    ! Condicionamento de Aarseth
    CASE("sorteio_aarseth")
      CALL condicionar_aarseth(G, massas, posicoes, momentos)

    ! Outro caso: erro
    CASE DEFAULT
      ERROR STOP "!! Modo de condicionamento desconhecido. !!"
  END SELECT

  ! Outputs com informacoes sobre os dados iniciais sorteados
  CALL condicionamento_outputs(G, massas, posicoes, momentos)
END SUBROUTINE

SUBROUTINE condicionamento_outputs (G, massas, posicoes, momentos)
  REAL(pf), INTENT(IN) :: G, massas(:), posicoes(:,:), momentos(:,:)
  
  ! informacoes para exibir
  REAL(pf) :: potencial, cinetica
  REAL(pf) :: lintot(3), angtot(3)
  REAL(pf) :: inercia, dilatacao
  REAL(pf) :: anitenine ! anisotropia do tensor de inercia

  potencial = energia_potencial(G, massas, posicoes)
  cinetica  = energia_cinetica(massas, momentos)
  lintot    = momentoLinear_total(momentos)
  angtot    = momento_angular_total(posicoes, momentos)
  
  inercia   = momento_inercia(massas, posicoes)
  dilatacao = momento_dilatacao(posicoes, momentos)

  anitenine = anisotropia_tensor_inercia(massas, posicoes)

  WRITE(*,*) ""
  WRITE(*,'(a)') "  > valores iniciais condicionados!"
  WRITE(*,*) "     * V   = ", potencial
  WRITE(*,*) "     * T   = ", cinetica
  WRITE(*,*) "     * E   = ", potencial + cinetica
  WRITE(*,*) "     * J   = ", angtot
  WRITE(*,*) "     * P   = ", lintot
  WRITE(*,*) "     * I   = ", inercia
  WRITE(*,*) "     * D   = ", dilatacao
  WRITE(*,*) "     * A_I = ", anitenine

END SUBROUTINE

END MODULE condicoesIniciais