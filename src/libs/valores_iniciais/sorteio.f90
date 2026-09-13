! *****************************************************************
!! SORTEIO DE VALORES INICIAIS
!
! Objetivos:
!   Funcoes para a geracao de valores iniciais
! 
! Modificado:
!   17 de dezembro de 2025
! 
! Autoria:
!   oap
! 
MODULE sorteio_mod
  USE tipos
  USE aleatorio
  IMPLICIT NONE

  TYPE :: sorteio_vetor
    CHARACTER(LEN=:), ALLOCATABLE :: distribuicao
    CHARACTER(LEN=:), ALLOCATABLE :: regiao
    REAL(pf) :: intervalo(2)
    LOGICAL  :: normalizado ! para uso nas massas
  END TYPE 

  TYPE :: sorteio

    ! Quantidade de corpos
    INTEGER :: N

    ! Constante de gravitacao
    REAL(pf) :: G

    ! Amortecedor
    REAL(pf) :: amortecedor

    ! Modo de sorteio
    CHARACTER(LEN=:), ALLOCATABLE :: modo

    ! Integrais primeiras desejadas
    REAL(pf) :: ed ! Energia desejada
    REAL(pf) :: jd(3) ! Angular total desejado
    REAL(pf) :: pd(3) ! Linear total desejado

    TYPE(sorteio_vetor) :: massas
    TYPE(sorteio_vetor) :: posicoes
    TYPE(sorteio_vetor) :: momentos

  END TYPE sorteio

CONTAINS

! ************************************************************
!! Gera valores
!
! Objetivos:
!   Gera valores iniciais aleatorios.
!
! Modificado:
!   16 de dezembro de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE gerar_valores (N, sorteio_infos, massas, posicoes, momentos)
  TYPE(sorteio), POINTER, INTENT(IN) :: sorteio_infos
  INTEGER, INTENT(IN)         :: N
  REAL(pf), INTENT(INOUT)     :: posicoes(:,:), momentos(:,:), massas(:)
  REAL(pf), ALLOCATABLE       :: intervalo(:)

  ! Gera massas
  WRITE (*,'(A)', ADVANCE='no') '    * gerando massas'
  massas = gerar_massas(N, sorteio_infos % massas)

  ! Gera as posições
  WRITE (*,'(A)', ADVANCE='no') '    * gerando posicoes'  
  posicoes = gerar_vetores3d(N, sorteio_infos % posicoes)

  ! Gera os momentos
  WRITE (*,'(A)', ADVANCE='no') '    * gerando momentos'
  momentos = gerar_vetores3d(N, sorteio_infos % momentos)

END SUBROUTINE gerar_valores

! ************************************************************
!! Gera vetores 3d
!
! Objetivos:
!   Gera vetores 3d utilizados para as posicoes e velocidades.
!
! Modificado:
!   17 de dezembro de 2025
!
! Autoria:
!   oap
! 
FUNCTION gerar_vetores3d (N, sorteio_infos) RESULT(vetores)

  INTEGER, INTENT(IN)         :: N
  TYPE(sorteio_vetor), INTENT(INOUT) :: sorteio_infos
  REAL(pf), DIMENSION(N,3)    :: vetores
  CHARACTER(:), ALLOCATABLE   :: distribuicao, regiao
  REAL(pf), ALLOCATABLE     :: intervalo(:)
  REAL(pf)                  :: raio, vmin, vmax, distmin, centro(3)

  distribuicao = sorteio_infos % distribuicao
  regiao = sorteio_infos % regiao

  WRITE (*,'(A)') ' ('//distribuicao//')'
  
  ! Intervalo de sorteio
  intervalo = sorteio_infos % intervalo
  vmin = intervalo(1)
  vmax = intervalo(2)
  raio = 0.5_pf*(vmax - vmin)
  centro(:) = 0.5_pf*(vmin + vmax)
  distmin = 0.0_pf
  ! Distancia minima
  IF (SIZE(intervalo) == 3) THEN
    distmin = intervalo(3)
  ENDIF

  SELECT CASE (TRIM(distribuicao))
    ! Uniforme (0,1)
    CASE ("uniforme"); CALL uniforme(vetores, N, distmin, vmin, vmax, regiao)
    ! Normal (0,1)
    CASE ("normal"); CALL normal(vetores, N, distmin, regiao, raio, centro)
    ! Cauchy
    CASE ("cauchy"); CALL cauchy(vetores, N, distmin, regiao, raio, centro)
  END SELECT

END FUNCTION gerar_vetores3d

! ************************************************************
!! Gera vetor de massas
!
! Objetivos:
!   Gera vetor de massas conforme um intervalo informado.
!
! Modificado:
!   17 de dezembro de 2025
!
! Autoria:
!   oap
! 
FUNCTION gerar_massas (N, sorteio_infos) RESULT(massas)

  INTEGER, INTENT(IN)           :: N
  TYPE(sorteio_vetor)           :: sorteio_infos
  REAL(pf), DIMENSION(N)        :: massas
  REAL(pf), ALLOCATABLE         :: intervalo(:)
  CHARACTER(LEN=:), ALLOCATABLE :: distribuicao
  REAL(pf)               :: vet_min(N)

  intervalo = sorteio_infos % intervalo
  vet_min = intervalo(1)

  ! 1 - Se a opcao "normalizadas" for TRUE, faz m = 1/N
  IF (sorteio_infos % normalizado) THEN
    massas = 1.0_pf / N
    WRITE(*,*)
  
  ! 2 - Se nao for normalizada mas forem iguais, determina
  ELSE IF (intervalo(1) == intervalo(2)) THEN
    massas = intervalo(1)
    WRITE(*,*)
    
  ELSE
    ! AQUI PRECISA APLICAR A DISTRIBUICAO DESEJADA. POR ENQUANTO SO TEM A UNIFORME
    distribuicao = sorteio_infos % distribuicao
    WRITE (*,'(A)') ' ('//TRIM(distribuicao)//')'

    IF (distribuicao == "uniforme") THEN
      CALL RANDOM_SEED()
      CALL RANDOM_NUMBER(massas)
      
      ! Agora condiciona no intervalo
      massas = massas * (intervalo(2) - intervalo(1)) + vet_min
    ENDIF
  ENDIF

END FUNCTION gerar_massas

END MODULE sorteio_mod