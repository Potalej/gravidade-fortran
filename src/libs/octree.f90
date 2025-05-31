! ************************************************************
!! ARVORE BINARIA
!
! Objetivos:
!   Geracao de uma arvore octonaria a partir de um vetor Nxk,
!   sendo N a quantidade de corpos, e de um indice 1<=i<=k.
!
! Modificado:
!   16 de marco de 2025
!
! Autoria:
!   oap
!  
MODULE octree
  USE tipos
  IMPLICIT NONE
  PRIVATE
  PUBLIC gerar_octree, arvore_octo, verificar_colisao_octree, localizar

  TYPE :: arvore_octo
    REAL(pf)              :: pos(3), raio
    INTEGER               :: N, corpos=0, max_corpos=2
    INTEGER               :: capacidade=4
    LOGICAL               :: folha = .TRUE.
    INTEGER, ALLOCATABLE  :: indices(:)
    TYPE(arvore_octo), ALLOCATABLE :: filhos(:)
  END TYPE arvore_octo

CONTAINS

! ************************************************************
!! Gera octree
!
! Objetivos:
!   A partir de um conjunto de massas e posicoes, gera um octree
!
! Modificado:
!   16 de marco de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE gerar_octree (arvore, N, massas, posicoes)

    TYPE(arvore_octo), INTENT(INOUT), ALLOCATABLE :: arvore
    INTEGER, INTENT(IN)  :: N
    REAL(pf), INTENT(IN) :: massas(N), posicoes(N,3)
    INTEGER  :: a
    REAL(pf) :: qc(3), max_dist

    IF (ALLOCATED(arvore)) THEN
        DEALLOCATE(arvore)
        ALLOCATE(arvore)
    ELSE
        ALLOCATE(arvore)
    ENDIF

    ! Media das posicoes para centralizar o octree
    qc = (/0.0_pf,0.0_pf,0.0_pf/)
    DO a = 1, N
        qc = qc + massas(a) * posicoes(a,:)
    END DO
    qc = qc / SUM(massas)

    ! Captura o corpo que esta mais distante
    max_dist = MAXVAL(ABS(posicoes)) + 1

    ! Cria o octree
    arvore % N = N
    arvore % raio = max_dist
    arvore % pos = qc
    
    ! Agora adiciona os corpos na arvore
    DO a = 1, N
        CALL adiciona_corpo(arvore, a, posicoes)
    END DO

END SUBROUTINE gerar_octree

! ************************************************************
!! Verifica colisao em um octree
!
! Objetivos:
!   Dadas uma arvore, um indice, as posicoes e os raios dos corpos,
!   verifica se o corpo `indice` colidiu com alguem.
!
! Modificado:
!   16 de marco de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE verificar_colisao_octree (arvore, indice, posicoes, raios, colisoes, colididos)

    TYPE(arvore_octo), INTENT(INOUT) :: arvore
    REAL(pf), INTENT(IN)   :: posicoes(arvore%N,3), raios(arvore%N)
    INTEGER, INTENT(IN)    :: indice
    INTEGER, INTENT(INOUT) :: colisoes(arvore%N-1), colididos
    INTEGER                :: i, ind, qntd_corpos
    LOGICAL                :: esta_no_cubo

    esta_no_cubo = ANY(arvore%indices==indice)

    ! Se for uma folha (corpo) e nao for o proprio indice, eh um corpo
    IF (arvore % folha) THEN
        DO i = 1, size(arvore%indices)
            IF (arvore%indices(i) == 0) THEN
                EXIT
            ELSE IF (arvore%indices(i) == indice) THEN
                CYCLE
            ENDIF

            ind = arvore % indices(i)
            IF (NORM2(posicoes(indice,:) - posicoes(ind,:)) <= raios(indice) + raios(ind)) THEN
                colididos = colididos + 1
                colisoes(colididos) = ind
            ENDIF
        END DO

    ! Se o indice fizer parte do cubo e tiver mais gente com ele, percorre os filhos
    ELSE IF (esta_no_cubo .AND. arvore % corpos > 1) THEN
        DO i=1,8
            IF (arvore%filhos(i)%corpos > 0) THEN
                CALL verificar_colisao_octree(arvore%filhos(i), indice, posicoes, raios, colisoes, colididos)
            ENDIF
        END DO

    ! Se o indice nao fizer parte do cubo mas houver interseccao, precisa ver o cubo
    ELSE IF (NORM2(posicoes(indice,:) - arvore%pos) <= raios(indice) + SQRT(3.0_pf) * arvore%raio) THEN

        ! Se o cubo for um corpo
        IF (arvore % folha) THEN
            DO i = 1, size(arvore%indices)
                IF (arvore%indices(i) == 0) THEN
                    EXIT
                ELSE IF (arvore%indices(i) == indice) THEN
                    CYCLE
                ENDIF

                ind = arvore % indices(i)
                IF (NORM2(posicoes(indice,:) - posicoes(ind,:)) <= raios(indice) + raios(ind)) THEN
                    colididos = colididos + 1
                    colisoes(colididos) = ind
                ENDIF
            END DO
        
        ! Se nao, percorre os filhos
        ELSE
            DO i=1,8
                IF (arvore%filhos(i)%corpos > 0) THEN
                    CALL verificar_colisao_octree(arvore%filhos(i), indice, posicoes, raios, colisoes, colididos)
                ENDIF
            END DO
        ENDIF
    ENDIF

END SUBROUTINE verificar_colisao_octree

! ************************************************************
!! Identifica em qual folha de um galho um corpo esta
!
! Objetivos:
!   Identificao de em qual folha de um galho um dado corpo esta.
!
! Modificado:
!   16 de marco de 2025
!
! Autoria:
!   oap
! 
FUNCTION qual_octo (arvore, indice, posicoes)
    
    TYPE(arvore_octo), INTENT(INOUT) :: arvore
    INTEGER  :: indice
    REAL(pf) :: posicoes(arvore%N,3)
    INTEGER  :: posx, posy, posz, qual_octo
    
    posx=0
    posy=1
    posz=0

    ! Verifica se esta no polo sul
    IF (posicoes(indice,3) < arvore % pos(3)) THEN
        posz = 4
    ENDIF
    ! Verifica se esta no sul
    IF (posicoes(indice,2) < arvore % pos(2)) THEN
        posy = 3
    ENDIF
    ! Verifica se esta no leste
    IF (posicoes(indice,1) >= arvore % pos(1) ) THEN
        posx = 1
    ENDIF

    qual_octo = posx + posy + posz

END FUNCTION qual_octo

! ************************************************************
!! Adiciona indices
!
! Objetivos:
!   Adiciona um indice em uma folha ou galho.
!
! Modificado:
!   16 de marco de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE adiciona_indice (arvore, indice)
    TYPE(arvore_octo), INTENT(INOUT) :: arvore
    INTEGER  :: indice
    INTEGER  :: capacidade, capacidade_antiga
    INTEGER, ALLOCATABLE :: temp(:)

    ! Se ainda nao tiver corpos
    IF (.NOT. ALLOCATED(arvore%indices)) THEN
        ALLOCATE(arvore%indices(2))
        arvore%capacidade=2
        arvore%indices = 0
        arvore%indices(1) = indice
        arvore%corpos = 1
    
    ! Verifica se ja contem o indice
    ELSE IF (ANY(arvore % indices == indice)) THEN
        RETURN

    ! Se ja tiver mas tiver atingido a capacidade maxima, aumenta
    ELSE IF(arvore%corpos == arvore%capacidade) THEN
        ! Aloca a variavel temporaria e armazena os dados
        ALLOCATE(temp(arvore%capacidade))
        temp = arvore%indices
        capacidade_antiga = arvore%capacidade
        
        ! Agora apaga os dados antigos
        DEALLOCATE(arvore%indices)
        
        ! E aloca novamente com o dobro de espacos
        arvore%capacidade = 2 * capacidade_antiga
        ALLOCATE(arvore%indices(arvore%capacidade))
        arvore%indices(1:capacidade_antiga) = temp
        
        DEALLOCATE(temp)

        ! Adiciona o novo elemento
        arvore%indices(capacidade_antiga+1) = indice
        arvore%corpos = arvore%corpos + 1

    ! Se nao, so adiciona
    ELSE
        arvore%corpos = arvore%corpos + 1
        arvore%indices(arvore%corpos) = indice
    ENDIF
    
END SUBROUTINE adiciona_indice

! ************************************************************
!! Adiciona um corpo em uma arvore
!
! Objetivos:
!   Dada uma arvore, adiciona um corpo nela.
!
! Modificado:
!   16 de marco de 2025
!
! Autoria:
!   oap
! 
RECURSIVE SUBROUTINE adiciona_corpo (arvore, indice, posicoes)

    TYPE(arvore_octo), INTENT(INOUT) :: arvore
    INTEGER  :: indice
    INTEGER  :: filho
    REAL(pf) :: posicoes(arvore%N,3)

    ! Salva o indice
    CALL adiciona_indice(arvore, indice)

    ! Se for uma folha e precisar subdividir
    IF (arvore % folha .AND. arvore%corpos == arvore%max_corpos) THEN
        CALL subdividir(arvore, posicoes)

    ! Se nao for uma folha, precisa adicionar na folha certa
    ELSE IF (.NOT. arvore % folha) THEN
        filho = qual_octo(arvore, indice, posicoes)
        CALL adiciona_corpo (arvore % filhos(filho), indice, posicoes)
    ENDIF

END SUBROUTINE adiciona_corpo

! ************************************************************
!! Subdivisao de um galho
!
! Objetivos:
!   Dada uma folha, a transforma em galho criando as subdivisoes
!   necessarias e espalhando os corpos presentes nos devidos
!   locais.
!
! Modificado:
!   16 de marco de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE subdividir (arvore, posicoes)

    TYPE(arvore_octo), INTENT(INOUT) :: arvore
    REAL(pf), INTENT(IN) :: posicoes(arvore%N,3)
    REAL(pf) :: subraio, xmin, xmax, ymin, ymax, zmin, zmax
    INTEGER :: i

    ! Metade do raio
    subraio = arvore%raio/2

    ! Cria 8 filhos
    ALLOCATE(arvore%filhos(8))
    arvore%filhos(:)%N = arvore%N
    arvore%filhos(:)%raio = subraio
    arvore%filhos(:)%folha = .TRUE.

    ! novas posicoes
    xmin = arvore%pos(1) - subraio
    xmax = arvore%pos(1) + subraio
    ymin = arvore%pos(2) - subraio
    ymax = arvore%pos(2) + subraio
    zmin = arvore%pos(3) - subraio
    zmax = arvore%pos(3) + subraio

    ! Define as posicoes de cada filho
    arvore%filhos(1)%pos = (/xmin, ymax, zmax/) ! NON
    arvore%filhos(2)%pos = (/xmax, ymax, zmax/) ! NEN
    arvore%filhos(3)%pos = (/xmin, ymin, zmax/) ! SON
    arvore%filhos(4)%pos = (/xmax, ymin, zmax/) ! SEN
    arvore%filhos(5)%pos = (/xmin, ymax, zmin/) ! NOS
    arvore%filhos(6)%pos = (/xmax, ymax, zmin/) ! NES
    arvore%filhos(7)%pos = (/xmin, ymin, zmin/) ! SOS
    arvore%filhos(8)%pos = (/xmax, ymin, zmin/) ! SES

    ! Agora eh um galho
    arvore%folha = .FALSE.

    ! Insere os corpos nos filhos
    DO i=1, arvore%corpos
        CALL adiciona_corpo(arvore, arvore%indices(i), posicoes)
    ENDDO
END SUBROUTINE subdividir

! ************************************************************
!! Localiza um indice em uma arvore
!
! Objetivos:
!   Localiza recursivamente um indice em uma arvore. Serve mais
!   para debug.
!
! Modificado:
!   16 de marco de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE localizar (arvore, indice, n, posicoes)
    TYPE(arvore_octo), INTENT(INOUT) :: arvore
    INTEGER, INTENT(IN)  :: indice, N
    REAL(pf), INTENT(IN) :: posicoes(N,3)
    INTEGER :: qual, filho

    IF (.NOT. arvore % folha) THEN
        DO filho = 1, 8
            IF (arvore % filhos(filho) % corpos > 0) THEN
                IF (ANY(arvore % filhos(filho) % indices == indice)) THEN
                    WRITE(*,*) filho
                    CALL localizar(arvore % filhos(filho), indice, n, posicoes)
                ENDIF
            ENDIF
        END DO
    ENDIF

END SUBROUTINE

END MODULE octree