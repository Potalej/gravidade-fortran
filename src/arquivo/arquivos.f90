! *****************************************************************
!! ARQUIVO
!
! Objetivos:
!   Este arquivo contem helpers para arquivos gerais, como para a 
!   abertura, leitura e escrita.
!   
!   O objeto `arquivo` tem as propriedades de criacao, escrita, 
!   fechamento, criacao de nome e criacao de formato.
!   
! Modificado:
!   26 de maio de 2024
! 
! Autoria:
!   oap
! 
MODULE arquivos
  USE, INTRINSIC :: iso_fortran_env, only: pf=>real64

  ! OPENMP
  USE OMP_LIB

  IMPLICIT NONE
  PRIVATE
  PUBLIC arquivo, ler_csv, criar_dir, salvar_sorteio, espacosVazios, capturar_unidade, diretorio_out

  ! classe de arquivo
  TYPE :: arquivo

  ! ids dos arquivos
  INTEGER :: idarqdata, idarqinfo, idarqbkp, qntdCorpos_int, dimensao_int
  ! nome dos arquivos, qntd de corpos, formato e dimensao
  CHARACTER(:), allocatable :: dirarq, nomearqdata, nomearqinfo, nomearqbkp, formato, formatoMassas, qntdCorpos, dimensao
  ! extensao
  CHARACTER(4) :: extensao = '.csv'
  ! diretorio padrao (fora da pasta build)
  CHARACTER(11) :: dir = "./out/data/"
  CHARACTER(5) :: dir_out = "./out"
  CHARACTER(4) :: dir_data = "data"
  
  CONTAINS
    PROCEDURE :: criar, escrever, fechar, nomeArquivo, criarFormato, escrever_massas, escrever_cabecalho, &
                 diretorio_data, inicializar_arquivo_info, atualizar_arquivo_info, arquivo_bkp, excluir_bkp
  END TYPE

CONTAINS

! ************************************************************
!! Diretorio "out"
!
! Objetivos:
!   Verifica se existe o diretorio "out". Se nao existir, cria.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE diretorio_out ()
  IMPLICIT NONE
  LOGICAL :: existe = .TRUE.
  ! verifica se existe o diretorio out
  INQUIRE(file="./out", exist=existe)
  IF (.NOT. existe) THEN
    CALL criar_dir("./out")
  ENDIF
END SUBROUTINE diretorio_out

! ************************************************************
!! Diretorio "data"
!
! Objetivos:
!   Verifica se existe o diretorio "data". Se nao existir, cria.
!
! Modificado:
!   26 de maio de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE diretorio_data (self)
  IMPLICIT NONE
  CLASS(arquivo), INTENT(INOUT) :: self
  LOGICAL :: existe = .TRUE.
  CALL diretorio_out()
  ! verifica se existe o diretorio padrao
  INQUIRE(file=TRIM(self % dir), exist=existe)
  IF (.NOT. existe) THEN
    CALL criar_dir(TRIM(self % dir_data), TRIM(self % dir_out))
  ENDIF
END SUBROUTINE diretorio_data

! ************************************************************
!! Nome do arquivo
!
! Objetivos:
!   Cria o nome do diretorio baseado no dia corrente e contando
!   a partir de 1, ou seja, nao havendo pastas do mesmo dia
!   eh criado um AAAAMMDD_01, e caso haja eh criado AAAAMMDD_02, 
!   AAAAMMDD_03, etc.
!
! Modificado:
!   26 de maio de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE nomeArquivo (self)

  IMPLICIT NONE
  CLASS(arquivo), INTENT(INOUT) :: self

  ! para iterar e nao repetir arquivo
  INTEGER :: i = 1
  CHARACTER(3) :: numero
  LOGICAL :: existe

  ! para capturar a data
  CHARACTER(8) :: datahoje

  ! verifica se existe o diretorio padrao
  CALL self % diretorio_data()

  ! Por padrao, existe
  existe = .TRUE.

  ! em string
  CALL DATE_AND_TIME(datahoje)

  DO WHILE (existe)
    WRITE(numero, '(I3.3)') i
    i = i + 1

    ! cria nomes
    self % dirarq      = TRIM(datahoje)//"_"//TRIM(numero)
    self % nomearqdata = self % dir//self % dirarq//"/data"//self % extensao
    self % nomearqinfo = self % dir//self % dirarq//"/info.txt"
    self % nomearqbkp  = self % dir//self % dirarq//"/bkp.txt"

    ! verifica se existe
    INQUIRE(file=TRIM(self % dir//self % dirarq), exist=existe)
  END DO

END SUBROUTINE nomeArquivo

! ************************************************************
!! Formato do arquivo
!
! Objetivos:
!   Para uma determinada quantidade de particulas e dimensoes,
!   eh criada a formatacao para transformar o array em uma 
!   string corretamente. SENDo N:= qntdCorpos e D := dimensao,
!   por exemplo, deve ser gerada a seguinte formatacao:
!   '(2(N(D(F25.7,:,","))))'
!   enquanto para as massas:
!   '(N(F25.7,:,","))'
! 
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE criarFormato (self, qntdCorpos, dimensao)

  IMPLICIT NONE
  CLASS(arquivo), INTENT(INOUT) :: self
  INTEGER, INTENT(IN)           :: qntdCorpos, dimensao
  CHARACTER                     :: formato
  INTEGER                       :: i = 1

  ! salva a quantidade de corpos e dimensao
  self % qntdCorpos = espacosVazios(qntdCorpos)
  self % dimensao = espacosVazios(dimensao)

  self % qntdCorpos_int = qntdCorpos
  self % dimensao_int = dimensao

  self % formato = '(2(' // self % dimensao // '(' // self % qntdCorpos // '(F25.13, :, ","))))'
  self % formatoMassas = '(' // self % qntdCorpos // '(F25.7, :, ","))'

END SUBROUTINE criarFormato

! ************************************************************
!! Criacao do arquivo
!
! Objetivos:
!   Cria um arquivo .csv que se configura para fazer a formatacao
!   de uma determinada quantidade de particulasd e uma determinada
!   quantidade de dimensoes. O `idarq` eh utilizado como
!   identificador unico do arquivo.
!
! Modificado:
!   26 de maio de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE criar (self, qntdCorpos, dimensao)

  IMPLICIT NONE
  CLASS(arquivo), INTENT(INOUT) :: self
  INTEGER, INTENT(IN)           :: qntdCorpos, dimensao
  INTEGER(kind=4)               :: idarqdata, idarqinfo, idarqbkp
  LOGICAL                       :: existe ! para verificar se ja existe ou nao

  WRITE (*, '(a)') 'CRIAR ARQUIVO PARA SALVAR PLOT:'

  ! cria formatacao
  CALL self % criarFormato(qntdCorpos, dimensao)
  WRITE (*, '(a)') '  > formato : ' // self % formato

  ! criacao do nome do arquivo
  CALL self % nomeArquivo()
  WRITE (*,'(a)') '  > diretorio de saida: ' // self % dirarq

  ! cria o diretorio
  CALL criar_dir(self % dirarq, self % dir)

  ! cria o arquivo csv
  CALL capturar_unidade(idarqdata)
  self % idarqdata = idarqdata
  OPEN(idarqdata, file = self % nomearqdata, status='new')

  ! cria o arquivo info
  CALL capturar_unidade(idarqinfo)
  self % idarqinfo = idarqinfo
  OPEN(idarqinfo, file = self % nomearqinfo, status='new')

  ! cria o arquivo bkp
  CALL capturar_unidade(idarqbkp)
  self % idarqbkp = idarqbkp
  OPEN(idarqbkp, file = self % nomearqbkp, status='new')

  WRITE (*,'(a)') '  > arquivos criados!'
  WRITE (*,*)

END SUBROUTINE criar

! ************************************************************
!! Escrita das massas
!
! Objetivos:
!   Salva as massas no arquivo de acordo com o formato criado.
!
! Modificado:
!   26 de maio de 2024
!
! Autoria:
!   oap
! 
SUBROUTINE escrever_massas (self, massas)

  IMPLICIT NONE
  CLASS(arquivo), INTENT(IN) :: self
  REAL(pf), INTENT(IN)       :: massas(:)

  ! salva
  WRITE (self % idarqdata, self % formatoMassas) massas

END SUBROUTINE escrever_massas


! ************************************************************
!! Escrita do cabecalho
!
! Objetivos:
!   Salva no comeco do arquivo algumas informacoes da simulacao,
!   como tamanho do passo, valor de G, etc.
!
! Modificado:
!   26 de maio de 2024
!
! Autoria:
!   oap
!
SUBROUTINE escrever_cabecalho (self, h, G, massas)

  IMPLICIT NONE
  CLASS(arquivo), INTENT(IN) :: self
  REAL(pf), INTENT(IN)       :: massas(:)
  REAL(pf)                   :: h, G

  ! Salva h
  WRITE (self % idarqdata, "(F25.7, :, ',')") h

  ! Salva G
  WRITE (self % idarqdata, "(F25.7, :, ',')") G

  ! Salva as massas
  CALL self%escrever_massas(massas)

END SUBROUTINE escrever_cabecalho

! ************************************************************
!! Escrita do arquivo
!
! Objetivos:
!   Escreve o array no arquivo, conforme formato.
!
! Modificado:
!   26 de maio de 2024
!
! Autoria:
!   oap
!
SUBROUTINE escrever (self, array)

  IMPLICIT NONE
  CLASS(arquivo), INTENT(IN) :: self
  REAL(pf), INTENT(IN)       :: array(2,self % qntdCorpos_int,self % dimensao_int)

  ! salva 
  WRITE (self % idarqdata, self % formato) array

END SUBROUTINE escrever

! ************************************************************
!! Fechamento dos arquivos
!
! Objetivos:
!   Fecha os arquivos.
!
! Modificado:
!   26 de maio de 2024
!
! Autoria:
!   oap
!
SUBROUTINE fechar (self)

  IMPLICIT NONE
  CLASS(arquivo), INTENT(IN) :: self
  
  CLOSE(self % idarqdata)
  CLOSE(self % idarqinfo)

END SUBROUTINE fechar

! ************************************************************
!! Criacao do arquivo de informacoes
!
! Objetivos:
!   Cria o arquivo de informacoes e preenche com o principal
!
! Modificado:
!   25 de maio de 2024
!
! Autoria:
!   oap
!
SUBROUTINE inicializar_arquivo_info (self, N, metodo, G, h, cor, col)

  IMPLICIT NONE
  CLASS(arquivo), INTENT(IN)   :: self
  INTEGER, INTENT(IN)          :: N
  REAL(pf), INTENT(IN)         :: G, h
  CHARACTER(len=*), INTENT(IN) :: metodo
  LOGICAL                      :: cor, col ! correcao, colisao

  WRITE (self % idarqinfo, '(*(g0,1x))') "# gravidade-fortran v0.1"
  WRITE (self % idarqinfo, *) 

  WRITE (self % idarqinfo, '(*(g0,1x))') "# configuracoes"
  WRITE (self % idarqinfo, '(*(g0,1x))') "-- corpos: ", N
  WRITE (self % idarqinfo, '(*(g0,1x))') "-- metodo: ", metodo
  WRITE (self % idarqinfo, '(*(g0,1x))') "-- G: ", G
  WRITE (self % idarqinfo, '(*(g0,1x))') "-- h: ", h
  WRITE (self % idarqinfo, '(*(g0,1x))') "-- passos: " 
  WRITE (self % idarqinfo, '(*(g0,1x))') "-- correcao: ", cor
  WRITE (self % idarqinfo, '(*(g0,1x))') "-- colisoes: ", col

  WRITE (self % idarqinfo, *)

  WRITE (self % idarqinfo, '(*(g0,1x))') "# simulacao: "
  WRITE (self % idarqinfo, '(*(g0,1x))') "-- duracao: "
  CLOSE(self % idarqinfo)

END SUBROUTINE inicializar_arquivo_info

! ************************************************************
!! Atualiza arquivo de informacoes
!
! Objetivos:
!   Atualiza o arquivo de informacoes, alterando a quantidade de
!   passos e a duracao.
!   Como o Fortran nao permite a alteracao de linhas especificas
!   sem apagar o restante do arquivo, esta funcao precisa ler
!   e armazenar o arquivo ja existente para em seguida reescrever.
!   Assim, eh uma funcao lenta, sendo chamada somente no fim
!   da simulacao. A informacao eh previamente guardada no bkp.
!
! Modificado:
!   25 de maio de 2024
!
! Autoria:
!   oap
!
SUBROUTINE atualizar_arquivo_info (self, qntd_passos, duracao)

  IMPLICIT NONE
  CLASS(arquivo), INTENT(INOUT)   :: self
  INTEGER, INTENT(IN)          :: qntd_passos
  REAL(pf), INTENT(IN)         :: duracao
  INTEGER :: linha_qntd_passos = 8, linha_duracao = 13
  INTEGER :: i = 1, tamanho_arquivo = 13, nova_unidade
  LOGICAL :: mudou_qntd_passos = .FALSE., mudou_duracao = .FALSE.
  CHARACTER(len=50), allocatable :: infos(:)

  allocate(infos(tamanho_arquivo))

  CALL capturar_unidade(nova_unidade)
  self % idarqinfo = nova_unidade

  OPEN(newunit=self%idarqinfo,file=self%nomearqinfo,status='old',action='readwrite')
  
  DO i = 1, tamanho_arquivo
    READ (self % idarqinfo, '(A)') infos(i)
  END DO

  ! volta para o comeco do arquivo
  REWIND(self % idarqinfo)

  ! DO WHILE (.NOT. mudou_qntd_passos .AND. .NOT. mudou_duracao)
  DO i = 1, 13
    ! se estiver na linha do tf
    IF (i == linha_qntd_passos) THEN
      WRITE (self % idarqinfo, '(*(g0,1x))') "-- passos: ", qntd_passos
      mudou_qntd_passos = .TRUE.
    ! se estiver na linha da duracao
    ELSE IF (i == linha_duracao) THEN
      WRITE (self % idarqinfo, '(*(g0,1x))') "-- duracao: ", duracao
      mudou_duracao = .TRUE.
    ! se nao
    ELSE
      WRITE (self % idarqinfo, '(*(g0,1x))') infos(i)
    ENDIF
  END DO

END SUBROUTINE atualizar_arquivo_info

! ************************************************************
!! Atualizacao do arquivo de backup
!
! Objetivos:
!   Atualiza o registro no arquivo de backup (bkp).
!
! Modificado:
!   26 de maio de 2024
!
! Autoria:
!   oap
!
SUBROUTINE arquivo_bkp (self, qntd_passos, duracao)
  IMPLICIT NONE
  CLASS(arquivo), INTENT(INOUT)   :: self
  INTEGER, INTENT(IN)          :: qntd_passos
  REAL(pf), INTENT(IN)         :: duracao

  REWIND(self%idarqbkp)
  WRITE(self%idarqbkp, '(*(g0,1x))') qntd_passos, duracao
END SUBROUTINE arquivo_bkp

! ************************************************************
!! Exclusao do arquivo de backup
!
! Objetivos:
!   Com o fim da simulacao nao eh mais necessario ter um arquivo
!   de backup, entao ele eh excluido.
!
! Modificado:
!   26 de maio de 2024
!
! Autoria:
!   oap
!
SUBROUTINE excluir_bkp (self)
  IMPLICIT NONE
  CLASS(arquivo), INTENT(INOUT)   :: self

  CLOSE(self % idarqbkp, status='delete') ! fecha o arquivo
END SUBROUTINE excluir_bkp

! ************************************************************
!! Leitura de arquivo CSV
!
! Objetivos:
!   Le os dados salvos em um arquivo csv no formato deste script.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
SUBROUTINE ler_csv (nome, h, G, massas, R, P)

  CHARACTER(len=*), INTENT(IN)         :: nome
  REAL(pf), INTENT(INOUT)              :: G, h
  REAL(pf), allocatable, INTENT(INOUT) :: R(:,:,:), P(:,:,:), massas(:)
  CHARACTER(len=100000)                :: massas_string
  INTEGER :: iu, i, qntdLinhas = 0, io, qntdCorpos = 0
  REAL(pf) :: t0, tf

  WRITE(*, '(a)') "LER_CSV:"
  WRITE(*, '(a)') "  > arquivo: " // TRIM(nome)

  OPEN(newunit=iu,file=nome,status='old',action='read')

  ! Captura o tamanho do passo
  READ(iu, *) h

  ! Captura a gravidade
  READ(iu, *) G

  ! captura a string de massas
  READ(iu,'(A)') massas_string    
  
  ! captura a quantidade de corpos a partir da quantidade de virgulas
  DO i = 1, len(massas_string)
    IF (massas_string(i:i) == ',') THEN
      qntdCorpos = qntdCorpos + 1
    ENDIF
  END DO
  qntdCorpos = qntdCorpos + 1 ! numero de virgulas = N - 1

  ! captura o numero de linhas do CSV
  DO WHILE (.TRUE.) 
    READ(iu,*, iostat=io)
    IF (io /= 0) exit
    qntdLinhas = qntdLinhas+1
  END DO
  
  REWIND(iu)
  ! Captura o tamanho do passo
  READ(iu, *) h

  ! Captura a gravidade
  READ(iu, *) G    
  
  ! captura as massas
  ALLOCATE(massas(qntdCorpos))
  READ(iu, *) massas

  ! aloca os tamanhos
  ALLOCATE(R(qntdLinhas,qntdCorpos,3))
  ALLOCATE(P(qntdLinhas,qntdCorpos,3))

  ! captura as posicoes e momentos
  t0 = omp_get_wtime()
  DO i = 1, qntdLinhas-1
    READ(iu,*) R(i,:,:),P(i,:,:)
  END DO
  tf = omp_get_wtime()
  
  CLOSE(iu)

  WRITE (*,'(a,F10.4,a)') "  > tempo de leitura: ", tf-t0, "s"
  WRITE (*,*)

END SUBROUTINE ler_csv

! ************************************************************
!! Criacao de diretorio
!
! Objetivos:
!   Cria um diretorio em algum lugar.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
SUBROUTINE criar_dir (dir, onde)

  IMPLICIT NONE
  CHARACTER(LEN=*) :: dir
  CHARACTER(LEN=*),OPTIONAL :: onde
  CHARACTER(LEN=LEN(dir)) :: res
  CHARACTER(:), ALLOCATABLE :: comando
  INTEGER :: i
  res = dir
  ! Remove o "./" se tiver
  DO i = 1, LEN(dir)
    IF (dir(i:i) == "/" .OR. dir(i:i) == ".") THEN
      res(i:i) = " "
    ENDIF
  END DO

  ALLOCATE(CHARACTER(3+LEN(onde)+10+LEN(res)) :: comando)
  comando = "cd "//onde//" && mkdir "// TRIM(res)

  CALL SYSTEM(comando)

  DEALLOCATE(comando)

END SUBROUTINE criar_dir

! ************************************************************
!! Salvamento de preset de sorteio
!
! Objetivos:
!   Escreve um preset sorteado como um preset de valores 
!   iniciais.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
SUBROUTINE salvar_sorteio (onde,subdir,arquivo,nome,G,massas,R,P,t0,tf,timestep,metodo,cor,corme,cormnt,colisoes,pas)

  IMPLICIT NONE
  CHARACTER(LEN=*)      :: onde, subdir, arquivo, metodo, nome
  CHARACTER(LEN=256)    :: dir_arquivo 
  CHARACTER(LEN=3)      :: num_arquivo
  LOGICAL               :: cor, colisoes, diretorio_existe, arquivo_existe
  REAL(pf)              :: G, timestep
  INTEGER               :: t0, tf
  REAL(pf),allocatable  :: massas(:), R(:,:), P(:,:)
  INTEGER               :: pas ! Passos antes de salvar
  INTEGER               :: u, i, arq_i
  REAL(pf)              :: corme ! MARGEM DE ERRO DO CORRETOR
  INTEGER               :: cormnt ! MAXIMO DE NUMERO DE TENTATIVAS DO CORRETOR

  WRITE(*,'(a)') 'SALVAR SORTEIO:'

  ! Verifica se o diretorio desejado existe
  INQUIRE(file=onde//TRIM(subdir), exist=diretorio_existe)
  IF (.NOT. diretorio_existe) THEN
    CALL criar_dir (subdir, onde)
  ENDIF

  ! Agora verifica se o arquivo ja existe
  dir_arquivo = TRIM(onde//subdir) // TRIM(arquivo)
  INQUIRE(file=TRIM(dir_arquivo), exist=arquivo_existe)
  IF (arquivo_existe) THEN
    arq_i = 1
    DO WHILE (arquivo_existe)
      WRITE(num_arquivo, '(I3.3)') arq_i
      INQUIRE(file=TRIM(dir_arquivo)//"_"//TRIM(num_arquivo)//".txt", exist=arquivo_existe)
    END DO
    dir_arquivo = TRIM(dir_arquivo)//"_"//TRIM(num_arquivo)//".txt"
  ENDIF

  WRITE(*,'(a)') '  > arquivo: ' // TRIM(dir_arquivo)

  ! Abre um arquivo
  OPEN(newunit=u,file=dir_arquivo)

  WRITE(u,'(*(g0,1x))') "! Configs"
  WRITE(u,'(*(g0,1x))') "modo vi"
  WRITE(u,'(*(g0,1x))') "nome ", nome
  WRITE(u,'(*(g0,1x))') "integrador ", metodo
  WRITE(u,'(*(g0,1x))') "timestep ", timestep
  WRITE(u,'(*(g0,1x))') "passos ", pas ! Passos antes de salvar
  WRITE(u,'(*(g0,1x))') "t0 ", t0
  WRITE(u,'(*(g0,1x))') "tf ", tf
  WRITE(u,'(*(g0,1x))') "colisoes ", colisoes

  WRITE(u,*)

  WRITE(u,'(*(g0,1x))') "! Opcoes do corretor"
  WRITE(u,'(*(g0,1x))') "corretor ", cor
  WRITE(u,'(*(g0,1x))') "margem_erro ", corme
  WRITE(u,'(*(g0,1x))') "max_num_tentativas ", cormnt

  WRITE(u,*)

  WRITE(u,'(*(g0,1x))') "! Valores do problema"
  WRITE(u,'(*(g0,1x))') "N ", SIZE(massas)
  WRITE(u,'(*(g0,1x))') "G ", G
  
  WRITE(u,*) 

  WRITE(u,'(*(g0,1x))') "! Massas"
  DO i = 1, SIZE(massas)
    WRITE(u,'(*(g0,1x))') massas(i)
  END DO

  WRITE(u,*)

  WRITE(u,'(*(g0,1x))') "! Posicoes"
  DO i = 1, SIZE(massas)
    WRITE(u,'(*(g0,1x,","))') R(i,:) 
  END DO

  WRITE(u,*)

  WRITE(u,'(*(g0,1x))') "! Momentos"
  DO i = 1, SIZE(massas)
    WRITE(u,'(*(g0,1x,","))') P(i,:) 
  END DO

  CLOSE(u)

  WRITE(*,'(a)') '  > arquivo salvo!'
  WRITE(*,*)

END SUBROUTINE salvar_sorteio

! ************************************************************
!! Remocao de espacos vazios
!
! Objetivos:
!   Salva as massas no arquivo de acordo com o formato criado.
!
! Modificado:
!   15 de marco de 2024
!
! Autoria:
!   oap
!
FUNCTION espacosVazios (valor)

  IMPLICIT NONE
  INTEGER, INTENT(IN)           :: valor
  CHARACTER(7)                  :: valor_str
  CHARACTER(:), ALLOCATABLE     :: valor_str_parcial, espacosVazios
  INTEGER                       :: i = 1

  ! transforma o valor em string
  WRITE(valor_str, '(I7)') valor

  ! alinha a esquerda para facilitar
  valor_str = ADJUSTL(valor_str)

  ! onde ficara salvo
  valor_str_parcial = ""
  
  ! elimina os caracteres vazios
  DO WHILE (.TRUE.)
    IF (valor_str(i:i).eq." ") THEN
      i = 1
      exit
    ELSE
      valor_str_parcial = valor_str_parcial // valor_str(i:i)
      i = i + 1
    ENDIF     
  END DO

  ! aloca a string para poder salvar
  ALLOCATE( CHARACTER(LEN_TRIM(valor_str_parcial)) :: espacosVazios)
  ! enfim, salva
  espacosVazios = TRIM(valor_str_parcial)

END FUNCTION espacosVazios

!*****************************************************************************
!! Retorna uma unidade FORTRAN que esteja livre
!
!  Objetivos:
!
!    Uma unidade de FORTRAN "livre" eh um inteiro entre 1 e 99 que nao esta
!    associado a nenhum dispositivo I/O, e eh utilizado para abrir arquivos.
!    Se a unidade eh nula, entao nao ha nenhuma unidade FORTRAN livre.
!    
!    Os numeros 5, 6 e 9 sao reservados, entao nunca sao retornados.
!    
!    O codigo foi baseado na biblioteca GNUFOR de John Burkardt.
!
!  Modificado:
!
!    02 de fevereiro de 2024
!
!  Autoria:
!
!    oap
!
SUBROUTINE capturar_unidade ( iunit )
  IMPLICIT NONE
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ios
  INTEGER ( kind = 4 ) iunit
  LOGICAL lopen
  iunit = 0
  DO i = 1, 99
    IF (i /= 5 .and. i /= 6 .and. i /= 9) THEN
      INQUIRE ( unit = i, opened = lopen, iostat = ios )
      IF ( ios == 0 ) THEN
        IF ( .not. lopen ) THEN
          iunit = i
          RETURN
        ENDIF
      ENDIF
    ENDIF
  END DO
  RETURN
END SUBROUTINE capturar_unidade
END MODULE arquivos