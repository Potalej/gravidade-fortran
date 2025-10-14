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
!   Utiliza o JSON-Fortran.
!   
! Modificado:
!   08 de agosto de 2025
! 
! Autoria:
!   oap
! 
MODULE arquivos
  USE tipos
  USE diretorio
  USE string_utils
  USE json_utils_mod

  IMPLICIT NONE
  PUBLIC arquivo, ler_csv, capturar_unidade

  TYPE :: arquivo
    !> Ids dos arquivos data, info e backup
    INTEGER :: id_arq_data, id_arq_info, id_arq_bkp
    
    !> Informacoes sobre a simulacao
    INTEGER :: qntd_corpos_int, dimensao_int

    !> Nomes dos arquivos
    CHARACTER(:), ALLOCATABLE :: dir_arq, nome_arq_data, nome_arq_info, nome_arq_bkp
    
    !> Formatos
    CHARACTER(:), ALLOCATABLE :: formato, formato_massas, qntd_corpos, dimensao

    !> Diretorios padrao
    CHARACTER(:), ALLOCATABLE :: dir_out

    CONTAINS
      !> Rotinas e funcoes
      PROCEDURE :: definir_diretorio_saida,& ! Define o diretorio de saida
                   gerar_nome_diretorio,& ! Gera o nome do diretorio de saida
                   criar_formatos,&       ! Cria os formatos para massas e dados
                   criar_data,&           ! Cria is arquivos de dados de saida
                   escrever_data,&        ! Escreve arrays no data.csv
                   fechar,&               ! Fecha os arquivos de saida
                   excluir_bkp,&          ! Exclui o arquivo de backup
                   escrever_cabecalho_data,&  ! Escreve o cabecalho do data.csv
                   inicializar_arquivo_info,& ! Inicializa arquivo de informacoes
                   atualizar_arquivo_info,&   ! Atualiza o arquivo de informacoes
                   atualizar_arquivo_bkp      ! Atualiza o arquivo de backup
  END TYPE arquivo

CONTAINS

SUBROUTINE definir_diretorio_saida (self, dir_param)
  CLASS(arquivo), INTENT(INOUT) :: self
  CHARACTER(LEN=*), INTENT(INOUT), OPTIONAL :: dir_param
  CHARACTER(:), ALLOCATABLE :: dir
  
  IF (PRESENT(dir_param)) THEN
    ALLOCATE(CHARACTER(LEN(TRIM(dir_param))) :: dir)
    dir = dir_param
  ELSE
    ALLOCATE(CHARACTER(5) :: dir)
    dir = "./out"
  ENDIF
  
  IF (ALLOCATED(self % dir_out)) DEALLOCATE(self % dir_out)
  ALLOCATE(CHARACTER(LEN(TRIM(dir))) :: self % dir_out)
  self % dir_out = TRIM(dir)
END SUBROUTINE definir_diretorio_saida

! ************************************************************
!! Nome da pasta de saida
!
! Objetivos:
!   Cria o nome do diretorio baseado no dia corrente e contando
!   a partir de 1, ou seja, nao havendo pastas do mesmo dia
!   eh criado um AAAAMMDD_01, e caso haja eh criado AAAAMMDD_02, 
!   AAAAMMDD_03, etc.
!
! Modificado:
!   26 de maio de 2024 (criado)
!   10 de outubro de 2025 (modificado)
!
! Autoria:
!   oap
! 
SUBROUTINE gerar_nome_diretorio (self)
  CLASS(arquivo), INTENT(INOUT) :: self
  CHARACTER(3) :: numero
  INTEGER :: i
  LOGICAL :: existe
  
  ! Para capturar a data de hoje
  CHARACTER(8) :: data_hoje

  CALL DATE_AND_TIME(data_hoje)

  ! Verifica se existe o diretorio de saida
  IF (.NOT. ALLOCATED(self % dir_out)) CALL self % definir_diretorio_saida()
  CALL diretorio_data(self % dir_out)

  ! Descobrindo iterativamente qual o nome do arquivo
  existe = .TRUE.
  i = 1
  DO WHILE (existe)
    WRITE(numero, '(I3.3)') i
    i = i + 1

    ! Cria nomes
    self % dir_arq = TRIM(data_hoje)//"_"//TRIM(numero)
    self % nome_arq_data = self%dir_out // "/data/" // self%dir_arq // "/data.bin"
    self % nome_arq_info = self%dir_out // "/data/" // self%dir_arq // "/info.txt"
    self % nome_arq_bkp  = self%dir_out // "/data/" // self%dir_arq // "/bkp.txt"
    
    ! Verifica se a pasta existe
    INQUIRE(file=TRIM(self%dir_out // "/data/" // self%dir_arq), exist=existe)
  END DO

END SUBROUTINE gerar_nome_diretorio

! ************************************************************
!! Formato do arquivo
!
! Objetivos:
!   Para uma determinada quantidade de particulas e dimensoes,
!   eh criada a formatacao para transformar o array em uma 
!   string corretamente. SENDo N:= qntd_corpos e D := dimensao,
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
SUBROUTINE criar_formatos (self, qntd_corpos, dimensao)

  IMPLICIT NONE
  CLASS(arquivo), INTENT(INOUT) :: self
  INTEGER, INTENT(IN)           :: qntd_corpos, dimensao

  ! salva a quantidade de corpos e dimensao
  self % qntd_corpos = int_para_string(qntd_corpos)
  self % dimensao = int_para_string(dimensao)

  self % qntd_corpos_int = qntd_corpos
  self % dimensao_int = dimensao

  self % formato = '(2(' // self % dimensao // '(' // self % qntd_corpos // '(F25.13, :, ","))))'
  self % formato_massas = '(' // self % qntd_corpos // '(F25.7, :, ","))'

END SUBROUTINE criar_formatos

! ************************************************************
!! Criacao dos arquivos de dados de saida
!
! Objetivos:
!   Cria os arquivos de saida da simulacao dentro do diretorio
!   devido.
!
! Modificado:
!   10 de outubro de 2025
!
! Autoria:
!   oap
! 
SUBROUTINE criar_data (self, qntd_corpos, dimensao)

  CLASS(arquivo), INTENT(INOUT) :: self
  INTEGER, INTENT(IN)           :: qntd_corpos, dimensao
  INTEGER(kind=4)               :: id_arq_data, id_arq_info, id_arq_bkp

  WRITE (*, '(a)') 'CRIAR ARQUIVO PARA SALVAR PLOT:'

  ! cria formatacao
  CALL self % criar_formatos(qntd_corpos, dimensao)
  WRITE (*, '(a)') '  > formato : ' // self % formato

  ! criacao do nome do diretorio
  CALL self % gerar_nome_diretorio()
  WRITE (*,'(a)') '  > diretorio de saida: ' // self % dir_arq

  ! cria o diretorio
  CALL criar_dir(self % dir_arq, self % dir_out // "/data")

  ! cria o arquivo "data"
  CALL capturar_unidade(id_arq_data)
  self % id_arq_data = id_arq_data
  OPEN(id_arq_data, file = self % nome_arq_data, status='new', access='stream')

  ! cria o arquivo info
  CALL capturar_unidade(id_arq_info)
  self % id_arq_info = id_arq_info
  OPEN(id_arq_info, file = self % nome_arq_info, status='new')

  ! cria o arquivo bkp
  CALL capturar_unidade(id_arq_bkp)
  self % id_arq_bkp = id_arq_bkp
  OPEN(id_arq_bkp, file = self % nome_arq_bkp, status='new')

  WRITE (*,'(a)') '  > arquivos criados!'
  WRITE (*,*)

END SUBROUTINE criar_data

! ************************************************************
!! Escrita do cabecalho no arquivo "data"
!
! Objetivos:
!   Salva no comeco do arquivo algumas informacoes da simulacao,
!   como tamanho do passo, valor de G, etc.
!
! Modificado:
!   10 de outubro de 2025
!
! Autoria:
!   oap
!
SUBROUTINE escrever_cabecalho_data (self, h, G, massas)

  IMPLICIT NONE
  CLASS(arquivo), INTENT(IN) :: self
  REAL(pf), INTENT(IN)       :: massas(:)
  REAL(pf)                   :: h, G

  ! Salva h, G, N
  ! WRITE (self % id_arq_data, "(F25.7, :, ',')") h
  WRITE (self % id_arq_data) h, G, SIZE(massas)

  ! Salva as massas
  WRITE (self % id_arq_data) massas

END SUBROUTINE escrever_cabecalho_data

! ************************************************************
!! Escrita no arquivo "data"
!
! Objetivos:
!   Escreve um array no arquivo "data", conforme formato.
!
! Modificado:
!   10 de outubro de 2025
!
! Autoria:
!   oap
!
SUBROUTINE escrever_data (self, array)

  IMPLICIT NONE
  CLASS(arquivo), INTENT(IN) :: self
  REAL(pf), INTENT(IN)       :: array(2,self%qntd_corpos_int,self%dimensao_int)

  ! salva 
  WRITE (self % id_arq_data) array

END SUBROUTINE escrever_data

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
  
  CLOSE(self % id_arq_data)
  CLOSE(self % id_arq_info)

END SUBROUTINE fechar

! ************************************************************
!! Criacao do arquivo de informacoes
!
! Objetivos:
!   Cria o arquivo de informacoes e preenche com o principal
!
! Modificado:
!   24 de julho de 2025
!
! Autoria:
!   oap
!
SUBROUTINE inicializar_arquivo_info (self, infos, version_string, precisao)
  CLASS(arquivo), INTENT(IN)   :: self
  TYPE(json_value), POINTER    :: infos
  INTEGER          :: N, t0, tf, checkpoints ! checkpoints
  REAL(pf)         :: G, h, soft
  LOGICAL          :: corrigir, colidir ! correcao
  CHARACTER(len=:), ALLOCATABLE :: metodo, colidir_modo ! colisao
  REAL(pf)         :: corme, colmd
  INTEGER          :: cormnt
  LOGICAL          :: paralelo, gpu ! forcas paralelas, gpu
  CHARACTER(20)    :: data_hora_str
  CHARACTER(LEN=*) :: version_string, precisao

  REAL(pf) :: densidade
  LOGICAL :: encontrado

  CALL json % get(infos, "N", N)
  G = json_get_float(infos, "G")
  
  ! Integracao
  metodo = json_get_string(infos, 'integracao.metodo')
  CALL json % get(infos, 'integracao.t0', t0)
  CALL json % get(infos, 'integracao.tf', tf)
  h = json_get_float(infos, "integracao.timestep")
  soft = json_get_float(infos, "integracao.amortecedor")
  CALL json % get(infos, "integracao.checkpoints", checkpoints)

  ! Correcao
  CALL json % get(infos, 'correcao.corrigir', corrigir)
  corme = json_get_float(infos, 'correcao.margem_erro')
  CALL json % get(infos, 'correcao.max_num_tentativas', cormnt)

  ! Colisao
  CALL json % get(infos, 'colisoes.colidir', colidir)
  colidir_modo = json_get_string(infos, 'colisoes.metodo')
  densidade = json_get_float(infos, 'colisoes.densidade')

  ! Paralelizacao
  CALL json % get(infos, 'paralelo', paralelo)
  CALL json % get(infos, 'gpu', gpu, encontrado)
  IF (.NOT. encontrado) gpu = .FALSE.

  ! Data e hora de inicio
  CALL data_hora_string(data_hora_str)

  WRITE (self % id_arq_info, '(*(g0,1x))') "# gravidade-fortran v"//version_string//'_'//precisao
  WRITE (self % id_arq_info, *) 

  WRITE (self % id_arq_info, '(*(g0,1x))') "# configuracoes"
  WRITE (self % id_arq_info, '(*(g0,1x))') "-- corpos: ", N
  WRITE (self % id_arq_info, '(*(g0,1x))') "-- metodo: ", metodo
  WRITE (self % id_arq_info, '(*(g0,1x))') "-- G: ", G
  WRITE (self % id_arq_info, '(*(g0,1x))') "-- h: ", h
  WRITE (self % id_arq_info, '(*(g0,1x))') "-- amortecimento: ", soft
  WRITE (self % id_arq_info, '(*(g0,1x))') "-- checkpoints: ", checkpoints
  WRITE (self % id_arq_info, '(*(g0,1x))') "-- total passos: " 
  WRITE (self % id_arq_info, '(*(g0,1x))') "-- t0: ", t0
  WRITE (self % id_arq_info, '(*(g0,1x))') "-- tf: ", tf
  WRITE (self % id_arq_info, '(*(g0,1x,1x))') "-- paralelizacao: ", paralelo, gpu
  
  IF (corrigir) THEN
    WRITE (self % id_arq_info, '(*(g0,1x))') "-- correcao: ", corrigir
    WRITE (self % id_arq_info, '(*(g0,1x))') "-- correcao margem erro: ", corme
    WRITE (self % id_arq_info, '(*(g0,1x))') "-- correcao max num tent.: ", cormnt
  ELSE
    WRITE (self % id_arq_info, '(*(g0,1x))') "-- correcao: ", corrigir
    WRITE (self % id_arq_info, '(*(g0,1x))') "-- correcao margem erro: n/a"
    WRITE (self % id_arq_info, '(*(g0,1x))') "-- correcao max num tent.: n/a"
  ENDIF

  IF (colidir) THEN
    WRITE (self % id_arq_info, '(*(g0,1x))') "-- colisoes: ", colidir, colidir_modo
    WRITE (self % id_arq_info, '(*(g0,1x))') "-- densidade: ", densidade
  ELSE
    WRITE (self % id_arq_info, '(*(g0,1x))') "-- colisoes: ", colidir
    WRITE (self % id_arq_info, '(*(g0,1x))') "-- densidade: n/a"
  ENDIF

  WRITE (self % id_arq_info, *)

  WRITE (self % id_arq_info, '(*(g0,1x))') "# simulacao: "
  WRITE (self % id_arq_info, '(*(g0,1x))') "-- inicio: ", data_hora_str
  WRITE (self % id_arq_info, '(*(g0,1x))') "-- duracao: "

  CLOSE(self % id_arq_info)
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
!   24 de julho de 2025
!
! Autoria:
!   oap
!
SUBROUTINE atualizar_arquivo_info (self, qntd_passos, duracao)
  CLASS(arquivo), INTENT(INOUT)   :: self
  INTEGER, INTENT(IN)          :: qntd_passos
  REAL(pf64), INTENT(IN)         :: duracao
  INTEGER :: linha_qntd_passos = 10, linha_duracao = 22
  INTEGER :: i = 1, tamanho_arquivo = 22, nova_unidade
  LOGICAL :: mudou_qntd_passos = .FALSE., mudou_duracao = .FALSE.
  CHARACTER(len=50), allocatable :: infos(:)

  allocate(infos(tamanho_arquivo))

  CALL capturar_unidade(nova_unidade)
  self % id_arq_info = nova_unidade

  OPEN(newunit=self%id_arq_info,file=self%nome_arq_info,status='old',action='readwrite')
  
  DO i = 1, tamanho_arquivo
    READ (self % id_arq_info, '(A)') infos(i)
  END DO

  ! volta para o comeco do arquivo
  REWIND(self % id_arq_info)

  ! DO WHILE (.NOT. mudou_qntd_passos .AND. .NOT. mudou_duracao)
  DO i = 1, 22
    ! se estiver na linha do tf
    IF (i == linha_qntd_passos) THEN
      WRITE (self % id_arq_info, '(*(g0,1x))') "-- passos: ", qntd_passos
      mudou_qntd_passos = .TRUE.
    ! se estiver na linha da duracao
    ELSE IF (i == linha_duracao) THEN
      WRITE (self % id_arq_info, '(*(g0,1x))') "-- duracao: ", duracao
      mudou_duracao = .TRUE.
    ! se nao
    ELSE
      WRITE (self % id_arq_info, '(*(g0,1x))') infos(i)
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
SUBROUTINE atualizar_arquivo_bkp (self, qntd_passos, duracao)
  CLASS(arquivo), INTENT(INOUT)   :: self
  INTEGER, INTENT(IN)             :: qntd_passos
  REAL(pf64), INTENT(IN)          :: duracao

  REWIND(self%id_arq_bkp)
  WRITE(self%id_arq_bkp, '(*(g0,1x))') qntd_passos, duracao
END SUBROUTINE atualizar_arquivo_bkp

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

  CLOSE(self % id_arq_bkp, status='delete') ! fecha o arquivo
END SUBROUTINE excluir_bkp


!*****************************************************************************
!! Retorna uma unidade FORTRAN que esteja livre
!
!  Objetivos:
!    Uma unidade de FORTRAN "livre" eh um inteiro entre 1 e 99 que nao esta
!    associado a nenhum dispositivo I/O, e eh utilizado para abrir arquivos.
!    Se a unidade eh nula, entao nao ha nenhuma unidade FORTRAN livre.
!    
!    Os numeros 5, 6 e 9 sao reservados, entao nunca sao retornados.
!    
!    O codigo foi baseado na biblioteca GNUFOR de John Burkardt.
!
!  Modificado:
!    02 de fevereiro de 2024
!
!  Autoria:
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
  INTEGER :: t_rate, t0, tf

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
  CALL SYSTEM_CLOCK(count_rate=t_rate)
  CALL SYSTEM_CLOCK(t0)
  DO i = 1, qntdLinhas-1
    READ(iu,*) R(i,:,:),P(i,:,:)
  END DO
  CALL SYSTEM_CLOCK(tf)
  
  CLOSE(iu)

  WRITE (*,'(a,F10.4,a)') "  > tempo de leitura: ", REAL(tf-t0)/REAL(t_rate), "s"
  WRITE (*,*)

END SUBROUTINE ler_csv

END MODULE arquivos