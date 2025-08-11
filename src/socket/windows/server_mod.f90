! ************************************************************
!! Sockets e servidor
!
! Modificado:
!   11 de agosto de 2025
!
! Autoria:
!   oap
! 
MODULE server_mod
    USE iso_c_binding
    USE wsa_mod
    IMPLICIT NONE

    INTERFACE
        FUNCTION socket(domain, TYPE, protocol) bind(C)
            IMPORT :: c_int
            INTEGER(c_int) :: socket
            INTEGER(c_int), value :: domain, TYPE, protocol
        END FUNCTION

        FUNCTION connect(sockfd, addr, addrlen) bind(C)
            IMPORT :: c_int, c_ptr
            INTEGER(c_int) :: connect
            INTEGER(c_int), value :: sockfd, addrlen
            TYPE(c_ptr), value :: addr
        END FUNCTION

        FUNCTION send(sockfd, buf, len, flags) bind(C)
            IMPORT :: c_int, c_ptr
            INTEGER(c_int) :: send
            INTEGER(c_int), value :: sockfd, len, flags
            TYPE(c_ptr), value :: buf
        END FUNCTION

        FUNCTION htons(port) bind(C)
            IMPORT :: c_int
            INTEGER(c_int), value :: port
            INTEGER(c_int) :: htons
        END FUNCTION

        FUNCTION inet_addr(cp) bind(C)
            IMPORT :: c_char, c_int
            CHARACTER(kind=c_char), dimension(*) :: cp
            INTEGER(c_int) :: inet_addr
        END FUNCTION
    END INTERFACE

    TYPE, bind(C) :: sockaddr_in
        INTEGER(c_short) :: sin_family
        INTEGER(c_short) :: sin_port
        INTEGER(c_int)   :: sin_addr
        CHARACTER(kind=c_char) :: sin_zero(8)
    END TYPE sockaddr_in

contains

    FUNCTION criar_socket () RESULT(sockfd)
        INTEGER(c_int) :: sockfd

        sockfd = socket(2_c_int, 1_c_int, 0_c_int)
        IF (sockfd < 0) THEN
            WRITE(*,'(A)') "[SOCKET] Erro ao criar socket"
            CALL wsa_cleanup()
            STOP
        ENDIF
    END FUNCTION

    SUBROUTINE conectar_servidor (sockfd, status)
        INTEGER(c_int) :: sockfd
        TYPE(sockaddr_in), TARGET :: server_addr
        INTEGER(c_int) :: res_conexao
        LOGICAL :: status

        ! Configurando o servidor
        server_addr%sin_family = 2_c_short
        server_addr%sin_port = htons(50007_c_int)
        server_addr%sin_addr = inet_addr("127.0.0.1"//c_null_char)
        server_addr%sin_zero = c_null_char

        ! Criando a conexao
        res_conexao = connect(sockfd, c_loc(server_addr), 16_c_int)
        IF (res_conexao /= 0) THEN
            WRITE(*,'(A)') "[SOCKET] Erro ao conectar"
            CALL wsa_cleanup()
            status = .FALSE.
        ELSE
            status = .TRUE.
        ENDIF
    END SUBROUTINE

    SUBROUTINE enviar_cabecalho (sockfd, N)
        INTEGER(c_int) :: sockfd
        INTEGER(c_int) :: res_envio
        INTEGER(c_int), INTENT(IN), TARGET :: N

        ! Enviando os dados
        res_envio = send(sockfd, c_loc(N), int(c_sizeof(N), c_int), 0_c_int)
        CALL sleep(1) ! Tempo para o python receber os dados

        IF (res_envio <= 0) THEN
            WRITE(*,*) "[SOCKET] Erro no envio do cabecalho (", res_envio, ")"
            CALL wsa_cleanup()
            STOP 0
        ENDIF
    END SUBROUTINE

    SUBROUTINE enviar_dados (sockfd, t, N, R)
        INTEGER(c_int) :: sockfd
        CHARACTER(len=100), target :: linha
        CHARACTER(len=65536), target :: buffer
        REAL(c_double), target :: t, R(:,:)
        INTEGER(c_int) :: res_envio
        INTEGER :: i, N

        ! Montando os dados
        WRITE(linha, '(F16.4)') t
        buffer = TRIM(linha)//CHAR(10)
        DO i = 1, N
            WRITE(linha, '(F16.4,1X,F16.4,1X,F16.4)') R(i,1), R(i,2), R(i,3)
            buffer = TRIM(buffer)//TRIM(linha)//CHAR(10)
        END DO

        ! Enviando os dados
        res_envio = send(sockfd, c_loc(buffer), len(trim(buffer)), 0_c_int)
    END SUBROUTINE

END MODULE server_mod