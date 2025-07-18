! ************************************************************
!! Conexao com o servidor
!
! Modificado:
!   11 de julho de 2025
!
! Autoria:
!   oap
! 
MODULE conexao
    USE iso_c_binding
    USE server_mod
    USE tipos

    PUBLIC conexao_socket
    PRIVATE

    TYPE :: conexao_socket
        INTEGER :: N
        INTEGER :: sockfd

        CONTAINS
            PROCEDURE :: inicializar_plot_tempo_real, &
                        encerrar_conexao, &
                        enviar
    END TYPE

CONTAINS

SUBROUTINE inicializar_plot_tempo_real (self, N)
    CLASS (conexao_socket) :: self
    INTEGER, INTENT(IN) :: N
    LOGICAL :: status

    self % N = N

    ! Criando um socket
    self % sockfd = criar_socket()

    ! Criando um servidor e conectando
    CALL conectar_servidor(self % sockfd, status)

    ! Se o servidor nao tiver dado certo, avisa e encerra o programa
    IF (.NOT. status) THEN
        WRITE (*,'(A)') "[SOCKET] Nao foi possivel conectar ao servidor."
        WRITE (*,'(A)') "[SOCKET] Verifique se o servidor esta ligado e tente novamente."
        STOP 0
    ENDIF

    ! Se tiver dado certo, precisa enviar o cabecalho (valor de N)
    CALL enviar_cabecalho(self % sockfd, self % N)
END SUBROUTINE

SUBROUTINE encerrar_conexao (self)
    CLASS (conexao_socket) :: self
END SUBROUTINE

SUBROUTINE enviar (self, R, P)
    CLASS (conexao_socket) :: self
    REAL(pf) :: R(self % N,3), P(self % N,3)
    REAL(c_double) :: pos(self % N,3), mom(self % N,3)

    pos = REAL(R, kind=c_double)
    mom = REAL(P, kind=c_double)

    CALL enviar_dados(self % sockfd, self % N, pos)
END SUBROUTINE

END MODULE conexao