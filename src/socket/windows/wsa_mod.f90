! ************************************************************
!! Windows Sockets
!
! Modificado:
!   11 de julho de 2025
!
! Autoria:
!   oap
! 
MODULE wsa_mod
    USE iso_c_binding
    IMPLICIT NONE

    INTERFACE
        FUNCTION WSAStartup (wVersionRequested, lpWSAData) bind (C, name="WSAStartup")
            IMPORT :: c_int, c_short, c_ptr
            INTEGER(c_int) :: WSAStartup
            INTEGER(c_short), value :: wVersionRequested
            TYPE(c_ptr), value :: lpWSAData
        END FUNCTION

        FUNCTION WSACleanup () bind (C, name="WSACleanup")
            IMPORT :: c_int
            INTEGER(c_int) :: WSACleanup
        END FUNCTION
    end interface

    TYPE, bind(C) :: WSADATA
        INTEGER(c_short) :: wVersion
        INTEGER(c_short) :: wHighVersion
        CHARACTER(kind=c_char) :: data(512)
    END TYPE WSADATA

CONTAINS
    SUBROUTINE wsa_startup ()
        INTEGER(c_short), parameter :: wsa_version = int(514, c_short)
        TYPE(WSADATA), target :: wsaData_var
        INTEGER(c_int) :: res_ws
        
        ! Limpa qualquer resquicio que possa ter
        CALL wsa_cleanup()

        res_ws = WSAStartup(wsa_version, c_loc(wsaData_var))
        IF (res_ws /= 0) THEN
            WRITE(*,'(A)') "[WSA] Erro no WSAStartup"
            STOP 0
        END IF
    END SUBROUTINE

    SUBROUTINE wsa_cleanup ()
        INTEGER(c_int) :: res_ws
        res_ws = WSACleanup()
    END SUBROUTINE
END MODULE wsa_mod