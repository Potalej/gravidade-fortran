! Hamiltoniano
! 
! Aqui constam funções ligadas ao hamiltoniano, ou seja, ao cálculo da
! energia do sistema.
! 
! Funções
! 
! = energia_cinetica (m, P) result(soma_ec)
!   Calcula a energia cinética do sistema informados os vetores de massas
!   e momentos lineares.
! 
! = energia_potencial (m, R) result(soma_ep)
!   Calcula a energia potencial do sistema informados os vetores de massas
!   e posições.
! 
! = energia_total (m, R, P) result(soma_et)
!   Calcula a energia total H = T - V do sistema, onde T é a energia cinética
!   e V = -U, onde U é a energia potencial do sistema.

module hamiltoniano

  implicit none
  private
  public energia_cinetica, energia_potencial, energia_total

contains

  ! Energia cinética
  real function energia_cinetica (m, P) result(soma_ec)
    implicit none
    real, allocatable :: m(:), P(:,:)
    integer :: a
    
    soma_ec = 0
    do a = 1, size(m)
      soma_ec = soma_ec + norm2(P(a,:))**2 / (2 * m(a))
    end do 
  end function

  ! Energia potencial
  real function energia_potencial (m, R) result(soma_ep)
    implicit none
    real, allocatable :: m(:), R(:,:)
    real :: distancia
    integer :: a = 1, b = 1

    soma_ep = 0
    do a = 2, size(m)
      do b = 1, a-1
        distancia = norm2(R(a,:) - R(b,:))
        soma_ep = soma_ep + m(a) * m(b) / distancia
      end do
    end do
    
    soma_ep = - soma_ep
  end function

  ! Energia total
  real function energia_total (m, R, P) result(soma)
    implicit none
    real, allocatable :: m(:), R(:,:), P(:,:)
    
    soma = 0
    soma = energia_cinetica(m, P) + energia_potencial(m, R)
  end function

end module hamiltoniano