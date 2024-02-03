! Hamiltoniano
! 
! Aqui constam funcoes ligadas ao hamiltoniano, ou seja, ao calculo da
! energia do sistema.
! 
! Funcoes
! 
! = energia_cinetica (m, P) result(soma_ec)
!   Calcula a energia cinetica do sistema informados os vetores de massas
!   e momentos lineares.
! 
! = energia_potencial (m, R) result(soma_ep)
!   Calcula a energia potencial do sistema informados os vetores de massas
!   e posicoes.
! 
! = energia_total (m, R, P) result(soma_et)
!   Calcula a energia total H = T - V do sistema, onde T eh a energia cinética
!   e V = -U, onde U eh a energia potencial do sistema.
! 

module hamiltoniano
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  implicit none
  private
  public energia_cinetica, energia_potencial, energia_total, energia_total_aproximada

contains

  ! Energia cinetica
  function energia_cinetica (m, P)
    implicit none
    real(pf) :: m(:), P(:,:), energia_cinetica
    integer :: a
    
    energia_cinetica = 0
    do a = 1, size(m)
      energia_cinetica = energia_cinetica + norm2(P(a,:))**2 / (2 * m(a))
    end do 
  end function

  ! Energia potencial
  function energia_potencial (G, m, R)
    implicit none
    real(pf) :: m(:), R(:,:)
    real(pf) :: distancia, energia_potencial, G
    integer :: a = 1, b = 1

    energia_potencial = 0.0_pf
    do a = 2, size(m)
      do b = 1, a-1
        distancia = norm2(R(a,:) - R(b,:))
        energia_potencial = energia_potencial + m(a) * m(b) / distancia
      end do
    end do
    
    energia_potencial = - G*energia_potencial
  end function

  ! Energia total
  function energia_total (G, m, R, P)
    implicit none
    real(pf) :: m(:), R(:,:), P(:,:)
    real(pf) :: energia_total, G
    energia_total = energia_cinetica(m, P) + energia_potencial(G, m, R)
  end function

  ! Energia total aproximada
  ! Baseado na relacao de Lagrange-Jacobi, a segunda derivada do momento de inercia
  ! fornece uma aproximacao para a diferenca entre a energia total e a potencial.
  ! Com alguns rearranjos eh possivel estimar a energia total.
  function energia_total_aproximada (qntd_passos, h, G, m, R_ant, P_ant, R_pos, P_pos)
    implicit none
    real(pf) :: m(:), R_ant(:,:), R_pos(:,:), P_ant(:,:), P_pos(:,:)
    real(pf) :: G, h
    real(pf) :: D_ant = 0.0_pf, D_pos = 0.0_pf, T = 0.0_pf, der_D
    real(pf) :: energia_total_aproximada
    integer   :: a, qntd_passos

    ! Derivada do momento de dilatacao D = sum ri dot pi
    do a = 1, size(m)
      D_ant = D_ant + DOT_PRODUCT(R_ant(a,:), P_ant(a,:))
      D_pos = D_pos + DOT_PRODUCT(R_pos(a,:), P_pos(a,:))
    end do

    ! energia cinetica
    T = energia_cinetica(m, P_pos)

    ! Derivada de D
    der_D = (D_pos - D_ant)/(h*qntd_passos)

    WRITE (*,*) 'D_pos=', D_pos, '| D_ant=', D_ant, '| der_D=', der_D, '| T=', T

    energia_total_aproximada = der_D - T
  end function

end module hamiltoniano