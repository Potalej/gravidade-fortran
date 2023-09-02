! Funcoes auxiliares gerais
! 
! Funcoes
! 
! = function produto_vetorial (u,v)
! Calcula o produto vetorial entre dois vetores u, v em R^3
! 
! = real function determinante (M)
! Calcula o determinante de uma matriz 3x3 dada.
! 
! = function tensor_inercia (m, R)
! Calcula o tensor de inercia dadas a massas e posicao de
! uma particula.
! 
! = function tensor_inercia_geral (massas, posicoes)
! Calcula o tensor de inercia geral dadas as massas e posicoes
! das particulas, sendo este a soma dos tensores de inercia
! individuais.
! 
! = function sistema_linear3 (A,b)
! Resolve um sistema de tres equacoes lineares na forma Ax = b
! usando o metodo de Cramer.
! 
! = function centro_massas (massas, posicoes)
! Calcula o centro de massas do sistema dadas as massas e posicoes
! das particulas.
! 
! = function momentoLinear_total (momentos)
! Calcula o momento linear total do sistema dados os momentos
! lineares individuais das particulas.
! 

module auxiliares
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  implicit none
  private
  public produto_vetorial, tensor_inercia_geral, sistema_linear3, centro_massas, momentoLinear_total

contains

  ! produto vetorial
  function produto_vetorial (u, v)

    implicit none
    real(pf), dimension(:), intent(in) :: u, v
    real(pf), dimension(3)             :: produto_vetorial

    produto_vetorial(:) = 0.0_pf

    produto_vetorial(1) =  u(2)*v(3)-v(2)*u(3)
    produto_vetorial(2) = -u(1)*v(3)+v(1)*u(3)
    produto_vetorial(3) =  u(1)*v(2)-v(1)*u(2)
  
  end function produto_vetorial

  ! determinante de uma matriz 3x3
  real function determinante (M) result (det)

    real(pf), dimension(3,3), intent(in) :: M
    det = M(1,1)*(-M(3,2)*M(2,3)) - M(1,2)*(M(2,1)*M(3,3)-M(2,3)*M(3,1)) + M(1,3)*(M(2,1)*M(3,2)-M(2,2)*M(3,1))
    det = det + M(1,1)*M(2,2)*M(3,3)

  end function determinante

  ! tensor de inercia
  function tensor_inercia (m, R)

    implicit none
    real(pf), dimension(3), intent(in) :: R
    real(pf), intent(in)               :: m
    real(pf), dimension(3,3) :: tensor_inercia
    integer              :: a, b

    do a = 1, 3
      do b = 1, 3
        if (a /= b) then
          tensor_inercia(a,b) = m * R(a) * R(b)
        end if
      end do 
    end do

    tensor_inercia(1,1) = - m * (R(2)**2 + R(3)**2)
    tensor_inercia(2,2) = - m * (R(1)**2 + R(3)**2)
    tensor_inercia(3,3) = - m * (R(1)**2 + R(2)**2)

  end function tensor_inercia

  function tensor_inercia_geral (massas, posicoes)

    implicit none
    real(pf), dimension(:), intent(in) :: massas
    integer :: a
    real(pf), dimension(size(massas),3) :: posicoes
    real(pf), dimension(3,3) :: tensor_inercia_geral

    tensor_inercia_geral(:,:) = 0.0_pf
    
    do a = 1, size(massas)
      tensor_inercia_geral = tensor_inercia_geral + tensor_inercia(massas(a), posicoes(a,:))
    end do   

  end function tensor_inercia_geral

  ! resolve um sistema de 3 equacoes lineares usando o metodo de Cramer
  function sistema_linear3 (A, b)

    implicit none
    real(pf), dimension(3,3), intent(inout) :: A
    real(pf), dimension(3), intent(in)   :: b
    real(pf), dimension(3)               :: sistema_linear3
    real(pf)                             :: detA, c=10**3

    A = A*(1/c)

    ! calcula a determinante de A
    detA = determinante(A) * (c**3)
    A = A*c

    ! aplica o metodo
    sistema_linear3(1) = determinante((/b, A(2,:), A(3,:)/))/detA
    sistema_linear3(2) = determinante((/A(1,:), b , A(3,:)/))/detA
    sistema_linear3(3) = determinante((/A(1,:), A(2,:), b/))/detA

  end function

  ! centro de massas
  function centro_massas (massas, posicoes)

    implicit none
    real(pf), intent(in)   :: massas(:), posicoes(:,:)
    real(pf), dimension(3) :: centro_massas
    integer            :: a

    centro_massas(:) = 0.0_pf

    do a = 1, size(massas)
      centro_massas = centro_massas + massas(a) * posicoes(a,:)
    end do
    centro_massas = centro_massas / sum(massas)

  end function

  ! momento linear total
  function momentoLinear_total (momentos)

    implicit none
    real(pf), intent(in)   :: momentos(:,:)
    real(pf), dimension(3) :: momentoLinear_total
    integer            :: a

    do a = 1, size(momentos,1)
      momentoLinear_total = momentoLinear_total + momentos(a,:)
    end do

  end function momentoLinear_total

end module auxiliares