! Funções auxiliares gerais
! 
! Funções
! 
! = produto_vetorial (u,v)
!   Clacula o produto vetorial entre dois vetores u, v \in R^3

module auxiliares

  implicit none
  private
  public produto_vetorial, tensor_inercia_geral, sistema_linear3

contains

  ! produto vetorial
  function produto_vetorial (u, v)

    implicit none
    real, dimension(3), intent(in) :: u, v
    real, dimension(3)             :: produto_vetorial

    produto_vetorial(:) = 0

    produto_vetorial(1) =  u(2)*v(3)-v(2)*u(3)
    produto_vetorial(2) = -u(1)*v(3)+v(1)*u(3)
    produto_vetorial(3) =  u(1)*v(2)-v(1)*u(2)
  
  end function produto_vetorial

  ! determinante de uma matriz 3x3
  real function determinante (M) result (det)

    real, dimension(3,3), intent(in) :: M

    det = M(1,1)*(M(2,2)*M(3,3)-M(3,2)*M(2,3)) - M(1,2)*(M(2,1)*M(3,3)-M(2,3)*M(3,1)) + M(1,3)*(M(2,1)*M(3,2)-M(2,2)*M(3,1))

  end function determinante

  ! tensor de inércia
  function tensor_inercia (m, R)

    implicit none
    real, dimension(3), intent(in) :: R
    real, intent(in)               :: m
    real, dimension(3,3) :: tensor_inercia
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
    real, dimension(:), intent(in) :: massas
    integer :: a
    real, dimension(size(massas),3) :: posicoes
    real, dimension(3,3) :: tensor_inercia_geral

    tensor_inercia_geral(:,:) = 0
    
    do a = 1, size(massas)
      tensor_inercia_geral = tensor_inercia_geral + tensor_inercia(massas(a), posicoes(a,:))
    end do   

  end function tensor_inercia_geral

  ! resolve um sistema de 3 equações lineares usando o método de Cramer
  function sistema_linear3 (A, b)

    implicit none
    real, dimension(3,3), intent(in) :: A
    real, dimension(3), intent(in)   :: b
    real, dimension(3)               :: sistema_linear3
    real                             :: detA

    ! calcula a determinante de A
    detA = determinante(A)

    ! aplica o método
    sistema_linear3(1) = (1/detA) * determinante((/b, A(2,:), A(3,:)/))
    sistema_linear3(2) = (1/detA) * determinante((/A(1,:), b , A(3,:)/))
    sistema_linear3(3) = (1/detA) * determinante((/A(1,:), A(2,:), b/))

  end function

end module auxiliares