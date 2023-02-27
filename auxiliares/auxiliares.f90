! Funções auxiliares gerais
! 
! Funções
! 
! = produto_vetorial (u,v)
!   Clacula o produto vetorial entre dois vetores u, v \in R^3

module auxiliares

  implicit none
  private
  public produto_vetorial

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

end module auxiliares