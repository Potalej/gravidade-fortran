! Angular
! 
! Aqui constam funcoes ligadas ao momento angular do sistema.
! 
! = function angular_individual (Ra, Pa)
! Calcula o momento angular de uma particula dadas as suas coordenadas
! de posicao e momento linear
!
! = function angular_geral (R, P)
! Calcula o momento angular do sistema inteiro, sendo entao a soma do momento
! angular individual de cada particula
!

module angular
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use auxiliares

  implicit none
  private
  public angular_geral

contains

  ! calcula o momento angular de uma partícula
  function angular_individual (Ra, Pa)
    implicit none
    real(pf), dimension(:), intent(in) :: Ra, Pa
    real(pf), dimension(3)             :: angular_individual
    
    angular_individual = produto_vetorial(Ra, Pa)

  end function angular_individual

  ! calcula o momento angular do sistema inteiro
  function angular_geral (R, P)
    implicit none
    real(pf), dimension(:,:), intent(in) :: R, P
    real(pf), dimension(3)               :: angular_geral
    integer                          :: a

    angular_geral = (/0.0_pf,0.0_pf,0.0_pf/)
    do a = 1, size(R,1)
      angular_geral = angular_geral + angular_individual(R(a,:), P(a,:))      
    end do


  end function angular_geral

end module angular