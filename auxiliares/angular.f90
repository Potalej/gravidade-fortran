! Angular
! 
! Aqui constam funções ligadas ao momento angular do sistema.
! 
! Funções
! 

module angular
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use auxiliares

  implicit none
  private
  public angular_geral, angular_correcao

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

  ! ajuste do momento angular
  subroutine angular_correcao (m, R, P, J0)
    implicit none
    real(pf), dimension(:), intent(in)      :: m
    real(pf), dimension(:,:), intent(inout) :: R, P
    real(pf), dimension(3), intent(in)      :: J0
    integer :: a

    ! gradientes e normas
    real(pf), dimension(2, size(R,1), 3) :: gradJx, gradJy, gradJz
    real(pf), dimension(3)                     :: J, difJ, normas2
    real(pf)                                   :: Rx, Ry, Rz, Px, Py, Pz
    
    gradJx(:,:,:) = 0.0_pf
    gradJy(:,:,:) = 0.0_pf
    gradJz(:,:,:) = 0.0_pf

    normas2(:) = 0.0_pf

    do a = 1, size(m)

      Rx = R(a, 1)
      Ry = R(a, 2)
      Rz = R(a, 3)
      
      Px = P(a, 1)
      Py = P(a, 2)
      Pz = P(a, 3)

      gradJx(1,a,:) = (/0.0_pf, Pz, -Py/)
      gradJx(2,a,:) = (/0.0_pf, -Rz, Ry/)

      gradJy(1,a,:) = (/-Pz, 0.0_pf, Px/)
      gradJy(2,a,:) = (/Rz, 0.0_pf, -Rx/)

      gradJz(1,a,:) = (/Py, -Px, 0.0_pf/)
      gradJz(2,a,:) = (/-Ry, Rx, 0.0_pf/)

      normas2(1) = normas2(1) + Pz**2 + Rz**2 + Py**2 + Ry**2
      normas2(2) = normas2(2) + Pz**2 + Rz**2 + Px**2 + Rx**2
      normas2(3) = normas2(3) + Py**2 + Ry**2 + Px**2 + Rx**2

    end do 

    J = angular_geral(R, P)

    difJ = J - J0

    gradJx = gradJx * difJ(1)/normas2(1)
    gradJy = gradJy * difJ(2)/normas2(2)
    gradJz = gradJz * difJ(3)/normas2(3)

    do a = 1, size(m)

      R(a,:) = R(a,:) - gradJx(1,a,:) - gradJy(1,a,:) - gradJz(1,a,:)
      P(a,:) = P(a,:) - gradJx(2,a,:) - gradJy(2,a,:) - gradJz(2,a,:)

    end do

  end subroutine angular_correcao

end module angular