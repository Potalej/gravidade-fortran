! Correcao
!
! Calculo da correcao numerica via matriz jacobiana e utilizando
! das propriedades simpleticas do espaco de fases do sistema
! hamiltoniano.
! 

module correcao
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  use hamiltoniano
  use auxiliares
  use angular

  implicit none
  external dgesv
  private
  public corrigir

contains

  ! aplica a correcao
  subroutine corrigir (G, massas, posicoes, momentos, grads, gradsT, vetorCorrecao)
    implicit none
    real(pf), intent(in) :: G, massas(:)
    real(pf), intent(inout) :: posicoes(:,:), momentos(:,:), grads(:,:),gradsT(:,:),vetorCorrecao(:)
    integer  :: N, a, INFO, b
    integer :: pivos(10)
    real(pf) :: JJt(10, 10), vetG(10)
    ! real(pf), allocatable :: grads(:,:), gradsT(:,:), vetorCorrecao(:)

    N = size(massas)
    
    ! calcula a mariz normal
    JJt = matriz_normal(G, massas, posicoes, momentos, N, grads)

    ! vetor G
    vetG = Gx(G, massas, posicoes, momentos, N)

    ! resolve o sistema
    ! r8mat_fs(10, JJt, vetG, vetorCorrecao)
    call dgesv(10, 1, JJt, 10, pivos, vetG, 10, INFO)

    if (INFO == 0) then
      ! se tiver solucao, aplica a correcao
      do a = 1, 10
        grads(a,:) = vetG(a)*grads(a,:)
      end do
      gradsT = transpose(grads)

      do a = 1, 6*N
        vetorCorrecao(a) = sum(gradsT(a,:))
      end do

      ! aplica a correcao
      do a = 1, N
        posicoes(a,1) = posicoes(a,1) + vetorCorrecao(6*a-5)
        posicoes(a,2) = posicoes(a,2) + vetorCorrecao(6*a-4)
        posicoes(a,3) = posicoes(a,3) + vetorCorrecao(6*a-3)
        momentos(a,1) = momentos(a,1) + vetorCorrecao(6*a-2)
        momentos(a,2) = momentos(a,2) + vetorCorrecao(6*a-1)
        momentos(a,3) = momentos(a,3) + vetorCorrecao(6*a-0)
      end do
    end if

  end subroutine corrigir

  ! calcula o lado direito (G(x))
  function Gx (G, massas, posicoes, momentos, N)
    implicit none
    integer  :: N
    real(pf) :: G, massas(N), posicoes(N,3), momentos(N,3), Gx(10), J(3), Rcm(3), P(3)

    ! energia total
    Gx(1) = energia_total(G, massas, posicoes, momentos)

    ! momento angular
    J = angular_geral(posicoes, momentos)
    Gx(2) = J(1)
    Gx(3) = J(2)
    Gx(4) = J(3)

    ! momento linear
    P = momentoLinear_total(momentos)
    Gx(5) = P(1)
    Gx(6) = P(2)
    Gx(7) = P(3)

    ! centro de massas
    Rcm = centro_massas(massas, posicoes)
    Gx(8) = Rcm(1)
    Gx(9) = Rcm(2)
    Gx(10) = Rcm(3)

    Gx = - Gx
  
  end function Gx

  ! matriz normal
  function matriz_normal (G, massas, posicoes, momentos, N, grads)

    implicit none
    integer  :: N, gi, gj
    real(pf) :: G, massas(N), posicoes(N,3), momentos(N,3)
    real(pf),intent(inout) :: grads(10, 6*N)
    real(pf) :: matriz_normal(10,10)

    grads(1,:)=gradiente_energia(G, massas, posicoes, momentos, N)
    grads(2,:)=gradiente_angularX(posicoes, momentos, N)
    grads(3,:)=gradiente_angularY(posicoes, momentos, N)
    grads(4,:)=gradiente_angularZ(posicoes, momentos, N)
    grads(5,:)=gradiente_linearX(N)
    grads(6,:)=gradiente_linearY(N)
    grads(7,:)=gradiente_linearZ(N)
    grads(8,:)=gradiente_centroMassasX(massas, N)
    grads(9,:)=gradiente_centroMassasY(massas, N)
    grads(10,:)=gradiente_centroMassasZ(massas, N)

    do gi = 1, 10
      do gj = 1, 10
        matriz_normal(gi, gj) = DOT_PRODUCT(grads(gi,:), grads(gj,:))
      end do
    end do

  end function matriz_normal

  ! gradiente da energia
  function gradiente_energia (G, massas, Rs, Ps, N)

    implicit none
    integer  :: N, a, b
    real(pf) :: distancia3, distancia(3)
    real(pf) :: G, massas(N), Rs(N,3), Ps(N,3), gradiente_energia(1:6*N)

    gradiente_energia(:) = 0.0_pf

    do a = 1, N
      gradiente_energia(6*a-2) = Ps(a,1)/massas(a)
      gradiente_energia(6*a-1) = Ps(a,2)/massas(a)
      gradiente_energia(6*a-0) = Ps(a,3)/massas(a)
      ! gradiente_energia(a,:) = [ 0.0_pf, 0.0_pf, 0.0_pf, Ps(a,1)/massas(a), Ps(a,2)/massas(a), Ps(a,3)/massas(a) ]
      
      do b = 1, N
        if (b /= a) then
          distancia = Rs(b,:)-Rs(a,:)
          distancia3 = norm2(distancia)**3
          distancia = (massas(b)/distancia3) * distancia

          gradiente_energia(6*a-5) = gradiente_energia(6*a-5)+distancia(1)
          gradiente_energia(6*a-4) = gradiente_energia(6*a-4)+distancia(2)
          gradiente_energia(6*a-3) = gradiente_energia(6*a-3)+distancia(3)

          ! gradiente_energia(a,1) = gradiente_energia(a,1)+massas(b)*distancia(1)
          ! gradiente_energia(a,2) = gradiente_energia(a,2)+massas(b)*distancia(2)
          ! gradiente_energia(a,3) = gradiente_energia(a,3)+massas(b)*distancia(3)
        end if
      end do

      gradiente_energia(6*a-5) = gradiente_energia(6*a-5) * (-G*massas(a))
      gradiente_energia(6*a-4) = gradiente_energia(6*a-4) * (-G*massas(a))
      gradiente_energia(6*a-3) = gradiente_energia(6*a-3) * (-G*massas(a))

      ! gradiente_energia(a,1) = gradiente_energia(a,1) * (-G*massas(a))
      ! gradiente_energia(a,2) = gradiente_energia(a,2) * (-G*massas(a))
      ! gradiente_energia(a,3) = gradiente_energia(a,3) * (-G*massas(a))
    end do
    
  end function gradiente_energia

  ! gradiente do momento angular (em X)
  function gradiente_angularX (posicoes, momentos, N)
    implicit none
    integer  :: N, a
    real(pf) :: gradiente_angularX(6*N), posicoes(N,3), momentos(N,3)

    do a = 1, N
      gradiente_angularX(6*a-5) = 0
      gradiente_angularX(6*a-4) = momentos(a,3)
      gradiente_angularX(6*a-3) = - momentos(a,2)
      gradiente_angularX(6*a-2) = 0
      gradiente_angularX(6*a-1) = - posicoes(a,3)
      gradiente_angularX(6*a-0) = posicoes(a,2)
    end do

  end function gradiente_angularX

  ! gradiente do momento angular (em Y)
  function gradiente_angularY (posicoes, momentos, N)
    implicit none
    integer  :: N, a
    real(pf) :: posicoes(N,3), momentos(N,3), gradiente_angularY(6*N)

    do a = 1, N
      gradiente_angularY(6*a-5) = - momentos(a,3)
      gradiente_angularY(6*a-4) = 0
      gradiente_angularY(6*a-3) = momentos(a,1)
      gradiente_angularY(6*a-2) = posicoes(a,3)
      gradiente_angularY(6*a-1) = 0
      gradiente_angularY(6*a-0) = - posicoes(a,1)
    end do

  end function gradiente_angularY

  ! gradiente do momento angular (em Z)
  function gradiente_angularZ (posicoes, momentos, N)
    implicit none
    integer  :: N, a
    real(pf) :: posicoes(N,3), momentos(N,3), gradiente_angularZ(6*N)

    do a = 1, N
      gradiente_angularZ(6*a-5) = momentos(a,2)
      gradiente_angularZ(6*a-4) = - momentos(a,1)
      gradiente_angularZ(6*a-3) = 0
      gradiente_angularZ(6*a-2) = - posicoes(a,2)
      gradiente_angularZ(6*a-1) = posicoes(a,1)
      gradiente_angularZ(6*a-0) = 0
    end do

  end function gradiente_angularZ

  ! gradiente do momento linear (em X)
  function gradiente_linearX (N)
    implicit none
    integer :: N, a
    real(pf) :: gradiente_linearX(6*N)

    do a = 1, N
      gradiente_linearX(6*a-5) = 0.0_pf
      gradiente_linearX(6*a-4) = 0.0_pf
      gradiente_linearX(6*a-3) = 0.0_pf
      gradiente_linearX(6*a-2) = 1.0_pf
      gradiente_linearX(6*a-1) = 0.0_pf
      gradiente_linearX(6*a-0) = 0.0_pf
    end do

  end function gradiente_linearX

  ! gradiente do momento linear (em Y)
  function gradiente_linearY (N)
    implicit none
    integer :: N, a
    real(pf) :: gradiente_linearY(6*N)

    do a = 1, N
      gradiente_linearY(6*a-5) = 0.0_pf
      gradiente_linearY(6*a-4) = 0.0_pf
      gradiente_linearY(6*a-3) = 0.0_pf
      gradiente_linearY(6*a-2) = 0.0_pf
      gradiente_linearY(6*a-1) = 1.0_pf
      gradiente_linearY(6*a-0) = 0.0_pf
    end do

  end function gradiente_linearY

  ! gradiente do momento linear (em Z)
  function gradiente_linearZ (N)
    implicit none
    integer :: N, a
    real(pf) :: gradiente_linearZ(6*N)

    do a = 1, N
      gradiente_linearZ(6*a-5) = 0.0_pf
      gradiente_linearZ(6*a-4) = 0.0_pf
      gradiente_linearZ(6*a-3) = 0.0_pf
      gradiente_linearZ(6*a-2) = 0.0_pf
      gradiente_linearZ(6*a-1) = 0.0_pf
      gradiente_linearZ(6*a-0) = 1.0_pf
    end do

  end function gradiente_linearZ

  ! gradiente do centro de massas (em X)
  function gradiente_centroMassasX (massas, N)
    implicit none
    integer  :: N, a
    real(pf) :: M
    real(pf) :: massas(N), gradiente_centroMassasX(6*N)

    M = sum(massas)

    do a = 1, N
      gradiente_centroMassasX(6*a-5) = massas(a)/M
      gradiente_centroMassasX(6*a-4) = 0.0_pf
      gradiente_centroMassasX(6*a-3) = 0.0_pf
      gradiente_centroMassasX(6*a-2) = 0.0_pf
      gradiente_centroMassasX(6*a-1) = 0.0_pf
      gradiente_centroMassasX(6*a-0) = 0.0_pf
    end do

  end function gradiente_centroMassasX

  ! gradiente do centro de massas (em Y)
  function gradiente_centroMassasY (massas, N)
    implicit none
    integer  :: N, a
    real(pf) :: M
    real(pf) :: massas(N), gradiente_centroMassasY(6*N)

    M = sum(massas)

    do a = 1, N
      gradiente_centroMassasY(6*a-5) = 0.0_pf
      gradiente_centroMassasY(6*a-4) = massas(a)/M
      gradiente_centroMassasY(6*a-3) = 0.0_pf
      gradiente_centroMassasY(6*a-2) = 0.0_pf
      gradiente_centroMassasY(6*a-1) = 0.0_pf
      gradiente_centroMassasY(6*a-0) = 0.0_pf
    end do

  end function gradiente_centroMassasY

  ! gradiente do centro de massas (em Z)
  function gradiente_centroMassasZ (massas, N)
    implicit none
    integer  :: N, a
    real(pf) :: M
    real(pf) :: massas(N), gradiente_centroMassasZ(6*N)

    M = sum(massas)

    do a = 1, N
      gradiente_centroMassasZ(6*a-5) = 0.0_pf
      gradiente_centroMassasZ(6*a-4) = 0.0_pf
      gradiente_centroMassasZ(6*a-3) = massas(a)/M
      gradiente_centroMassasZ(6*a-2) = 0.0_pf
      gradiente_centroMassasZ(6*a-1) = 0.0_pf
      gradiente_centroMassasZ(6*a-0) = 0.0_pf
    end do

  end function gradiente_centroMassasZ

end module correcao