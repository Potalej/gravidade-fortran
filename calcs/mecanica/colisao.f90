! Funcoes relativas a colisao de particulas
! 
! Funcoes
! 
! = subroutine verificar_e_colidir (m, R, P)
! Passa por cada par de particulas avaliando as distancias e
! a taxa de variacao desta, verificando se trata-se de caso
! de colisao ou nao.
! 
! = subourtine colidir (ma, Ra, Pa, mb, Rb, Pb)
! Faz o calculo da colisao perfeitamente elastica entre
! duas particulas A e B dadas.
! 

module colisao
  use, intrinsic :: iso_fortran_env, only: pf=>real64
  implicit none
  private
  public verificar_e_colidir, colidir

contains

  subroutine verificar_e_colidir (m, R, P)
    implicit none
    real(pf) :: m(:), R(:,:), P(:,:), max_aprox = 3
    integer :: a, b
    do a = 2, size(m)
      do b = 1, a-1
        if (norm2(R(b,:)-R(a,:)) <= max_aprox) then
          ! Agora verifica o sinal da derivada da distancia
          ! e estiver negativo, eh porque estao se aproximando
          if (DOT_PRODUCT(R(b,:) - R(a,:), P(b,:)-P(a,:)) < 0) then
            call colidir (m(a), R(a,:), P(a,:), m(b), R(b,:), P(b,:))
            WRITE (*,*) 'colidiu'
          endif
        endif 
      end do
    end do

  end subroutine verificar_e_colidir

  subroutine colidir (ma, Ra, Pa, mb, Rb, Pb)
    implicit none
    real(pf) :: ma, mb, Ra(3), Pa(3), Rb(3), Pb(3)
    real(pf) :: ua(3), ub(3), Normal(3), Normal_(3), u1, u2, ua_p(3), ub_p(3)

    ! separa as velocidades
    ua = Pa/ma
    ub = Pb/mb

    ! vetor normal e normal unitario
    Normal = Rb - Ra
    Normal_ = Normal/norm2(Normal)

    ! calcula a componente tangente
    u1 = DOT_PRODUCT(ua, Normal_)
    u2 = DOT_PRODUCT(ub, Normal_)
    
    ! calcula as componentes do plano
    ua_p = ua - u1*Normal
    ub_p = ub - u2*Normal

    ! obtem as novas velocidades
    Pa = ma * (ua_p + (u1*(ma-mb)+2*mb*u2)/(ma+mb) * Normal_)
    Pb = mb * (ub_p + (u2*(mb-ma)+2*ma*u1)/(ma+mb) * Normal_)

  end subroutine colidir

end module colisao