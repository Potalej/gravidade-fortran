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
! 
! = energia_gradiente (m, R, P)
!   Calcula o vetor gradiente da energia total do sistema e retorna dois 
!   vetores de derivadas parciais.
! 
! = energia_derivada_posicao (m, R)
!   Calcula as derivadas parciais da energia total em relação às posições.
! 
! = energia_derivada_momentos (m, P)
!   Calcula as derivadas pariciais da eneria total em relação aos momentos.
! 

module hamiltoniano

  implicit none
  private
  public energia_cinetica, energia_potencial, energia_total, energia_gradiente, energia_correcao

contains

  ! Energia cinética
  real function energia_cinetica (m, P) result(soma_ec)
    implicit none
    real :: m(:), P(:,:)
    integer :: a
    
    soma_ec = 0
    do a = 1, size(m)
      soma_ec = soma_ec + norm2(P(a,:))**2 / (2 * m(a))
    end do 
  end function

  ! Energia potencial
  real function energia_potencial (m, R) result(soma_ep)
    implicit none
    real :: m(:), R(:,:)
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
    real :: m(:), R(:,:), P(:,:)
    
    soma = 0
    soma = energia_cinetica(m, P) + energia_potencial(m, R)
  end function

  ! Gradiente da energia total
  function energia_gradiente (m, R, P)
    implicit none
    real, intent(in) :: m(:), R(:,:), P(:,:)

    ! onde será armazenado o gradiente (2, N, n)
    real, dimension( 2, size(R), size(R,1) ) :: energia_gradiente

    ! derivada em relação às posições
    energia_gradiente(1,:,:) = energia_derivada_posicao(m, R)

    ! derivada em relação aos momentos
    energia_gradiente(2,:,:) = energia_derivada_momentos(m, P)

  end function energia_gradiente

  ! Derivadas parciais da energia total em relação às posições
  function energia_derivada_posicao (m, R)
    implicit none
    real, intent(in)           :: m(:), R(:,:)

    real, dimension( size(R), size(R,1) )   :: energia_derivada_posicao ! vetor de derivadas
    integer                                 :: a, b ! iteradores
    real                                    :: rab3 ! para salvar o cubo da distancia entre dois corpos

    ! define como zero para começar
    energia_derivada_posicao(:,:) = 0

    do a = 2, size(m)
      do b = 1, a-1
        ! cubo da distancia entre os corpos
        rab3 = norm2( R(b,:) - R(a,:) )**3

        ! salva a primeira componente (a)
        energia_derivada_posicao(a,:) = energia_derivada_posicao(a,:) - m(b)*m(a)*( R(b,:) - R(a,:) )/rab3

        ! salva a segunda componente (b)
        energia_derivada_posicao(b,:) = energia_derivada_posicao(b,:) + energia_derivada_posicao(a,:)
      end do 
    end do 
  end function energia_derivada_posicao

  ! Derivadas parcais da energia total em relação aos momentos
  function energia_derivada_momentos (m, P)
    implicit none
    real, intent(in)         :: m(:), P(:,:)

    real, dimension( size(P), size(P,1) ) :: energia_derivada_momentos ! vetor de derivadas
    integer                               :: a ! iterador

    do a = 1, size(m)
      energia_derivada_momentos(a,:) = P(a,:) / m(a)
    end do

  end function energia_derivada_momentos

  ! Correção utilizando o gradiente de energia
  subroutine energia_correcao (m, R, P, E0)
    implicit none
    real, intent(in)                       :: m(:)
    real, intent(inout)                    :: R(:,:), P(:,:)
    real, intent(in)                       :: E0 ! energia total inicial

    real, dimension(2, size(R), size(R,1)) :: energia_grad
    real                                   :: norma_grad2, fator

    ! calcula a energia no momento atual
    E = energia_total(m, R, P)

    ! calcula o gradiente da energia
    energia_grad = energia_gradiente(m, R, P)

    ! norma do gradiente
    norma_grad2 = norm2(energia_grad)**2

    ! fator
    fator = (E-E0) / norma_grad2
    energia_grad = fator * energia_grad

    ! aplica a correção
    R = R - energia_grad(1,:,:)
    P = P - energia_grad(2,:,:)

  end subroutine energia_correcao

end module hamiltoniano