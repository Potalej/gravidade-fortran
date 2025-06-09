MODULE tipos
  IMPLICIT NONE
  INTEGER, PARAMETER, PUBLIC :: pf32  = SELECTED_REAL_KIND(6, 37)    ! precisao simples (~32 bits)
  INTEGER, PARAMETER, PUBLIC :: pf64  = SELECTED_REAL_KIND(15, 307)  ! precisao dupla (~64 bits)
  INTEGER, PARAMETER, PUBLIC :: pf128 = SELECTED_REAL_KIND(33, 4931) ! precisao quadrupla (~128 bits)

#ifdef REAL32
  INTEGER, PARAMETER, PUBLIC :: pf = pf32   !! Tipo real padrao [4 bytes]
#elif REAL64
  INTEGER, PARAMETER, PUBLIC :: pf = pf64   !! Tipo real padrao [8 bytes]
#elif REAL128
  INTEGER, PARAMETER, PUBLIC :: pf = pf128  !! Tipo real padrao [16 bytes]
#else
  INTEGER, PARAMETER, PUBLIC :: pf = pf64   !! Tipo real padrao se nao especificado [8 bytes]
#endif

END MODULE tipos