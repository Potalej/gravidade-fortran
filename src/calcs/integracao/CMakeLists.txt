# ============================================================================
# Arquivo CMake para os integradores
# ----------------------------------------------------------------------------
# Este arquivo eh chamado pelo CMakeLists.txt na raiz do programa, e contem
# o conjunto dos integradores disponiveis.
# Quando um integrador for adicionado, basta alterar aqui e recompilar que ele
# ja ficara disponivel para uso.
# ============================================================================

set(INTEGRADORES_TRADICIONAIS
  # Metodos de Euler
  "./src/calcs/integracao/euler/euler_exp.f90" # (explicito)
  "./src/calcs/integracao/euler/euler_imp.f90" # (implicito)

  # Metodos de Runge-Kutta
  "./src/calcs/integracao/rk/rungekutta2.f90" 
  "./src/calcs/integracao/rk/rungekutta3.f90" 
  "./src/calcs/integracao/rk/rungekutta4.f90" 
CACHE STRING "Lista de integradores tradicionais")

set(INTEGRADORES_SIMPLETICOS
  # Metodos Simpleticos basicos
  "./src/calcs/integracao/simpleticos/verlet.f90"     # (Verlet)
  "./src/calcs/integracao/simpleticos/euler_simp.f90" # (Euler Simpletico)

  # Metodos de Ruth
  "./src/calcs/integracao/simpleticos/ruth/ruth3.f90" # (ordem 3)
  "./src/calcs/integracao/simpleticos/ruth/ruth4.f90" # (ordem 4)

  # Metodos de Runge-Kutta-Nystrom
  "./src/calcs/integracao/simpleticos/rkn/rkn551.f90" # (ordem 5)
  "./src/calcs/integracao/simpleticos/rkn/rkn671.f90" # (ordem 6)

  # Metodos de Stormer-Verlet Compostos
  "./src/calcs/integracao/simpleticos/svc/svcp8s15.f90"  # (ordem 8)
  "./src/calcs/integracao/simpleticos/svc/svcp10s35.f90" # (ordem 10)
CACHE STRING "Lista de integradores tradicionais")

set(INTEGRADORES_MULTIPASSO
  # Adam-Bashforth
  "./src/calcs/integracao/multipasso/ab/ab2.f90"
  "./src/calcs/integracao/multipasso/ab/ab3.f90"
  "./src/calcs/integracao/multipasso/ab/ab4.f90"
  "./src/calcs/integracao/multipasso/ab/ab5.f90"
CACHE STRING "Lista de integradores multipasso")

set(INTEGRADORES_SOURCE
  ${INTEGRADORES_TRADICIONAIS}
  ${INTEGRADORES_SIMPLETICOS}
  ${INTEGRADORES_MULTIPASSO}
CACHE STRING "Lista de integradores")

# Listas formatadas para o Fortran
set(INTEGRADORES_USE "")
set(INTEGRADORES_CLASSES "")
set(INTEGRADORES_LISTA "")
set(INTEGRADORES_SELECT "")

foreach(src ${INTEGRADORES_SOURCE})
  # Captura o arquivo
  get_filename_component(fname ${src} NAME_WE) # ex: euler_exp
  string(TOUPPER ${fname} FNAME_UP) # caso queira usar maiúsculas

  # USE no script
  string(APPEND INTEGRADORES_USE "  USE ${fname}\n")

  # TYPE no script
  string(APPEND INTEGRADORES_CLASSES "  TYPE(integracao_${fname}), TARGET, SAVE :: INT_${fname}\n")

  # CASE do metodo no script
  string(APPEND INTEGRADORES_SELECT "    CASE ('${fname}'); integrador => INT_${fname}\n")
endforeach()

# Agora gera as strings
# Contador para facilitar a formatacao
set(CONTADOR 0)

# Primeiro, a string de tradicionais
list(LENGTH INTEGRADORES_TRADICIONAIS QUANTIDADE_INTEGRADORES)
string(APPEND INTEGRADORES_LISTA "WRITE(*,*)\n")
string(APPEND INTEGRADORES_LISTA "\n  WRITE(*,'(A)') '     - TRADICIONAIS:'\n")
foreach(src ${INTEGRADORES_TRADICIONAIS})
  # Captura o arquivo
  get_filename_component(fname ${src} NAME_WE)
  string(TOUPPER ${fname} FNAME_UP)

  # Adiciona no contador
  math(EXPR CONTADOR "${CONTADOR} + 1")

  # Adiciona na string
  if (CONTADOR EQUAL 1)
    string(APPEND INTEGRADORES_LISTA "  WRITE(*,'(A)',ADVANCE='NO') '         ${fname}, '\n")
  elseif(CONTADOR EQUAL QUANTIDADE_INTEGRADORES)
    string(APPEND INTEGRADORES_LISTA "  WRITE(*,'(A)') '${fname}.'\n")
  else()
    string(APPEND INTEGRADORES_LISTA "  WRITE(*,'(A)',ADVANCE='NO') '${fname}, '\n")
  endif()
endforeach()

string(APPEND INTEGRADORES_LISTA "WRITE(*,*)\n")

# Agora os simpleticos
set(CONTADOR 0)
list(LENGTH INTEGRADORES_SIMPLETICOS QUANTIDADE_INTEGRADORES)
string(APPEND INTEGRADORES_LISTA "\n  WRITE(*,'(A)') '     - SIMPLETICOS:'\n")
foreach(src ${INTEGRADORES_SIMPLETICOS})
  # Captura o arquivo
  get_filename_component(fname ${src} NAME_WE)
  string(TOUPPER ${fname} FNAME_UP)

  # Adiciona no contador
  math(EXPR CONTADOR "${CONTADOR} + 1")

  # Adiciona na string
  if (CONTADOR EQUAL 1)
    string(APPEND INTEGRADORES_LISTA "  WRITE(*,'(A)',ADVANCE='NO') '         ${fname}, '\n")
  elseif(CONTADOR EQUAL QUANTIDADE_INTEGRADORES)
    string(APPEND INTEGRADORES_LISTA "  WRITE(*,'(A)') '${fname}.'\n")
  else()
    string(APPEND INTEGRADORES_LISTA "  WRITE(*,'(A)',ADVANCE='NO') '${fname}, '\n")
  endif()
endforeach()

string(APPEND INTEGRADORES_LISTA "WRITE(*,*)\n")

# Agora os multipasso
set(CONTADOR 0)
list(LENGTH INTEGRADORES_MULTIPASSO QUANTIDADE_INTEGRADORES)
string(APPEND INTEGRADORES_LISTA "\n  WRITE(*,'(A)') '     - Multipasso:'\n")
foreach(src ${INTEGRADORES_MULTIPASSO})
  # Captura o arquivo
  get_filename_component(fname ${src} NAME_WE)
  string(TOUPPER ${fname} FNAME_UP)

  # Adiciona no contador
  math(EXPR CONTADOR "${CONTADOR} + 1")

  # Adiciona na string
  if (CONTADOR EQUAL 1)
    string(APPEND INTEGRADORES_LISTA "  WRITE(*,'(A)',ADVANCE='NO') '         ${fname}, '\n")
  elseif(CONTADOR EQUAL QUANTIDADE_INTEGRADORES)
    string(APPEND INTEGRADORES_LISTA "  WRITE(*,'(A)') '${fname}.'\n")
  else()
    string(APPEND INTEGRADORES_LISTA "  WRITE(*,'(A)',ADVANCE='NO') '${fname}, '\n")
  endif()
endforeach()

# Gera o arquivo a partir do template
configure_file(
  ${CMAKE_SOURCE_DIR}/src/calcs/integracao/integradores.f90.in
  ${CMAKE_BINARY_DIR}/src/calcs/integracao/integradores.f90
  @ONLY
)