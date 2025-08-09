# ============================================================================
# Arquivo CMake para os integradores
# ----------------------------------------------------------------------------
# Este arquivo eh chamado pelo CMakeLists.txt na raiz do programa, e contem
# o conjunto dos integradores disponiveis.
# Quando um integrador for adicionado, basta alterar aqui e recompilar que ele
# ja ficara disponivel para uso.
# ============================================================================

set(INTEGRADORES_SOURCE
  # Metodos de Euler
  "./src/calcs/integracao/euler/euler_exp.f90" # (explicito)
  "./src/calcs/integracao/euler/euler_imp.f90" # (implicito)

  # Metodos de Runge-Kutta
  "./src/calcs/integracao/rk/rungekutta2.f90" 
  "./src/calcs/integracao/rk/rungekutta3.f90" 
  "./src/calcs/integracao/rk/rungekutta4.f90" 

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

CACHE STRING "Lista de integradores")

# Listas formatadas para o Fortran
set(INTEGRADORES_USE "")
set(INTEGRADORES_CLASSES "")
set(INTEGRADORES_LISTA "")
set(INTEGRADORES_SELECT "")

# Contador para facilitar a formatacao
set(CONTADOR 0)
list(LENGTH INTEGRADORES_SOURCE QUANTIDADE_INTEGRADORES)

foreach(src ${INTEGRADORES_SOURCE})
  # Adiciona no contador
  math(EXPR CONTADOR "${CONTADOR} + 1")

  get_filename_component(fname ${src} NAME_WE) # ex: euler_exp
  string(TOUPPER ${fname} FNAME_UP) # caso queira usar maiúsculas

  # 1. USE
  string(APPEND INTEGRADORES_USE "  USE ${fname}\n")

  # 2. TYPE
  string(APPEND INTEGRADORES_CLASSES "  TYPE(integracao_${fname}), TARGET, SAVE :: INT_${fname}\n")

  # 3. LISTA
  if (CONTADOR EQUAL 1)
    string(APPEND INTEGRADORES_LISTA "  WRITE(*,'(A)',ADVANCE='NO') '    ${fname}, '\n")
  elseif(CONTADOR EQUAL QUANTIDADE_INTEGRADORES)
    string(APPEND INTEGRADORES_LISTA "  WRITE(*,'(A)') '${fname}.'\n")
  else()
    string(APPEND INTEGRADORES_LISTA "  WRITE(*,'(A)',ADVANCE='NO') '${fname}, '\n")
  endif()

  # 4. DEFINIR_METODO
  string(APPEND INTEGRADORES_SELECT "    CASE ('${fname}'); integrador => INT_${fname}\n")
endforeach()

# Gera o arquivo a partir do template
configure_file(
  ${CMAKE_SOURCE_DIR}/src/calcs/integracao/integradores.f90.in
  ${CMAKE_BINARY_DIR}/src/calcs/integracao/integradores.f90
  @ONLY
)