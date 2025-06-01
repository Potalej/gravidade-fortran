# ----------------------------------------------------------------------------
# Arquivo CMake para o Gravidade-Fortran
#
# Na raiz do diretorio, utilize:
#   > cmake -B build
# Para compilar com o Ninja, use:
#   > cmake -B build -G Ninja
# ----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.7)
project("gravidade" VERSION 0.4.2
        DESCRIPTION "Simulacoes de gravidade em Fortran!"
        LANGUAGES Fortran)
enable_language(Fortran)

# Data e hora de compilacao
string(TIMESTAMP BUILD_DATE "%Y-%m-%d")
string(TIMESTAMP BUILD_TIME "%H:%M:%S")

# Flags para compilacao
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -march=native -fopenmp")
message("\n***************************************")
message("#       gravidade-fortran v${CMAKE_PROJECT_VERSION}        #")
message("***************************************")

# ----------------------------------------------------------------------------
# Precisao
# ----------------------------------------------------------------------------
set(PRECISAO "64" CACHE STRING "Escolha a precisão: 32, 64 ou 128")
set_property(CACHE PRECISAO PROPERTY STRINGS 32 64 128)

# Escolha da precisão
if(PRECISAO STREQUAL "32")
  set(PF_KIND "selected_real_kind(6, 37)")
elseif(PRECISAO STREQUAL "64")
  set(PF_KIND "selected_real_kind(15, 300)")
elseif(PRECISAO STREQUAL "128")
  set(PF_KIND "selected_real_kind(33, 4931)")
else()
  message(FATAL_ERROR "Precisão inválida: ${PRECISAO}. Use 32, 64 ou 128")
endif()

# Gera o módulo tipos.f90 automaticamente
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/src/tipos.f90.in
  ${CMAKE_CURRENT_BINARY_DIR}/tipos.f90
  @ONLY
)

add_library(tipos
  ${CMAKE_CURRENT_BINARY_DIR}/tipos.f90
)

# ----------------------------------------------------------------------------
# Versionamento
# ----------------------------------------------------------------------------
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/version.f90.in
  ${CMAKE_CURRENT_BINARY_DIR}/version.f90
  @ONLY
)

# ----------------------------------------------------------------------------
# Executaveis utilizados no programa
# ----------------------------------------------------------------------------
add_executable(gravidade
  "main.f90"

  "./src/calcs/simulacao/simulacao_sorteio.f90"
  "./src/calcs/simulacao/simulacao_vi.f90"
  "./src/calcs/simulacao/simulacao.f90"

  "./src/arquivo/arquivos_json.f90"
  "./src/arquivo/arquivos.f90"
  "./src/arquivo/plot.f90"

  ${CMAKE_CURRENT_BINARY_DIR}/version.f90
)
target_include_directories(gravidade PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

# ----------------------------------------------------------------------------
# Procura pelo OpenBLAS
# ----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(OpenBLAS REQUIRED)

IF (NOT OpenBLAS_LIB)
  message(FATAL_ERROR "OpenBLAS nao encontrado!!!")
ENDIF()

# ----------------------------------------------------------------------------
# Procura pelo JSONFortran
# ----------------------------------------------------------------------------
find_package(JsonFortran REQUIRED)

include_directories(${JSONFortran_INCLUDE_DIRS})
IF (NOT JSONFortran_LIB)
  message(FATAL_ERROR "JSON Fortran nao encontrado!!!")
ENDIF()

# ----------------------------------------------------------------------------
# Bibliotecas locais
# ----------------------------------------------------------------------------

# JSON-FILE-MOD
add_library(json_utils_mod 
  "./src/arquivo/json_utils_mod.f90"
)
target_link_libraries(json_utils_mod PRIVATE tipos)

# ALEATORIO
add_library(aleatorio 
  "./src/libs/aleatorio.f90" 
)
target_link_libraries(aleatorio PRIVATE tipos)

# ARVORE BINARIA
add_library(octree 
  "./src/libs/octree.f90"
)
target_link_libraries(octree PRIVATE tipos)

# Mecanica
add_library(mecanica 
  "./src/calcs/mecanica.f90" 
  "./src/calcs/auxiliares.f90" 
  "./src/calcs/mecanica/correcao.f90" 
  "./src/calcs/mecanica/colisao.f90"
)
target_link_libraries(mecanica PRIVATE
  ${OpenBLAS_LIB} # utilizada no corretor
  tipos
)

# Condicionamento de valores iniciais
add_library(condicionamento
  "./src/calcs/condicionamento/condicoesIniciais.f90"  
  "./src/calcs/condicionamento/condicionamento.f90"  
)
target_link_libraries(condicionamento PRIVATE
  mecanica
  ${OpenBLAS_LIB}
  tipos
)

# Metodos de integracao
add_library(integracao 
  "./src/calcs/integracao/forcas.f90" 
  "./src/calcs/integracao/forcas_mi.f90" 
  
  "./src/calcs/integracao/integrador.f90" 
  "./src/calcs/integracao/rk/rungekutta.f90" 
  "./src/calcs/integracao/rk/rungekutta2.f90" 
  "./src/calcs/integracao/rk/rungekutta3.f90" 
  "./src/calcs/integracao/rk/rungekutta4.f90" 
  "./src/calcs/integracao/simpleticos/verlet.f90" 
  "./src/calcs/integracao/simpleticos/euler_simp.f90" 

  # Metodos de Ruth
  "./src/calcs/integracao/simpleticos/ruth/ruth3.f90" 
  "./src/calcs/integracao/simpleticos/ruth/ruth4.f90" 

  # Metodos de Runge-Kutta-Nystrom
  "./src/calcs/integracao/simpleticos/rkn/rkn551.f90" 
  "./src/calcs/integracao/simpleticos/rkn/rkn671.f90" 

  # Metodos de Stormer-Verlet Compostos
  "./src/calcs/integracao/simpleticos/svc/svcp8s15.f90" 
  "./src/calcs/integracao/simpleticos/svc/svcp10s35.f90" 
  
  # Euler (explicito, tradicional)
  "./src/calcs/integracao/euler_exp.f90" 
  # Euler (implicito, tradicional)
  "./src/calcs/integracao/euler_imp.f90" 
)
target_link_libraries(integracao PRIVATE
  mecanica # usa o corretor e o de colisoes
  tipos
)

# Utilitarios de OS
# Alguns recursos depENDem do sistema operacional, entao precisam ser
# filtrados antes da inclusao no programa.
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  add_library(os "./src/os/linux.f90")

elseif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  add_library(os "./src/os/windows.f90")
ENDIF()

# ----------------------------------------------------------------------------
# Procura pelo Python
# A exibicao dos graficos eh feita com o matplotlib, entao o Python precisa
# estar instalado.
# ----------------------------------------------------------------------------
find_program(PYTHON "python" QUIET)
IF(PYTHON)
  find_package(Python COMPONENTS Interpreter REQUIRED QUIET)

  # procura pelo matplotlib
  execute_process(
    COMMAND pip show matplotlib
    RESULT_VARIABLE EXIT_CODE
    OUTPUT_QUIET
  )

  IF(NOT EXIT_CODE EQUAL 0)
    # instala o matplotlib
    execute_process(
      COMMAND ${Python_EXECUTABLE} -m pip install matplotlib
      RESULT_VARIABLE EXIT_MATPLOTLIB
      OUTPUT_QUIET
    )

    IF(NOT EXIT_MATPLOTLIB EQUAL 0)
      message(FATAL_ERROR "A depENDencia 'matplotlib' nao foi encontrada e houve um erro em sua instalacao. Nao sera possivel criar graficos.")
    ENDIF()
  ENDIF()

  add_library(python_lib 
    "./src/libs/python.f90"
  )
  message(STATUS "Python localizado.")
ELSE()
  message(FATAL_ERROR "Python nao encontrado! Nao sera possivel criar graficos.")
ENDIF()

# ----------------------------------------------------------------------------
# Linka as bibliotecas no executavel
# ----------------------------------------------------------------------------
target_link_libraries(gravidade
  tipos
  os
  integracao
  mecanica
  python_lib
  condicionamento
  aleatorio
  octree
  ${JSONFortran_LIBRARIES}
  json_utils_mod
)

# ----------------------------------------------------------------------------
# Apos a compilacao faz uma copia do executavel para fora da pasta build
# ----------------------------------------------------------------------------
add_custom_command(TARGET gravidade POST_BUILD 
  COMMAND "${CMAKE_COMMAND}" -E copy 
     "$<TARGET_FILE:gravidade>"
     "../$<CONFIGURATION>/$<TARGET_FILE_NAME:gravidade>" 
  COMMENT "Copiando para o lado de fora")