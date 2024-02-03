# ----------------------------------------------------------------------------
# Arquivo CMake para o Gravidade-Fortran
#
# Na raiz do diretorio, utilize:
#   > cmake -B build
# Para compilar com o Ninja, use:
#   > cmake -B build -G Ninja
# ----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.7)
project("gravidade" VERSION 0.1
        DESCRIPTION "Simulacoes de gravidade em Fortran!"
        LANGUAGES Fortran)
enable_language(Fortran)

# Flags para compilacao
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -march=native -ffast-math -fopenmp -cpp")

# ----------------------------------------------------------------------------
# Executaveis utilizados no programa
# ----------------------------------------------------------------------------
add_executable(gravidade
  "main.f90"

  "./calcs/simulacao/simulacao_sorteio.f90"
  "./calcs/simulacao/simulacao_vi.f90"
  "./calcs/simulacao/simulacao.f90"

  "./arquivo/leitura.f90"
  "./arquivo/arquivos.f90"
  "./arquivo/plot.f90"
)

# ----------------------------------------------------------------------------
# Procura pelo OpenBLAS
# ----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(OpenBLAS REQUIRED)

# ----------------------------------------------------------------------------
# Bibliotecas locais
# ----------------------------------------------------------------------------

# Mecanica
add_library(mecanica 
  "./calcs/mecanica.f90" 
  "./calcs/auxiliares/auxiliares.f90" 
  "./calcs/mecanica/correcao.f90" 
  "./calcs/mecanica/colisao.f90" 
)
target_link_libraries(mecanica
  ${OpenBLAS_LIB} # utilizada no corretor
)

# Condicionamento de valores iniciais
add_library(condicionamento
  "./calcs/condicionamento/condicoesArtigo.f90"  
)
target_link_libraries(condicionamento
  mecanica
)

# Metodos de integracao
add_library(integracao 
  "./calcs/integracao/integrador.f90" 
  "./calcs/integracao/verlet.f90" 
  "./calcs/integracao/rungekutta.f90" 
  "./calcs/integracao/rungekutta4.f90" 
  "./calcs/integracao/rkf45.f90" 
)
target_link_libraries(integracao
  mecanica # usa o corretor e o de colisoes
)

# Utilitarios de OS
# Alguns recursos dependem do sistema operacional, entao precisam ser
# filtrados antes da inclusao no programa.
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  add_library(os "./os/linux.f90")

elseif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  add_library(os "./os/windows.f90")
endif()

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
      message(FATAL_ERROR "A dependencia 'matplotlib' nao foi encontrada e houve um erro em sua instalacao. Nao sera possivel criar graficos.")
    ENDIF()
  ENDIF()

  add_library(python_lib 
    "./libs/python.f90"
  )
  message(STATUS "Python localizado.")
ELSE()
  message(FATAL_ERROR "Python nao encontrado! Nao sera possivel criar graficos.")
ENDIF()

# ----------------------------------------------------------------------------
# Linka as bibliotecas no executavel
# ----------------------------------------------------------------------------
target_link_libraries(gravidade
  os
  integracao
  mecanica
  python_lib
  condicionamento
)

# ----------------------------------------------------------------------------
# Apos a compilacao faz uma copia do executavel para fora da pasta build
# ----------------------------------------------------------------------------
add_custom_command(TARGET gravidade POST_BUILD 
  COMMAND "${CMAKE_COMMAND}" -E copy 
     "$<TARGET_FILE:gravidade>"
     "../$<CONFIGURATION>/$<TARGET_FILE_NAME:gravidade>" 
  COMMENT "Copiando para o lado de fora")