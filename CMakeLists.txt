# ----------------------------------------------------------------------------
# Arquivo CMake para o Gravidade-Fortran
#
# Na raiz do diretorio, utilize:
#   > cmake -B build
# Para compilar com o Ninja, use:
#   > cmake -B build -G Ninja
# ----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.7)
project("gravidade" VERSION 0.3.0
        DESCRIPTION "Simulacoes de gravidade em Fortran!"
        LANGUAGES Fortran)
enable_language(Fortran)

# Flags para compilacao
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -march=native -fopenmp")
message("\n***************************************")
message("#       gravidade-fortran v${CMAKE_PROJECT_VERSION}        #")
message("***************************************")

# ----------------------------------------------------------------------------
# Executaveis utilizados no programa
# ----------------------------------------------------------------------------
add_executable(gravidade
  "main.f90"

  "./src/calcs/simulacao/simulacao_sorteio.f90"
  "./src/calcs/simulacao/simulacao_vi.f90"
  "./src/calcs/simulacao/simulacao.f90"

  "./src/arquivo/leitura.f90"
  "./src/arquivo/arquivos.f90"
  "./src/arquivo/plot.f90"
)

# ----------------------------------------------------------------------------
# Procura pelo OpenBLAS
# ----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(OpenBLAS REQUIRED)

IF (NOT OpenBLAS_LIB)
  message(FATAL_ERROR "OpenBLAS nao encontrado!!!")
ENDIF()

# ----------------------------------------------------------------------------
# Bibliotecas locais
# ----------------------------------------------------------------------------

# ALEATORIO
add_library(aleatorio 
  "./src/libs/aleatorio.f90" 
)

# ARVORE BINARIA
add_library(octree 
  "./src/libs/octree.f90"
)

# Mecanica
add_library(mecanica 
  "./src/calcs/mecanica.f90" 
  "./src/calcs/auxiliares.f90" 
  "./src/calcs/mecanica/correcao.f90" 
  "./src/calcs/mecanica/colisao.f90"
)
target_link_libraries(mecanica
  ${OpenBLAS_LIB} # utilizada no corretor
)

# Condicionamento de valores iniciais
add_library(condicionamento
  "./src/calcs/condicionamento/condicoesIniciais.f90"  
)
target_link_libraries(condicionamento
  mecanica
  ${OpenBLAS_LIB}
)

# Metodos de integracao
add_library(integracao 
  "./src/calcs/integracao/forcas.f90" 
  
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
target_link_libraries(integracao
  mecanica # usa o corretor e o de colisoes
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
  os
  integracao
  mecanica
  python_lib
  condicionamento
  aleatorio
  octree
)

# ----------------------------------------------------------------------------
# Apos a compilacao faz uma copia do executavel para fora da pasta build
# ----------------------------------------------------------------------------
add_custom_command(TARGET gravidade POST_BUILD 
  COMMAND "${CMAKE_COMMAND}" -E copy 
     "$<TARGET_FILE:gravidade>"
     "../$<CONFIGURATION>/$<TARGET_FILE_NAME:gravidade>" 
  COMMENT "Copiando para o lado de fora")