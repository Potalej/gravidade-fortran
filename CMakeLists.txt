# ============================================================================
# Arquivo CMake para o Gravidade-Fortran
#
# Na raiz do diretorio, utilize:
#   > cmake -B build
# Para compilar com o Ninja, use:
#   > cmake -B build -G Ninja
# ============================================================================

cmake_minimum_required(VERSION 3.15)
project("gravidade" VERSION 1.0.3
        DESCRIPTION "Simulacoes de gravidade em Fortran!"
        LANGUAGES Fortran)
enable_language(Fortran)

# Data e hora de compilacao
string(TIMESTAMP BUILD_DATE "%Y-%m-%d")
string(TIMESTAMP BUILD_TIME "%H:%M:%S")

# ============================================================================
# Compilacao
# ----------------------------------------------------------------------------
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -march=native -fopenmp")
# set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fcheck=all -Wall -fbacktrace -O3 -march=native -fopenmp") # Debug
message("\n******************************************")
message("#        gravidade-fortran v${CMAKE_PROJECT_VERSION}        #")
message("******************************************")
# ============================================================================

# ============================================================================
# GProf
# ----------------------------------------------------------------------------
option(GPROF "Habilitar criacao de perfil gprof" OFF)

if (GPROF)
  message(STATUS "[debug] GNU Profiler (gprof) ativado")
  add_compile_options(-pg)
  add_link_options(-pg)
endif()
# ============================================================================

# ============================================================================
# GPU
# ----------------------------------------------------------------------------
option(USAR_GPU "Utilizar GPU (OpenMP offload)" OFF)

if (USAR_GPU)
  message(STATUS "Compilando para comportar GPU.")
  set(USAR_GPU_BOOL ".TRUE.")
  add_compile_options(-foffload=nvptx-none)
  add_compile_definitions(USAR_GPU)
else()
  set(USAR_GPU_BOOL ".FALSE.")
endif()
# ============================================================================

# ============================================================================
# Precisao
# ----------------------------------------------------------------------------
set(PRECISAO "64" CACHE STRING "Escolha a precisão: 32, 64 ou 128")
set_property(CACHE PRECISAO PROPERTY STRINGS 32 64 128)

# Define flags de pre-processador conforme a precisão
if(PRECISAO STREQUAL "32")
  add_compile_definitions(REAL32)
elseif(PRECISAO STREQUAL "64")
  add_compile_definitions(REAL64)
elseif(PRECISAO STREQUAL "128")
  add_compile_definitions(REAL128)
else()
  message(FATAL_ERROR "Precisão inválida: ${PRECISAO}. Use 32, 64 ou 128")
endif()
# ============================================================================

# ============================================================================
# Versionamento
# ----------------------------------------------------------------------------
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/version.f90.in
  ${CMAKE_CURRENT_BINARY_DIR}/version.f90
  @ONLY
)
add_library(versao
  ${CMAKE_CURRENT_BINARY_DIR}/version.f90
)
# ============================================================================

# ============================================================================
# Executaveis utilizados no programa
# ----------------------------------------------------------------------------
add_executable(gravidade
  "main.f90"
)
target_link_libraries(versao)
target_include_directories(gravidade PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
# ============================================================================

# ============================================================================
# Procura pelo OpenBLAS
# ----------------------------------------------------------------------------
message("\n> Biblioteca externa: OpenBLAS")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(OpenBLAS REQUIRED)

IF (NOT OpenBLAS_LIB)
  message(FATAL_ERROR "OpenBLAS nao encontrado!!!")
ENDIF()
# ============================================================================

# ============================================================================
# Inclui o JSONFortran
# ----------------------------------------------------------------------------
message("\n> Biblioteca externa: JSON-Fortran")
option(FETCHJSONFORTRAN "Compilar o json-fortran localmente" OFF)

if (FETCHJSONFORTRAN OR PRECISAO EQUAL "32" OR PRECISAO EQUAL "128" OR NOT DEFINED ENV{CONDA_PREFIX})
  message(STATUS "Baixando do repositório para compilar localmente. Isso pode demorar um pouco...")
  
  include(FetchContent)

  FetchContent_Declare(
    jsonfortran
    GIT_REPOSITORY https://github.com/jacobwilliams/json-fortran.git
    GIT_TAG 9.2.0
    GIT_SHALLOW TRUE
  )

  set(JSON_REAL_KIND "REAL${PRECISAO}" CACHE INTERNAL "Tipo da variavel real do JSON-Fortran")
  
  FetchContent_MakeAvailable(jsonfortran)

  target_include_directories(jsonfortran
    INTERFACE
      ${jsonfortran_BINARY_DIR}
  )

elseif (DEFINED ENV{CONDA_PREFIX})
  message(STATUS "Ambiente conda ativo, verificando instalacao...")
  set(JSONFORTRAN_LIB "$ENV{CONDA_PREFIX}/lib/libjsonfortran.so")
  set(JSONFORTRAN_INC "$ENV{CONDA_PREFIX}/include")
  
  if (EXISTS "${JSONFORTRAN_LIB}")
    add_library(jsonfortran::jsonfortran SHARED IMPORTED)
    set_target_properties(jsonfortran::jsonfortran PROPERTIES
      IMPORTED_LOCATION "${JSONFORTRAN_LIB}"
      INTERFACE_INCLUDE_DIRECTORIES "${JSONFORTRAN_INC}"
    )
    add_library(jsonfortran ALIAS jsonfortran::jsonfortran)
    message(STATUS "JSON-Fortran habilitado via conda!")
  else()
    message(FATAL_ERROR "JSON-Fortran nao encontrado no conda! Instale no conda ou compile localmente atraves da flag -DFETCHJSONFORTRAN=ON")
  endif()
endif()
# ============================================================================

# ============================================================================
# Utilidades
# ----------------------------------------------------------------------------
# Esta contem rotinas e funcoes dos antigos modulos MECANICA e AUXILIARES,
# ou seja, contem tanto rotinas para resolucao de sistemas lineares, calculo
# de autovalores e afins, como tambem rotinas para calculo de energia total,
# momento angular, tensor de inercia e afins.
# ----------------------------------------------------------------------------
message("\n> Biblioteca externa: utilidades")

# diretorio
set(LOCAL_LIB_DIR_UT "${CMAKE_CURRENT_SOURCE_DIR}/libs/ncorpos_utilidades")
add_subdirectory(${LOCAL_LIB_DIR_UT} utilidades)
message(STATUS "Versão da biblioteca utilidades: ${NCORPOS_UTILIDADES_PROJECT_VERSION}")
# ============================================================================

# ============================================================================
#> Valores iniciais
# ----------------------------------------------------------------------------
# Biblioteca para geracao de valores iniciais. Foi separada do programa
# principal, mas serve para gerar valores iniciais e condicioná-los se for
# o caso.
# ----------------------------------------------------------------------------
message("\n> Biblioteca externa: valores-iniciais")

# diretorio
set(LOCAL_LIB_DIR_VI "${CMAKE_CURRENT_SOURCE_DIR}/libs/ncorpos_valores_iniciais")
add_subdirectory(${LOCAL_LIB_DIR_VI} valores_iniciais)
message(STATUS "Versão da biblioteca valores-iniciais: ${NCORPOS_VI_PROJECT_VERSION}")
# ============================================================================

message("")

# ============================================================================
# =-=-=-=-=-=-=-=-=-=-=-=-=-=- Bibliotecas locais =-=-=-=-=-=-=-=-=-=-=-=-=-=-
# ============================================================================

# ============================================================================
#> Tipos
# ----------------------------------------------------------------------------
# Esta biblioteca contem os tipos de variaveis reais utilizadas globalmente,
# incluindo o "pf". Eh chamada por todas as bibliotecas e modulos.
# ----------------------------------------------------------------------------
add_library(tipos
"./src/tipos.F90"
)
# ============================================================================

# ============================================================================
#> json_utils_mod
# ----------------------------------------------------------------------------
# Esta biblioteca contem funcionalidades basicas do Json-Fortran e rotinas
# para facilitar o uso local de variaveis de tipo JSON_VALUE.
# ----------------------------------------------------------------------------
add_library(json_utils_mod 
  "./src/libs/json_utils_mod.f90"
)
target_link_libraries(json_utils_mod PRIVATE tipos)
target_link_libraries(json_utils_mod PRIVATE jsonfortran)
# ============================================================================

# ============================================================================
#> Conexao
# ----------------------------------------------------------------------------
# Esta biblioteca eh utilizada para fazer a transmissao de dados via socket
# quando se deseja visualizar simulacoes em tempo real. Eh chamada pelo
# modulo SIMULACAO principal.
#
# Como ela depende do sistema operacional utilizado, eh necessario fazer
# a verificacao devida.
# ----------------------------------------------------------------------------
IF (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  add_library(conexao_mod STATIC
    "./src/socket/linux/conexao.f90"
    "./src/socket/linux/server_mod.f90"
  )
# ----------------------------------------------------------------------------
# No caso do Windows, eh preciso adicionar as rotinas do WSA para que o
# socket funcione corretamente
# ----------------------------------------------------------------------------
ELSEIF (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  add_library(conexao_mod STATIC
    "./src/socket/windows/conexao.f90"
    "./src/socket/windows/server_mod.f90"
    "./src/socket/windows/wsa_mod.f90"
  )
  target_link_libraries(conexao_mod PRIVATE ws2_32)
ENDIF()
target_link_libraries(conexao_mod PRIVATE tipos)
# ============================================================================

# ============================================================================
#> Octree
# ----------------------------------------------------------------------------
# Esta biblioteca implementa a separacao de vetores 3D em uma octree (arvore
# "octanaria"). No momento, a arvore eh utilizada somente pelo modulo de
# colisoes, quando desejado.
# ----------------------------------------------------------------------------
add_library(octree 
  "./src/libs/octree.f90"
)
target_link_libraries(octree PRIVATE tipos)
# ============================================================================

# ============================================================================
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= Modulos =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# ============================================================================

# ============================================================================
#> Arquivos
# ----------------------------------------------------------------------------
# Este modulo contem as funcionalidades de arquivos .txt, .csv e .json que
# sao utilizados pelo programa.
# ----------------------------------------------------------------------------
add_library(arquivos
  # Principal (reune os dois)
  "./src/arquivo/arquivos_mod.f90"

  # Arquivos gerais
  "./src/arquivo/arquivos.f90"
  "./src/arquivo/diretorio.f90"
  "./src/arquivo/string_utils.f90"

  # Arquivos JSON
  "./src/arquivo/arquivos_json.f90"
)
target_link_libraries(arquivos PRIVATE 
  tipos
  json_utils_mod 
)
# ============================================================================

# ============================================================================
#> Corretor
# ----------------------------------------------------------------------------
# Contem rotinas para, dado um estado no espaco de fases e as integrais que
# deveriam estar sendo satisfeitas, aplicar uma projecao do estado para obter
# um estado que atenda - aproximadamente - tais integrais. Eh utilizado apos
# integracoes numericas, permitindo melhorar os resultados.
# Depende da biblioteca "utilidades"
# ----------------------------------------------------------------------------
add_library(corretor
  "./src/calcs/correcao.f90" 
)
target_link_libraries(corretor PRIVATE
  ${OpenBLAS_LIB} # utilizada no corretor
  utilidades
  tipos
)

# ============================================================================
#> Colisoes
# ----------------------------------------------------------------------------
# Contem rotinas para aplicar choques perfeitamente elasticos entre dois
# corpos. Contem um verificador direto de colisoes, mas tambem importa a
# biblioteca "OCTREE" e permite utilizar a arvore para verificar colisoes.
# ----------------------------------------------------------------------------
add_library(colisoes
  "./src/calcs/colisao.f90" 
)
target_link_libraries(colisoes PRIVATE
  tipos
  octree
)
# ============================================================================

# ============================================================================

# ============================================================================
#> Integracao
# ----------------------------------------------------------------------------
# Modulo que centraliza a integracao. Todos os metodos de integracao numerica
# sao adicionados em
#
# "src/calcs/integracao/CMakeLists.txt"
#
# e sao chamados aqui, com as funcoes de forca e dos modulos centralizadores.
# ----------------------------------------------------------------------------
add_subdirectory(src/calcs/integracao) # Metodos de integracao
add_library(integracao 
  "./src/calcs/forcas/forcas.F90" 
  "./src/calcs/forcas/forcas_md.F90" 
  "./src/calcs/forcas/forcas_mi.F90"
  
  "./src/calcs/integracao/integrador.f90" 
  ${CMAKE_CURRENT_BINARY_DIR}/src/calcs/integracao/integradores.f90

  ${INTEGRADORES_SOURCE}
)
target_link_libraries(integracao PRIVATE
  tipos
  json_utils_mod
)
# ============================================================================

# ============================================================================
# Procura pelo Python
# ----------------------------------------------------------------------------
# A exibicao dos graficos eh feita com o matplotlib, entao o Python precisa
# estar instalado. Se nao estiver, um aviso sera enviado mas a biblioteca sera
# adicionada de qualquer jeito, entao caso o Python e o matplotlib sejam
# instalados posteriormente, o programa funcionara.
# ----------------------------------------------------------------------------
add_library(pyplot "./src/libs/pyplot.f90")
target_link_libraries(pyplot PRIVATE 
  tipos
  arquivos
)
# para detectar dentro de ambientes conda
set(Python_FIND_STRATEGY LOCATION)
set(Python_FIND_VIRTUALENV FIRST)
find_program(PYTHON NAMES python python3 QUIET)

message("> Python")
IF(PYTHON)
  message(STATUS "Python localizado (${PYTHON})")
  find_package(Python COMPONENTS Interpreter REQUIRED QUIET)

  # procura pelo matplotlib
  execute_process(
    COMMAND ${Python_EXECUTABLE} -m pip show matplotlib
    RESULT_VARIABLE EXIT_CODE
    OUTPUT_QUIET
    ERROR_QUIET
  )

  # Se nao achar, so da o WARNING
  IF(NOT EXIT_CODE EQUAL 0 AND NOT EXIT_MATPLOTLIB EQUAL 0)
    message(WARNING "A dependencia 'matplotlib' nao foi encontrada. Nao sera possivel criar graficos.")
  ELSE()
    message(STATUS "matplotlib localizado.\n")
  ENDIF()
ELSE()
  message(WARNING "Python nao encontrado! Nao sera possivel criar graficos.")
ENDIF()
# ============================================================================

# ============================================================================
# Modulos para simulacao
# ----------------------------------------------------------------------------
# O modulo de simulacao centraliza todo o sistema e eh efetivamente o que faz
# o programa funcionar, depois do MAIN.
# 
# O modulo simulacao_sorteio contem as rotinas para sortear e salvar os dados,
# alem de simular se for desejado.
#
# O modulo simulacao_vi contem as rotinas para ler e rodar simulacoes com
# valores predefinidos em algum arquivo de valores iniciais.
# 
# Os outros dois modulos sao chamados por esses primeiros, unindo os processos
# de integracao numerica, colisoes, correcao, arquivos etc.
# ----------------------------------------------------------------------------
add_library(simulacao 
  "./src/calcs/simulacao/simulacao_sorteio.f90"
  "./src/calcs/simulacao/simulacao_vi.f90"
  "./src/calcs/simulacao/simulacao.f90"
  "./src/calcs/simulacao/simulador.f90"
)
target_link_libraries(simulacao PRIVATE
  tipos
  versao
  arquivos
  integracao
  conexao_mod
  utilidades
  valores_iniciais
  colisoes
  corretor
)
# ============================================================================

# ============================================================================
# Linka as bibliotecas no executavel
# ----------------------------------------------------------------------------
target_link_libraries(gravidade
  tipos
  pyplot
  simulacao
  integracao
)
# ============================================================================

# ============================================================================
# Apos a compilacao faz uma copia do executavel para fora da pasta build
# ----------------------------------------------------------------------------
add_custom_command(TARGET gravidade POST_BUILD 
  COMMAND "${CMAKE_COMMAND}" -E copy 
     "$<TARGET_FILE:gravidade>"
     "../$<CONFIGURATION>/$<TARGET_FILE_NAME:gravidade>" 
  COMMENT "Copiando para o lado de fora")
# ----------------------------------------------------------------------------
