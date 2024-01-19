cmake_minimum_required(VERSION 3.7)
project("gravidade" VERSION 0.1
        DESCRIPTION "Simulacoes de gravidade em Fortran!"
        LANGUAGES Fortran)
enable_language(Fortran)

#dependece lib
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(OpenBLAS REQUIRED)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -march=native -ffast-math -fopenmp")

add_executable(
  "gravidade" 
  "main.f90" 

  "./calcs/simulacao/simulacao_sorteio.f90"
  "./calcs/simulacao/simulacao_vi.f90"
  "./calcs/simulacao/simulacao.f90" 

  "./arquivo/leitura.f90"
  "./arquivo/arquivos.f90" 

  "./strings/strings.f90" 

  "./calcs/condicionamento/condicoesArtigo.f90" 

  "./calcs/auxiliares/auxiliares.f90" 
  "./calcs/auxiliares/analise.f90"

  "./calcs/mecanica/hamiltoniano.f90" 
  "./calcs/mecanica/angular.f90" 
  "./calcs/mecanica/correcao.f90" 
  "./calcs/mecanica/colisao.f90" 

  "./calcs/integracao/integrador.f90" 
  "./calcs/integracao/verlet.f90" 
  "./calcs/integracao/rungekutta.f90" 
  "./calcs/integracao/rungekutta4.f90" 
  "./calcs/integracao/rkf45.f90" 
)

target_link_libraries(gravidade
  ${OpenBLAS_LIB}
)