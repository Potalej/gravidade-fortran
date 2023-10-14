cmake_minimum_required(VERSION 3.7)
project("gravidade" LANGUAGES "Fortran")

#dependece lib
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(OpenBLAS REQUIRED)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -march=native -ffast-math -fopenmp")

add_executable(
  "gravidade" 
  "main.f90" 
  "./simulacao/simulacao.f90" 
  "./auxiliares/hamiltoniano.f90" 
  "./integracao/rungekutta4.f90" 
  "arquivo/arquivos.f90" 
  "auxiliares/strings.f90" 
  "auxiliares/angular.f90" 
  "auxiliares/auxiliares.f90" 
  "condicionamento/condicoesArtigo.f90" 
  "auxiliares/analise.f90" 
  "integracao/rungekutta.f90" 
  "integracao/rkf45.f90" 
  "auxiliares/correcao.f90" 
  "auxiliares/colisao.f90" 
  "integracao/verlet.f90" 
  "integracao/integrador.f90" 
)
target_link_libraries(gravidade
  ${OpenBLAS_LIB}
)