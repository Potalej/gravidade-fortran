cmake_minimum_required(VERSION 3.7)
project("gravidade" LANGUAGES "Fortran")
add_executable("gravidade" "main.f90" "./simulacao/simulacao.f90" "./auxiliares/hamiltoniano.f90" "./integracao/rungekutta4.f90" "arquivo/arquivos.f90" "auxiliares/strings.f90" "auxiliares/angular.f90" "auxiliares/auxiliares.f90" "condicionamento/condicoesArtigo.f90" "auxiliares/analise.f90" )