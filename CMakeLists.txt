cmake_minimum_required(VERSION 3.7)
project("gravidade" LANGUAGES "Fortran")
add_executable("gravidade" "main.f90" "./simulacao/simulacao.f90" "./auxiliares/hamiltoniano.f90")