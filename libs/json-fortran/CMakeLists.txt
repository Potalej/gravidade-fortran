# lib/json-fortran/CMakeLists.txt
# json-fortran 9.0.3 (fev/2025)

add_library(jsonfortran STATIC
  json_file_module.F90
  json_kinds.F90
  json_module.F90
  json_parameters.F90
  json_string_utilities.F90
  json_value_module.F90
)

# Agora sim, o alvo jsonfortran já existe
set_target_properties(jsonfortran PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod
)

target_include_directories(jsonfortran PUBLIC
  ${CMAKE_BINARY_DIR}/mod
  ${CMAKE_CURRENT_SOURCE_DIR}
)
