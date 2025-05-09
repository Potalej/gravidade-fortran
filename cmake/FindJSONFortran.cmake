find_path(JSONFortran_INCLUDE_DIR json_module.mod
    PATHS
        # Ambiente
        ENV JSONFORTRAN_DIR
        # Linux
        /usr/include
        /usr/include/json-fortran
        /usr/local/include
        /usr/local/include/json-fortran
        # MSYS2 (UCRT64)
        /ucrt64/include/json-fortran
        C:/msys64/ucrt64/include/json-fortran
        # Local
        D:/programas/msys64/ucrt64/include
)

find_library(JSONFortran_LIB NAMES jsonfortran
    PATHS
        # Ambiente
        ENV JSONFORTRAN_DIR
        # Linux
        /lib/
        /lib/json-fortran
        /lib64/
        /usr/lib
        /usr/lib/json-fortran
        /usr/lib64
        /usr/local/lib
        /usr/local/lib64
        # MSYS2 (UCRT64)
        /ucrt64/lib
        C:/msys64/ucrt64/lib
        # Local
        D:/programas/msys64/ucrt64/lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(JsonFortran DEFAULT_MSG JSONFortran_INCLUDE_DIR JSONFortran_LIB)

if(JsonFortran_FOUND)
    set(JSONFortran_LIBRARIES ${JSONFortran_LIB})
    set(JSONFortran_INCLUDE_DIRS ${JSONFortran_INCLUDE_DIR})
    message("-- JSONFortran_INCLUDE_DIR: ${JSONFortran_INCLUDE_DIR}")
    message("-- JSONFortran_LIB: ${JSONFortran_LIBRARIES}")
endif()