SET(Open_BLAS_INCLUDE_SEARCH_PATHS
    # Linux
    /usr/include
    /usr/include/openblas
    /usr/include/openblas-base
    /usr/local/include
    /usr/local/include/openblas
    /opt/OpenBLAS/include
    # MinGW (padrao)
    C:/MinGW/include
    C:/MinGW/include/openblas
    # MSYS2 (padrao)
    C:/msys64/usr/include
    C:/msys64/usr/include/openblas
    # MSYS2 + MINGW64 (padrao)
    C:/msys64/mingw64/usr/include
    C:/msys64/mingw64/usr/include/openblas

    D:/programas/msys64/mingw64/include/openblas
    
    # Ambiente (PATH)
    $ENV{OpenBLAS_HOME}
    $ENV{OpenBLAS_HOME}/include
)

SET(Open_BLAS_LIB_SEARCH_PATHS
    # Linux
    /lib/
    /lib/openblas-base
    /lib64/
    /usr/lib
    /usr/lib/openblas-base
    /usr/lib64
    /usr/local/lib
    /usr/local/lib64
    /opt/OpenBLAS/lib
    # MinGW (padrao)
    C:/MinGW/lib
    C:/MinGW/lib/openblas
    # MSYS2 (padrao)
    C:/msys64/usr/lib
    C:/msys64/usr/lib/openblas
    # MSYS2 + MINGW64 (padrao)
    C:/msys64/mingw64/usr/lib
    C:/msys64/mingw64/usr/lib/openblas

    D:/programas/msys64/mingw64/lib

    # Ambiente (PATH)
    $ENV{OpenBLAS}cd
    $ENV{OpenBLAS}/lib
    $ENV{OpenBLAS_HOME}
    $ENV{OpenBLAS_HOME}/lib
)

FIND_PATH(OpenBLAS_INCLUDE_DIR NAMES 
    cblas.h PATHS ${Open_BLAS_INCLUDE_SEARCH_PATHS}
)

FIND_LIBRARY(OpenBLAS_LIB NAMES 
    openblas PATHS ${Open_BLAS_LIB_SEARCH_PATHS}
)

message("-- OpenBLAS_INCLUDE_DIR: ${OpenBLAS_INCLUDE_DIR}")
messagE("-- OpenBLAS_LIB: ${OpenBLAS_LIB}")