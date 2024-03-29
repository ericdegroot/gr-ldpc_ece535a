INCLUDE(FindPkgConfig)
PKG_CHECK_MODULES(PC_LDPC_ECE535A ldpc_ece535a)

FIND_PATH(
    LDPC_ECE535A_INCLUDE_DIRS
    NAMES ldpc_ece535a/api.h
    HINTS $ENV{LDPC_ECE535A_DIR}/include
        ${PC_LDPC_ECE535A_INCLUDEDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    LDPC_ECE535A_LIBRARIES
    NAMES gnuradio-ldpc_ece535a
    HINTS $ENV{LDPC_ECE535A_DIR}/lib
        ${PC_LDPC_ECE535A_LIBDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/lib
          ${CMAKE_INSTALL_PREFIX}/lib64
          /usr/local/lib
          /usr/local/lib64
          /usr/lib
          /usr/lib64
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(LDPC_ECE535A DEFAULT_MSG LDPC_ECE535A_LIBRARIES LDPC_ECE535A_INCLUDE_DIRS)
MARK_AS_ADVANCED(LDPC_ECE535A_LIBRARIES LDPC_ECE535A_INCLUDE_DIRS)

