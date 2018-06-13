include(FindPackageHandleStandardArgs)

find_library(PSRDADA_LIBRARY psrdada)

# handle the QUIETLY and REQUIRED arguments and set PSRDADA_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(psrdada DEFAULT_MSG PSRDADA_LIBRARY )

mark_as_advanced( PSRDADA_LIBRARY )

set(PSRDADA_LIBRARIES ${PSRDADA_LIBRARY} )