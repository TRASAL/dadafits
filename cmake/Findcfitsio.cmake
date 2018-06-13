include(FindPackageHandleStandardArgs)

find_library(CFITSIO_LIBRARY cfitsio)

find_path(CFITSIO_INCLUDE_DIR fitsio.h)

find_package_handle_standard_args(cfitsio DEFAULT_MSG PSRDADA_LIBRARY CFITSIO_INCLUDE_DIR)

set(CFITSIO_LIBRARIES ${CFITSIO_LIBRARY} )
