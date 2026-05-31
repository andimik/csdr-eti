# Standalone csdr-eti tool
# This builds csdr-eti as a command-line executable without requiring csdr runtime

cmake_minimum_required(VERSION 3.0...3.27.7)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

project(csdr-eti-standalone VERSION 0.1.0)

enable_language(CXX C)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 99)

include(GNUInstallDirs)

find_package(Csdr REQUIRED)

# Use pkg-config for FFTW3 (same as main CMakeLists.txt)
include(FindPkgConfig)
pkg_check_modules(FFTW3 REQUIRED fftw3f)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)

# Build the standalone executable - EXCLUDE csdr-eti.cpp to avoid duplicate symbols
add_executable(csdr-eti-standalone
    ../src/csdr-eti-standalone.cpp
    ../src/dab.cpp
    ../src/dab_tables.c
    ../src/fic.cpp
    ../src/misc.cpp
    ../src/depuncture.cpp
    ../src/meta.cpp
    ../src/viterbi.c
)

target_compile_options(csdr-eti-standalone PRIVATE
    -Wall -Wextra -Wno-deprecated-declarations
    -Wno-unused-but-set-variable
    -Wno-type-limits
    -DVERSION=\"${PROJECT_VERSION}-standalone\"
)

target_link_libraries(csdr-eti-standalone
    Csdr::csdr
    ${FFTW3_LIBRARIES}
)

target_include_directories(csdr-eti-standalone PRIVATE
    ${FFTW3_INCLUDE_DIRS}
)

install(TARGETS csdr-eti-standalone
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
