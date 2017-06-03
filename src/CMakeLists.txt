# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

# Generate Release builds by default
if ( NOT CMAKE_BUILD_TYPE )
    set ( CMAKE_BUILD_TYPE "Release" CACHE STRING "Build Type" FORCE )
endif()

# Name of project and that it is C++ only.
project (sim_iCLIP CXX)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

LIST ( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/ )
find_package ( OpenMP REQUIRED )
include ( SeqAn )

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (SYSTEM ${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# Add executable and link against SeqAn dependencies.
add_executable (sim_iclip
                    sim_iclip.cpp
                    sim_iclip.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (sim_iclip ${SEQAN_LIBRARIES})
