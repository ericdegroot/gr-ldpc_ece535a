# Copyright 2011,2012 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.

########################################################################
# Setup library
########################################################################
include(GrPlatform) #define LIB_SUFFIX

include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})

list(APPEND ldpc_ece535a_sources
    ldpc_encoder_bc_impl.cc
    ldpc_decoder_cb_impl.cc
    image_sink_impl.cc
)

set(ldpc_ece535a_sources "${ldpc_ece535a_sources}" PARENT_SCOPE)
if(NOT ldpc_ece535a_sources)
	MESSAGE(STATUS "No C++ sources... skipping lib/")
	return()
endif(NOT ldpc_ece535a_sources)

add_library(gnuradio-ldpc_ece535a SHARED ${ldpc_ece535a_sources})
target_link_libraries(gnuradio-ldpc_ece535a ${Boost_LIBRARIES} ${GNURADIO_ALL_LIBRARIES} gnuradio-digital lapack lapacke)
set_target_properties(gnuradio-ldpc_ece535a PROPERTIES DEFINE_SYMBOL "gnuradio_ldpc_ece535a_EXPORTS")

if(APPLE)
    set_target_properties(gnuradio-ldpc_ece535a PROPERTIES
        INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib"
    )
endif(APPLE)

########################################################################
# Install built library files
########################################################################
install(TARGETS gnuradio-ldpc_ece535a
    LIBRARY DESTINATION lib${LIB_SUFFIX} # .so/.dylib file
    ARCHIVE DESTINATION lib${LIB_SUFFIX} # .lib file
    RUNTIME DESTINATION bin              # .dll file
)

########################################################################
# Build and register unit test
########################################################################
include(GrTest)

include_directories(${CPPUNIT_INCLUDE_DIRS})

list(APPEND test_ldpc_ece535a_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/test_ldpc_ece535a.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/qa_ldpc_ece535a.cc
)

add_executable(test-ldpc_ece535a ${test_ldpc_ece535a_sources})

target_link_libraries(
  test-ldpc_ece535a
  ${GNURADIO_RUNTIME_LIBRARIES}
  ${Boost_LIBRARIES}
  ${CPPUNIT_LIBRARIES}
  gnuradio-ldpc_ece535a
)

GR_ADD_TEST(test_ldpc_ece535a test-ldpc_ece535a)
