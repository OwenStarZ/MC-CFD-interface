cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project (h5rw_so)

set ( SRC_LIST
      ./src/h5rw.cpp
)

include_directories ( ./include )
include_directories ( /usr/local/hdf5/1.10.5/include )

set (LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

link_directories (/usr/local/hdf5/1.10.5/lib )

add_library (h5rw SHARED ${SRC_LIST})

target_link_libraries (h5rw hdf5)
