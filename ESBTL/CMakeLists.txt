# cmake file for ESBTL benchmarks
# ESBTL_DIR needs to be set as the ESBTL root directory
project( benchmarks )
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)
include_directories($ESBTL_DIR/include/)
add_executable(parse_pdb parse_pdb.cc)
