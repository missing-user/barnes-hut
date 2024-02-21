cmake_minimum_required(VERSION 3.13)
file(GLOB HEADER_LIST CONFIGURE_DEPENDS "../include/*.h")
file(GLOB LIB_SRC_FILES CONFIGURE_DEPENDS  "*.cpp")

add_library(barneshut STATIC ${LIB_SRC_FILES} ${HEADER_LIST})
target_compile_features(barneshut PUBLIC cxx_std_17)

target_include_directories(barneshut PUBLIC ../include)

source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})


if(OpenMP_CXX_FOUND)
    target_link_libraries(barneshut PUBLIC OpenMP::OpenMP_CXX)
endif()

target_include_directories(barneshut PUBLIC ../include)

# Boost (for commandline program options)
target_link_libraries(barneshut PRIVATE ${Boost_LIBRARIES})

# GLM
target_include_directories(barneshut PRIVATE ${GLM_INCLUDE_DIRS})
target_link_libraries(barneshut PRIVATE ${GLM_INCLUDE_DIRS})

# XSIMD
target_include_directories(barneshut PRIVATE ${xsimd_INCLUDE_DIRS})

# LIBMORTON
target_include_directories(barneshut PRIVATE libmorton::libmorton)

