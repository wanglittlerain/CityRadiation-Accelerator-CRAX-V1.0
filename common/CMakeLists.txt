cmake_minimum_required (VERSION 3.27.2)
project (common)
include_directories(
	./
	./include
	./thirdparty/nlohmann3.11.3/
	./thirdparty/boost_1_87_0/
	./thirdparty/eigen/
)

set(include "${PROJECT_SOURCE_DIR}/include/")
set(src "${PROJECT_SOURCE_DIR}/src/")
file(GLOB_RECURSE h_list "${include}/*.h")
file(GLOB_RECURSE hpp_list "${include}/*.hpp")
file(GLOB_RECURSE cc_list "${src}/*.cc")
set(CMAKE_EXE_LINKER_FLAGS " -static")
set(Head ${h_list} ${hpp_list})
set(Code ${cc_list} ${c_list})
if(MSVC)
	set(CMAKE_CXX_STANDARD "23")
	set(CMAKE_CONFIGURATION_TYPES "Debug;Release") 
	set(CMAKE_CXX_FLAGS "/Zi /W4 /WX /MP /EHsc /wd4819 /wd4505 /D_WINDOWS /D_SILENCE_CXX23_DENORM_DEPRECATION_WARNING") 
	set(CMAKE_CXX_FLAGS_RELEASE "/O2 /MD") 
	set(CMAKE_CXX_FLAGS_DEBUG "/Od /MDd /D_DEBUG") 
else()
	
	set(CMAKE_CXX_FLAGS_DEBUG "-Og -Wall -Wextra -Werror -Wfatal-errors -std=c++2a -ggdb -D_DEBUG")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -Wall -Wextra -Werror -Wfatal-errors -std=c++2a -fopenmp")
endif(MSVC)

source_group("" FILES ${Code})
source_group("" FILES ${Head})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${CMAKE_SOURCE_DIR}/../libs)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/../libs)
add_library(common STATIC ${Code} ${Head})
