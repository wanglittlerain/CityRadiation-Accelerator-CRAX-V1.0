cmake_minimum_required (VERSION 3.27.2)
project (common)

# Always include common headers
include_directories(
	./
	./include
)

# Additional include directories for backwards compatibility
if(NOT USE_AUTOMATED_DEPENDENCIES)
	include_directories(
		./thirdparty/nlohmann3.11.3/
		./thirdparty/boost_1_87_0/
		./thirdparty/eigen/
	)
endif()

set(include "${PROJECT_SOURCE_DIR}/include/")
set(src "${PROJECT_SOURCE_DIR}/src/")
file(GLOB_RECURSE h_list "${include}/*.h")
file(GLOB_RECURSE hpp_list "${include}/*.hpp")
file(GLOB_RECURSE cc_list "${src}/*.cc")
if(NOT APPLE)
	set(CMAKE_EXE_LINKER_FLAGS " -static")
endif()
set(Head ${h_list} ${hpp_list})
set(Code ${cc_list} ${c_list})
if(MSVC)
	set(CMAKE_CXX_STANDARD "23")
	set(CMAKE_CONFIGURATION_TYPES "Debug;Release") 
	set(CMAKE_CXX_FLAGS "/Zi /W4 /WX /MP /EHsc /wd4819 /wd4505 /D_WINDOWS /D_SILENCE_CXX23_DENORM_DEPRECATION_WARNING") 
	set(CMAKE_CXX_FLAGS_RELEASE "/O2 /MD") 
	set(CMAKE_CXX_FLAGS_DEBUG "/Od /MDd /D_DEBUG") 
elseif(APPLE)
	set(CMAKE_CXX_STANDARD "23")
	set(CMAKE_CXX_FLAGS_DEBUG "-Og -Wall -Wextra -Wfatal-errors -std=c++2b -ggdb -D_DEBUG -Wno-deprecated -Wno-deprecated-declarations")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -Wall -Wextra -Wfatal-errors -std=c++2b -Wno-deprecated -Wno-deprecated-declarations")
else()
	set(CMAKE_CXX_FLAGS_DEBUG "-Og -Wall -Wextra -Werror -Wfatal-errors -std=c++2a -ggdb -D_DEBUG")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -Wall -Wextra -Werror -Wfatal-errors -std=c++2a -fopenmp")
endif(MSVC)

source_group("" FILES ${Code})
source_group("" FILES ${Head})
add_library(common STATIC ${Code} ${Head})

# Link with modern dependency targets when using automated dependencies
if(USE_AUTOMATED_DEPENDENCIES)
	# Link common library with third-party targets
	if(TARGET Boost::headers)
		target_link_libraries(common PUBLIC Boost::headers)
	endif()
	if(TARGET Eigen3::Eigen)
		target_link_libraries(common PUBLIC Eigen3::Eigen)
	endif()
	if(TARGET nlohmann_json::nlohmann_json)
		target_link_libraries(common PUBLIC nlohmann_json::nlohmann_json)
	endif()
endif()
