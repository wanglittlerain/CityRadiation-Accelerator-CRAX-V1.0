cmake_minimum_required (VERSION 3.27.2)
project(crax)
include_directories(
	common/include/
	common/thirdparty/nlohmann3.11.3/
	common/thirdparty/boost_1_87_0/
	common/thirdparty/eigen/
)
if(MSVC)
	include_directories(common/thirdparty/arrow20.0.0)
	set(CMAKE_CXX_STANDARD "23")
	set(CMAKE_CONFIGURATION_TYPES "Debug;Release")
	set(CMAKE_CXX_FLAGS "/Zi /W4 /WX /MP /EHsc /wd4819 /wd4505 /wd4324 /wd4100 /D_WINDOWS /DNOMINMAX /D_SILENCE_CXX23_DENORM_DEPRECATION_WARNING") 
	set(CMAKE_CXX_FLAGS_RELEASE "/O2 /MD") 
	set(CMAKE_CXX_FLAGS_DEBUG "/Od /MDd /D_DEBUG")
elseif(APPLE)
	set(CMAKE_CXX_STANDARD "23")
	set(CMAKE_CXX_FLAGS_DEBUG "-Og -Wall -Wextra -Wfatal-errors -std=c++2b -ggdb -D_DEBUG -Wno-deprecated -Wno-deprecated-declarations")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -Wall -Wextra -Wfatal-errors -std=c++2b -Wno-deprecated -Wno-deprecated-declarations -Wno-unused-parameter")
else()
	set(CMAKE_CXX_FLAGS_DEBUG "-Og -Wall -Wextra -Werror -Wfatal-errors -std=c++2a -ggdb -D_DEBUG")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -Wall -Wextra -Werror -Wfatal-errors -std=c++2a -fopenmp")
endif(MSVC)

set(CMAKE_BUILD_TYPE "Release")

set(Libs ${PROJECT_SOURCE_DIR}/libs)
set(Bin ${PROJECT_SOURCE_DIR}/bin)

# Build common library first
add_subdirectory(common)
add_subdirectory(rad)
add_subdirectory(mesh)