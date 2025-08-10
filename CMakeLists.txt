cmake_minimum_required (VERSION 3.27.2)
project(crax)

# Add cmake module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Include new dependency management modules
include(FetchThirdParty)
include(FindArrowDynamic)

# Options for dependency management (defaults maintain backwards compatibility)
option(USE_AUTOMATED_DEPENDENCIES "Use automated dependency fetching instead of static folder" OFF)
option(USE_DYNAMIC_ARROW "Use dynamic Arrow linking (via system or pyarrow)" OFF)

# Setup dependencies based on options
if(USE_AUTOMATED_DEPENDENCIES)
    message(STATUS "Using automated dependency management")
    setup_third_party_dependencies()
else()
    message(STATUS "Using static dependency folder 'thirdparty'")
    include_directories(
        common/include/
        common/thirdparty/nlohmann3.11.3/
        common/thirdparty/boost_1_87_0/
        common/thirdparty/eigen/
    )
endif()

# Always include common headers
include_directories(common/include/)
# Setup Arrow dependency
setup_arrow_dependency()

if(MSVC)
	if(NOT USE_DYNAMIC_ARROW)
		# Traditional static Arrow for MSVC (backwards compatibility)
		include_directories(common/thirdparty/arrow20.0.0)
	endif()

	set(CMAKE_CXX_STANDARD "23")
	set(CMAKE_CONFIGURATION_TYPES "Debug;Release")
	set(CMAKE_CXX_FLAGS "/Zi /W4 /WX /MP /EHsc /wd4819 /wd4505 /wd4324 /wd4100 /D_WINDOWS /DNOMINMAX /D_SILENCE_CXX23_DENORM_DEPRECATION_WARNING") 
	set(CMAKE_CXX_FLAGS_RELEASE "/O2 /MD") 
	set(CMAKE_CXX_FLAGS_DEBUG "/Od /MDd /D_DEBUG")
elseif(APPLE)
	set(CMAKE_CXX_STANDARD "23")
	set(CMAKE_OSX_DEPLOYMENT_TARGET "13.3")
    message(STATUS "macOS deployment target: ${CMAKE_OSX_DEPLOYMENT_TARGET}")

    # Add deployment target specific flags
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmacosx-version-min=${CMAKE_OSX_DEPLOYMENT_TARGET}")

	set(CMAKE_CXX_FLAGS_DEBUG "-Og -Wall -Wextra -Wfatal-errors -std=c++2b -ggdb -D_DEBUG -Wno-deprecated -Wno-deprecated-declarations")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -Wall -Wextra -Wfatal-errors -std=c++2b -Wno-deprecated -Wno-deprecated-declarations -Wno-unused-parameter")
else()
	set(CMAKE_CXX_FLAGS_DEBUG "-Og -Wall -Wextra -Werror -Wfatal-errors -std=c++2a -ggdb -D_DEBUG")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -Wall -Wextra -Werror -Wfatal-errors -std=c++2a -fopenmp")
endif()

set(CMAKE_BUILD_TYPE "Release")

set(Bin ${PROJECT_SOURCE_DIR}/bin)

# Build common library first
add_subdirectory(common)
add_subdirectory(rad)
add_subdirectory(mesh)