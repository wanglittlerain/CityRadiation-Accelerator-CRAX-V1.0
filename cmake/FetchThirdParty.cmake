# FetchThirdParty.cmake
# Automated third-party dependency management for CRAX
# This module handles fetching and configuring external dependencies

cmake_minimum_required(VERSION 3.27.2)

include(FetchContent)

# Options for dependency management
option(USE_SYSTEM_DEPENDENCIES "Try to use system-installed dependencies first" OFF)
option(USE_STATIC_DEPENDENCIES "Use static third-party folder as fallback" ON)

# Dependency versions - update these when needed
set(BOOST_VERSION "1.87.0" CACHE STRING "Boost version to fetch")
set(EIGEN_VERSION "3.4.0" CACHE STRING "Eigen version to fetch") 
set(NLOHMANN_JSON_VERSION "3.11.3" CACHE STRING "nlohmann/json version to fetch")

# Function to setup third-party dependencies
function(setup_third_party_dependencies)
    message(STATUS "Setting up third-party dependencies...")
    message(STATUS "  USE_SYSTEM_DEPENDENCIES: ${USE_SYSTEM_DEPENDENCIES}")
    message(STATUS "  USE_STATIC_DEPENDENCIES: ${USE_STATIC_DEPENDENCIES}")

    # Setup Boost
    setup_boost_dependency()
    
    # Setup Eigen
    setup_eigen_dependency()
    
    # Setup nlohmann JSON  
    setup_nlohmann_json_dependency()
    
    message(STATUS "Third-party dependencies setup complete")
endfunction()

# Boost dependency setup
function(setup_boost_dependency)
    set(boost_found FALSE)
    
    # Try system Boost first if requested
    if(USE_SYSTEM_DEPENDENCIES)
        find_package(Boost ${BOOST_VERSION} QUIET)
        if(Boost_FOUND)
            message(STATUS "Using system Boost ${Boost_VERSION}")
            set(boost_found TRUE)
        endif()
    endif()
    
    # Try static folder if enabled and system not found
    if(NOT boost_found AND USE_STATIC_DEPENDENCIES)
        set(BOOST_STATIC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/common/thirdparty/boost_1_87_0")
        if(EXISTS "${BOOST_STATIC_PATH}/boost")
            message(STATUS "Using static Boost from: ${BOOST_STATIC_PATH}")
            # Create interface library for static Boost
            add_library(Boost::headers INTERFACE IMPORTED GLOBAL)
            target_include_directories(Boost::headers INTERFACE "${BOOST_STATIC_PATH}")
            set(boost_found TRUE)
        endif()
    endif()
    
    # Fetch from source as fallback - fetch complete Boost repository
    if(NOT boost_found)
        message(STATUS "Fetching Boost ${BOOST_VERSION} from source...")
        FetchContent_Declare(
            boost
            URL https://github.com/boostorg/boost/releases/download/boost-${BOOST_VERSION}/boost-${BOOST_VERSION}-b2-nodocs.tar.xz
        )
        FetchContent_MakeAvailable(boost)
        include_directories(${boost_SOURCE_DIR})
        message(STATUS "Boost fetched to: ${boost_SOURCE_DIR}")
    endif()
endfunction()

# Eigen dependency setup
function(setup_eigen_dependency)
    set(eigen_found FALSE)
    
    # Try system Eigen first if requested
    if(USE_SYSTEM_DEPENDENCIES)
        find_package(Eigen3 ${EIGEN_VERSION} QUIET)
        if(TARGET Eigen3::Eigen)
            message(STATUS "Using system Eigen ${EIGEN3_VERSION}")
            set(eigen_found TRUE)
        endif()
    endif()
    
    # Try static folder if enabled and system not found
    if(NOT eigen_found AND USE_STATIC_DEPENDENCIES)
        set(EIGEN_STATIC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/common/thirdparty/eigen")
        if(EXISTS "${EIGEN_STATIC_PATH}/Eigen")
            message(STATUS "Using static Eigen from: ${EIGEN_STATIC_PATH}")
            # Create interface library for static Eigen
            add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
            target_include_directories(Eigen3::Eigen INTERFACE "${EIGEN_STATIC_PATH}")
            set(eigen_found TRUE)
        endif()
    endif()
    
    # Fetch from source as fallback
    if(NOT eigen_found)
        message(STATUS "Fetching Eigen ${EIGEN_VERSION} from source...")
        
        # Suppress FetchContent_Populate deprecation warning
        cmake_policy(SET CMP0169 OLD)
        
        FetchContent_Declare(
            eigen
            URL https://gitlab.com/libeigen/eigen/-/archive/${EIGEN_VERSION}/eigen-${EIGEN_VERSION}.tar.gz
        )

        # Use FetchContent_Populate to prevent building of eigen
        FetchContent_GetProperties(eigen)
        if(NOT eigen_POPULATED)
            FetchContent_Populate(eigen)
            add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
            target_include_directories(Eigen3::Eigen INTERFACE "${eigen_SOURCE_DIR}")
        endif()
        
        message(STATUS "Eigen fetched to: ${eigen_SOURCE_DIR}")
    endif()
endfunction()

# nlohmann JSON dependency setup
function(setup_nlohmann_json_dependency)
    set(nlohmann_json_found FALSE)
    
    # Try system nlohmann_json first if requested
    if(USE_SYSTEM_DEPENDENCIES)
        find_package(nlohmann_json ${NLOHMANN_JSON_VERSION} QUIET)
        if(TARGET nlohmann_json::nlohmann_json)
            message(STATUS "Using system nlohmann_json ${nlohmann_json_VERSION}")
            set(nlohmann_json_found TRUE)
        endif()
    endif()
    
    # Try static folder if enabled and system not found
    if(NOT nlohmann_json_found AND USE_STATIC_DEPENDENCIES)
        set(NLOHMANN_STATIC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/common/thirdparty/nlohmann3.11.3")
        if(EXISTS "${NLOHMANN_STATIC_PATH}/nlohmann")
            message(STATUS "Using static nlohmann_json from: ${NLOHMANN_STATIC_PATH}")
            # Create interface library for static nlohmann_json
            add_library(nlohmann_json::nlohmann_json INTERFACE IMPORTED GLOBAL)
            target_include_directories(nlohmann_json::nlohmann_json INTERFACE "${NLOHMANN_STATIC_PATH}")
            set(nlohmann_json_found TRUE)
        endif()
    endif()
    
    # Fetch from source as fallback
    if(NOT nlohmann_json_found)
        message(STATUS "Fetching nlohmann_json ${NLOHMANN_JSON_VERSION} from source...")
        
        FetchContent_Declare(
            nlohmann_json
            URL https://github.com/nlohmann/json/releases/download/v${NLOHMANN_JSON_VERSION}/json.tar.xz 
        )
        
        set(JSON_BuildTests OFF CACHE INTERNAL "")
        FetchContent_MakeAvailable(nlohmann_json)
        
        message(STATUS "nlohmann_json fetched to: ${nlohmann_json_SOURCE_DIR}")
    endif()
endfunction()