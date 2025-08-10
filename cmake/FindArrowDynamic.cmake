# FindArrowDynamic.cmake
# Advanced Arrow detection with multiple source support including pyarrow
# This module finds Arrow libraries with preference for dynamic linking

cmake_minimum_required(VERSION 3.27.2)

# Options for Arrow detection
option(USE_DYNAMIC_ARROW "Use dynamic Arrow linking (via system or pyarrow)" OFF)
option(PREFER_PYARROW "Prefer pyarrow installation over system Arrow" ON)

# Function to setup Arrow dependencies
function(setup_arrow_dependency)
    message(STATUS "Setting up Arrow dependency...")
    message(STATUS "  USE_DYNAMIC_ARROW: ${USE_DYNAMIC_ARROW}")
    message(STATUS "  PREFER_PYARROW: ${PREFER_PYARROW}")

    set(arrow_found FALSE)
    set(arrow_source "none")
    
    if(USE_DYNAMIC_ARROW)
        # Try dynamic linking sources
        if(PREFER_PYARROW)
            detect_pyarrow_installation()
            if(TARGET Arrow::arrow_shared)
                set(arrow_found TRUE)
                set(arrow_source "pyarrow")
                message(STATUS "Arrow found via: ${arrow_source}")
                return()
            endif()
        endif()
        
        # Try system Arrow if pyarrow not found or not preferred
        if(NOT arrow_found)
            detect_system_arrow()
            if(TARGET Arrow::arrow_shared)
                set(arrow_found TRUE)
                set(arrow_source "system")
                message(STATUS "Arrow found via: ${arrow_source}")
                return()
            endif()
        endif()
        
        # If still not found and pyarrow wasn't tried, try it now
        if(NOT arrow_found AND NOT PREFER_PYARROW)
            detect_pyarrow_installation()
            if(TARGET Arrow::arrow_shared)
                set(arrow_found TRUE)
                set(arrow_source "pyarrow")
                message(STATUS "Arrow found via: ${arrow_source}")
                return()
            endif()
        endif()
    endif()
    
    # Fall back to static Arrow for Windows or if dynamic not found/requested
    if(NOT arrow_found)
        detect_static_arrow()
        if(TARGET Arrow::arrow_static OR DEFINED ARROW_STATIC_LIBS)
            set(arrow_found TRUE)
            set(arrow_source "static")
            message(STATUS "Arrow found via: ${arrow_source}")
            return()
        endif()
    endif()
    
    # Try fetching Arrow headers as final fallback
    if(NOT arrow_found)
        fetch_arrow_headers_fallback()
        if(TARGET Arrow::arrow_shared)
            set(arrow_found TRUE)
            set(arrow_source "fetched_headers")
            message(STATUS "Arrow found via: ${arrow_source}")
            return()
        endif()
    endif()
    
    # If no Arrow found, create empty target and warn
    message(WARNING "Arrow not found via any method")
    if(NOT TARGET Arrow::arrow_shared)
        add_library(Arrow::arrow_shared INTERFACE IMPORTED GLOBAL)
        message(STATUS "Created empty Arrow target - build may fail if Arrow is required")
    endif()
endfunction()

# Detect pyarrow installation
function(detect_pyarrow_installation)
    message(STATUS "Looking for pyarrow installation...")
    
    # Allow user to specify Python executable explicitly
    if(DEFINED PYARROW_PYTHON_EXECUTABLE AND EXISTS "${PYARROW_PYTHON_EXECUTABLE}")
        set(PYTHON_EXE "${PYARROW_PYTHON_EXECUTABLE}")
        message(STATUS "  Using user-specified Python: ${PYTHON_EXE}")
    else()
        # Use system Python discovery
        find_package(Python3 COMPONENTS Interpreter QUIET)
        if(Python3_FOUND)
            set(PYTHON_EXE "${Python3_EXECUTABLE}")
            message(STATUS "  Using system Python: ${PYTHON_EXE}")
        else()
            message(STATUS "  Python3 not found - cannot detect pyarrow")
            message(STATUS "  Tip: Set PYARROW_PYTHON_EXECUTABLE to specify Python path explicitly")
            return()
        endif()
    endif()
    
    # Check if pyarrow is installed
    execute_process(
        COMMAND ${PYTHON_EXE} -c "import pyarrow; print(pyarrow.get_library_dirs()[0])"
        OUTPUT_VARIABLE PYARROW_LIB_DIR
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE PYARROW_RESULT
    )
    
    if(PYARROW_RESULT EQUAL 0 AND EXISTS "${PYARROW_LIB_DIR}")
        message(STATUS "  Found pyarrow library directory: ${PYARROW_LIB_DIR}")
        
        # Get include directory
        execute_process(
            COMMAND ${PYTHON_EXE} -c "import pyarrow; print(pyarrow.get_include())"
            OUTPUT_VARIABLE PYARROW_INCLUDE_DIR
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE PYARROW_INCLUDE_RESULT
        )
        
        if(PYARROW_INCLUDE_RESULT EQUAL 0 AND EXISTS "${PYARROW_INCLUDE_DIR}")
            message(STATUS "  Found pyarrow include directory: ${PYARROW_INCLUDE_DIR}")
            
            # Get version
            execute_process(
                COMMAND ${Python3_EXECUTABLE} -c "import pyarrow; print(pyarrow.__version__)"
                OUTPUT_VARIABLE PYARROW_VERSION
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE
            )
            message(STATUS "  pyarrow version: ${PYARROW_VERSION}")
            
            # Look for Arrow libraries in the conda/pip environment
            # pyarrow is in lib/python3.x/site-packages/pyarrow, so env lib is ../../../../lib
            get_filename_component(ENV_LIB_DIR "${PYARROW_LIB_DIR}/../../../../lib" ABSOLUTE)
            if(EXISTS "${ENV_LIB_DIR}")
                message(STATUS "  Searching for Arrow libraries in environment: ${ENV_LIB_DIR}")
                find_arrow_libraries_in_environment("${ENV_LIB_DIR}")
                
                if(ARROW_CORE_LIB)
                    # Create imported target with proper Arrow library
                    add_library(Arrow::arrow_shared SHARED IMPORTED GLOBAL)
                    set_target_properties(Arrow::arrow_shared PROPERTIES
                        IMPORTED_LOCATION "${ARROW_CORE_LIB}"
                        INTERFACE_INCLUDE_DIRECTORIES "${PYARROW_INCLUDE_DIR}"
                        INTERFACE_COMPILE_DEFINITIONS "ARROW_STATIC=0"
                    )
                    
                    message(STATUS "  Created Arrow::arrow_shared target with environment libraries")
                    # Set parent scope variables immediately 
                    set(ARROW_FOUND TRUE PARENT_SCOPE)
                    set(ARROW_VERSION ${PYARROW_VERSION} PARENT_SCOPE) 
                    return()
                endif()
            endif()
            
            # Fallback: headers-only approach if libraries not found
            message(STATUS "  Arrow libraries not found in environment, using headers-only approach")
            add_library(Arrow::arrow_shared INTERFACE IMPORTED GLOBAL)
            set_target_properties(Arrow::arrow_shared PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${PYARROW_INCLUDE_DIR}"
                INTERFACE_COMPILE_DEFINITIONS "ARROW_STATIC=0;ARROW_PYTHON_RUNTIME_LINK=1"
            )
            
            # Set parent scope variables 
            set(ARROW_FOUND TRUE PARENT_SCOPE)
            set(ARROW_VERSION ${PYARROW_VERSION} PARENT_SCOPE) 
            return()
        endif()
    else()
        message(STATUS "  pyarrow not found or not installed")
    endif()
endfunction()

# Helper function to find main Arrow libraries in environment
function(find_arrow_libraries_in_environment env_lib_dir)
    # Look for the main Arrow library (libarrow.dylib, libarrow.so, etc.)
    set(arrow_lib_names
        "libarrow.dylib"
        "libarrow.so" 
        "libarrow.dll"
        "arrow.dll"
    )
    
    # Also try versioned names
    file(GLOB versioned_libs "${env_lib_dir}/libarrow.*.dylib" "${env_lib_dir}/libarrow.*.so")
    
    # Check unversioned names first
    foreach(lib_name ${arrow_lib_names})
        set(lib_path "${env_lib_dir}/${lib_name}")
        if(EXISTS "${lib_path}")
            set(ARROW_CORE_LIB ${lib_path} PARENT_SCOPE)
            message(STATUS "  Found main Arrow library: ${lib_path}")
            return()
        endif()
    endforeach()
    
    # Check versioned names
    if(versioned_libs)
        list(GET versioned_libs 0 first_versioned)
        set(ARROW_CORE_LIB ${first_versioned} PARENT_SCOPE)
        message(STATUS "  Found versioned Arrow library: ${first_versioned}")
        return()
    endif()
    
    message(STATUS "  No main Arrow library found in ${env_lib_dir}")
endfunction()

# Helper function to find library files in pyarrow installation (legacy)
function(find_library_files_in_pyarrow lib_dir)
    # This function is kept for compatibility but not used in current implementation
    message(STATUS "  Checking pyarrow internal libraries (not used for linking)")
endfunction()

# Detect system Arrow installation
function(detect_system_arrow)
    message(STATUS "Looking for system Arrow installation...")
    
    # Try CMake's find_package first
    find_package(Arrow QUIET)
    if(Arrow_FOUND)
        message(STATUS "  Found system Arrow via find_package")
        set(ARROW_FOUND TRUE PARENT_SCOPE)
        return()
    endif()
    
    # Try pkg-config as fallback
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(ARROW_PC QUIET arrow)
        if(ARROW_PC_FOUND)
            message(STATUS "  Found system Arrow via pkg-config")
            
            # Create imported target from pkg-config info
            add_library(Arrow::arrow_shared SHARED IMPORTED GLOBAL)
            set_target_properties(Arrow::arrow_shared PROPERTIES
                IMPORTED_LOCATION "${ARROW_PC_LIBRARIES}"
                INTERFACE_INCLUDE_DIRECTORIES "${ARROW_PC_INCLUDE_DIRS}"
                INTERFACE_COMPILE_OPTIONS "${ARROW_PC_CFLAGS_OTHER}"
            )
            
            set(ARROW_FOUND TRUE PARENT_SCOPE)
            set(ARROW_VERSION ${ARROW_PC_VERSION} PARENT_SCOPE)
            return()
        endif()
    endif()
    
    message(STATUS "  System Arrow not found")
endfunction()

# Detect static Arrow installation (Windows fallback)
function(detect_static_arrow)
    message(STATUS "Looking for static Arrow installation...")
    
    if(MSVC)
        set(ARROW_STATIC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/common/thirdparty/arrow20.0.0")
        if(EXISTS "${ARROW_STATIC_PATH}")
            message(STATUS "  Found static Arrow at: ${ARROW_STATIC_PATH}")
            
            # Set variables for MSVC static linking (matching existing pattern)
            set(ARROW_STATIC_LIBS "arrow" PARENT_SCOPE)
            set(ARROW_STATIC_INCLUDE_DIR "${ARROW_STATIC_PATH}" PARENT_SCOPE)
            set(ARROW_FOUND TRUE PARENT_SCOPE)
            set(ARROW_SOURCE "static" PARENT_SCOPE)
            return()
        endif()
    endif()
    
    message(STATUS "  Static Arrow not found or not applicable for this platform")
endfunction()

# Create appropriate Arrow targets based on source
function(create_arrow_targets source)
    if(source STREQUAL "static" AND MSVC)
        message(STATUS "Using static Arrow linking for MSVC")
        # For MSVC static linking, we rely on link_directories and library names
        # This matches the existing pattern in rad/CMakeLists.txt
    else()
        # Dynamic linking targets should already be created by detection functions
        if(TARGET Arrow::arrow_shared)
            message(STATUS "Arrow::arrow_shared target ready for use")
        else()
            message(WARNING "Expected Arrow::arrow_shared target not found")
        endif()
    endif()
endfunction()

# Fetch Arrow headers as final fallback when all other methods fail
function(fetch_arrow_headers_fallback)    
    # Include FetchContent if not already available
    include(FetchContent)
    
    # Arrow version to fetch (adjust as needed)
    set(ARROW_FETCH_VERSION "21.0.0" CACHE STRING "Arrow version to fetch for headers")
    
    # Try to fetch Arrow source for headers
    message(STATUS "  Fetching Arrow ${ARROW_FETCH_VERSION} headers from GitHub...")
    
    # Use ExternalProject for better isolation
    include(ExternalProject)
    
    set(ARROW_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/arrow-install")
    
    # Create install directories to avoid path errors
    file(MAKE_DIRECTORY "${ARROW_INSTALL_PREFIX}")
    file(MAKE_DIRECTORY "${ARROW_INSTALL_PREFIX}/include")
    file(MAKE_DIRECTORY "${ARROW_INSTALL_PREFIX}/lib")
    file(MAKE_DIRECTORY "${ARROW_INSTALL_PREFIX}/bin")
    
    # Build Arrow as external project with complete isolation
    message(STATUS "  Building Arrow as isolated external project...")
    ExternalProject_Add(
        arrow_external
        URL "https://github.com/apache/arrow/releases/download/apache-arrow-${ARROW_FETCH_VERSION}/apache-arrow-${ARROW_FETCH_VERSION}.tar.gz"
        SOURCE_SUBDIR cpp
        CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX=${ARROW_INSTALL_PREFIX}
            -DCMAKE_BUILD_TYPE=Release
            -DARROW_BUILD_STATIC=OFF
            -DARROW_BUILD_SHARED=ON
            -DARROW_BUILD_TESTS=OFF
            -DARROW_BUILD_EXAMPLES=OFF
            -DARROW_BUILD_BENCHMARKS=OFF
            -DARROW_WITH_BOOST=OFF
            -DARROW_DEPENDENCY_SOURCE=BUNDLED
            -DARROW_VERBOSE_THIRDPARTY_BUILD=OFF
        BUILD_ALWAYS FALSE
        INSTALL_DIR ${ARROW_INSTALL_PREFIX}
    )
    
    # Create imported target with pre-created directories
    add_library(Arrow::arrow_shared SHARED IMPORTED GLOBAL)
    
    # Set platform-specific library properties
    if(WIN32)
        # Windows needs both DLL location and import library
        set_target_properties(Arrow::arrow_shared PROPERTIES
            IMPORTED_LOCATION "${ARROW_INSTALL_PREFIX}/bin/arrow.dll"
            IMPORTED_IMPLIB "${ARROW_INSTALL_PREFIX}/lib/arrow.lib"
            INTERFACE_INCLUDE_DIRECTORIES "${ARROW_INSTALL_PREFIX}/include"
            INTERFACE_COMPILE_DEFINITIONS "ARROW_STATIC=0"
        )
    else()
        # Unix-like systems
        set_target_properties(Arrow::arrow_shared PROPERTIES
            IMPORTED_LOCATION "${ARROW_INSTALL_PREFIX}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}arrow${CMAKE_SHARED_LIBRARY_SUFFIX}"
            INTERFACE_INCLUDE_DIRECTORIES "${ARROW_INSTALL_PREFIX}/include"
            INTERFACE_COMPILE_DEFINITIONS "ARROW_STATIC=0"
        )
    endif()
    
    # Make sure Arrow is built before any target that uses it
    add_dependencies(Arrow::arrow_shared arrow_external)
    
    message(STATUS "  Arrow external build configured")
    set(ARROW_FOUND TRUE PARENT_SCOPE)
    set(ARROW_VERSION ${ARROW_FETCH_VERSION} PARENT_SCOPE)
    return()
    
    message(STATUS "  Arrow headers fetch fallback failed")
endfunction()