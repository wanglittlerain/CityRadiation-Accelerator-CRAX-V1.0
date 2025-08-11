include_directories(./)
file(GLOB_RECURSE CODELIST "${PROJECT_SOURCE_DIR}/rad/*.cc" "${PROJECT_SOURCE_DIR}/rad/*.h")
source_group("" FILES ${CODELIST})

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${Bin})
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${Bin})

add_executable(radiation ${CODELIST})

# Arrow dependency for radiation target
# 1. Try local static Arrow first
if(EXISTS "${CMAKE_SOURCE_DIR}/common/thirdparty/arrow20.0.0")
    include_directories("${CMAKE_SOURCE_DIR}/common/thirdparty/arrow20.0.0")
    set(ARROW_LIBS arrow)
    message(STATUS "Found local Arrow: Linking radiation with local static Arrow")
elseif(USE_DYNAMIC_ARROW)
    # 2. Try dynamic Arrow if requested	
    find_package(Arrow QUIET)
    if(TARGET Arrow::arrow_shared)
        set(ARROW_LIBS Arrow::arrow_shared)
        get_target_property(ARROW_LOCATION Arrow::arrow_shared LOCATION)
        message(STATUS "Linking radiation with Arrow::arrow_shared: ${ARROW_LOCATION}")
    else()
        message(WARNING "Dynamic Arrow requested but not found, linking radiation without Arrow")
        set(ARROW_LIBS "")
    endif()
else()
    # 3. No Arrow
    message(WARNING "No Arrow linked for radiation")
    set(ARROW_LIBS "")
endif()

if(APPLE)
	# macOS doesn't support -static-libgcc and -static-libstdc++
	target_link_libraries(radiation common ${ARROW_LIBS})
else()
	target_link_libraries(radiation common ${ARROW_LIBS} -static-libgcc -static-libstdc++)
endif()
