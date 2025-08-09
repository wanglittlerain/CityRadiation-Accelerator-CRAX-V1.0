include_directories(./)
file(GLOB_RECURSE CODELIST "${PROJECT_SOURCE_DIR}/rad/*.cc" "${PROJECT_SOURCE_DIR}/rad/*.h")
source_group("" FILES ${CODELIST})

# Link directory setup
link_directories(${Libs})

# Determine Arrow linking strategy based on configuration
if(USE_DYNAMIC_ARROW)
	message(STATUS "Using dynamic Arrow linking for radiation target")
	if(TARGET Arrow::arrow_shared)
		set(LinkLib PRIVATE Arrow::arrow_shared common)
		message(STATUS "  Linking with Arrow::arrow_shared")
	else()
		message(WARNING "  Dynamic Arrow requested but Arrow::arrow_shared target not available")
		set(LinkLib PRIVATE common)
	endif()
elseif(MSVC)
	# Traditional static Arrow linking for MSVC
	message(STATUS "Using static Arrow linking for MSVC radiation target")
	set(LinkLib PRIVATE arrow common)
else()
	# Traditional Unix behavior - try to find Arrow, fallback gracefully
	message(STATUS "Using traditional Arrow detection for Unix radiation target")
	find_package(Arrow QUIET)
	if(Arrow_FOUND)
		set(LinkLib PRIVATE Arrow::arrow_shared common)
		message(STATUS "  Found system Arrow, linking dynamically")
	else()
		message(WARNING "  Arrow not found via find_package, linking without Arrow")
		set(LinkLib PRIVATE common)
	endif()
endif()

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${Bin})
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${Bin})

add_executable(radiation ${CODELIST})
if(APPLE)
	# macOS doesn't support -static-libgcc and -static-libstdc++
	target_link_libraries(radiation ${LinkLib})
else()
	target_link_libraries(radiation ${LinkLib} -static-libgcc -static-libstdc++)
endif()
