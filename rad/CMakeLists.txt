include_directories(./)
file(GLOB_RECURSE CODELIST "${PROJECT_SOURCE_DIR}/rad/*.cc" "${PROJECT_SOURCE_DIR}/rad/*.h")
source_group("" FILES ${CODELIST})

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${Bin})
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${Bin})

add_executable(radiation ${CODELIST})
if(APPLE)
	# macOS doesn't support -static-libgcc and -static-libstdc++
	target_link_libraries(radiation common)
else()
	target_link_libraries(radiation common -static-libgcc -static-libstdc++)
endif()
