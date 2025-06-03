include_directories(./)
file(GLOB_RECURSE CODELIST "${PROJECT_SOURCE_DIR}/mesh/*.cc" "${PROJECT_SOURCE_DIR}/mesh/*.h")
source_group("" FILES ${CODELIST})
if(NOT APPLE)
	set(CMAKE_EXE_LINKER_FLAGS " -static")
endif()
link_directories(${Libs})
set(LinkLib common)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${Bin})
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${Bin})

add_executable(mesh-generation ${CODELIST})
if(APPLE)
	# macOS doesn't support -static-libgcc and -static-libstdc++
	target_link_libraries(mesh-generation ${LinkLib})
else()
	target_link_libraries(mesh-generation ${LinkLib} -static-libgcc -static-libstdc++)
endif()
