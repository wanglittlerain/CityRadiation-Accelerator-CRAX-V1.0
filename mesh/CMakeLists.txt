include_directories(./)
file(GLOB_RECURSE CODELIST "${PROJECT_SOURCE_DIR}/mesh/*.cc" "${PROJECT_SOURCE_DIR}/mesh/*.h")
source_group("" FILES ${CODELIST})
set(CMAKE_EXE_LINKER_FLAGS " -static")
link_directories(${Libs})
set(LinkLib common)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${Bin})
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${Bin})

add_executable(mesh-generation ${CODELIST})
target_link_libraries(mesh-generation ${LinkLib} -static-libgcc -static-libstdc++)
