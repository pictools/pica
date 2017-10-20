include_directories(
    ../include
	src)

add_executable(benchmark
    src/Parameters.h
    src/Parser.h
	
    src/Main.cpp
	src/Parser.cpp)
	
if (APPLE)
    target_link_libraries(benchmark
        -Wl,-all_load
        pica
        -Wl,-noall_load)
elseif (UNIX)
    target_link_libraries(benchmark
        -Wl,-whole-archive
        pica
        -Wl,-no-whole-archive)
else()
    target_link_libraries(benchmark
        pica)
endif()