include_directories(
    ../include
	src)

add_executable(benchmark
    src/Output.h
    src/Parameters.h
    src/Parser.h
    src/PerformanceTracker.h
    src/Simulation.h
    src/Timer.h

    src/Main.cpp
    src/Output.cpp
	src/Parser.cpp
    src/PerformanceTracker.cpp
    src/Simulation.cpp	
    src/Timer.cpp)
	
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