add_subdirectory(gtest)

include_directories(
    gtest/include
    ../include
	src)

add_executable(tests  
    src/math/FPTests.cpp
    src/math/VectorsTests.cpp	
    src/EnsembleTests.cpp
    src/FormfactorTests.cpp
	src/Main.cpp
	src/OpenMPHelperTests.cpp	
    src/ParticleMoveTests.cpp	
    src/ParticleSystemTests.cpp
    src/ParticleTests.cpp
    src/ScalarFieldTests.cpp
	src/TestingUtility.h
	src/TestingUtility.cpp
	src/TileInterpolationTests.cpp
	src/UtilityTests.cpp)
	
if (APPLE)
    target_link_libraries(tests
        -Wl,-all_load
        gtest
        pica
        -Wl,-noall_load)
elseif (UNIX)
    target_link_libraries(tests
        -Wl,-whole-archive
        gtest
        pica
        -Wl,-no-whole-archive)
else()
    target_link_libraries(tests
        gtest
        pica)
endif()