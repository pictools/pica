add_subdirectory(gtest)

include_directories(
    gtest/include
    ../include
	src)

add_executable(tests
    src/fieldInterpolation/FormfactorTests.cpp
	src/fieldInterpolation/TileInterpolationTests.cpp
    src/grid/ScalarFieldTests.cpp
    src/math/FPTests.cpp
    src/math/VectorsTests.cpp
    src/particle/ParticleTests.cpp
	src/threading/OpenMPHelperTests.cpp	
	src/utility/UtilityTests.cpp	
    src/EnsembleTests.cpp
    src/ParticleMoveTests.cpp	
    src/ParticleSystemTests.cpp
	src/Main.cpp
	src/TestingUtility.h
	src/TestingUtility.cpp)
	
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