add_subdirectory(gtest)

include_directories(
    gtest/include
    ../include
	src)

add_executable(tests
    src/fieldInterpolation/FormfactorTests.cpp
	src/fieldInterpolation/TileInterpolationTests.cpp
	src/fieldSolver/YeeSolverTests.cpp
    src/grid/ScalarFieldTests.cpp
    src/grid/YeeGridTests.cpp
    src/math/ConstantsTests.cpp
    src/math/DimensionTests.cpp
    src/math/FPTests.cpp
    src/math/VectorsTests.cpp
    src/particles/EnsembleTests.cpp
    src/particles/ParticleTests.cpp
    src/particles/ParticleSystemTests.cpp
    src/particlePush/ParticleMoveTests.cpp
	src/threading/OpenMPHelperTests.cpp
	src/utility/ArrayTests.cpp
	src/utility/UtilityTests.cpp
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