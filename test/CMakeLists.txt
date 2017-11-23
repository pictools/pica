add_subdirectory(gtest)

include_directories(
    gtest/include
    ../include
	src)

add_executable(tests
    src/fieldInterpolation/FieldInterpolatorTests.cpp
	src/fieldSolver/YeeSolverTests.cpp
    src/grid/YeeGridTests.cpp
    src/math/ConstantsTests.cpp
    src/math/DimensionTests.cpp
    src/math/FPTests.cpp
    src/math/VectorsTests.cpp
    src/particles/ParticleTests.cpp
    src/particles/ParticleArrayTests.cpp
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
        -Wl,-noall_load)
elseif (UNIX)
    target_link_libraries(tests
        -Wl,-whole-archive
        gtest
        -Wl,-no-whole-archive)
else()
    target_link_libraries(tests
        gtest)
endif()

set_target_properties(gtest PROPERTIES FOLDER Tests)
set_target_properties(tests PROPERTIES FOLDER Tests)