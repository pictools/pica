set(PICA_INCLUDE_DIR ../include)	
set(PICA_HEADER_DIR ${PICA_INCLUDE_DIR}/pica)

include_directories(
    ${PICA_INCLUDE_DIR})
	
add_library(pica
    ${PICA_HEADER_DIR}/currentDeposition/CurrentDeposition.h
    ${PICA_HEADER_DIR}/fieldInterpolation/Formfactor.h
    ${PICA_HEADER_DIR}/fieldInterpolation/TileInterpolator.h	
    ${PICA_HEADER_DIR}/fieldSolver/Fdtd.h
    ${PICA_HEADER_DIR}/fieldSolver/FieldSolver.h
    ${PICA_HEADER_DIR}/grid/Grid.h
    ${PICA_HEADER_DIR}/grid/ScalarField.h	
    ${PICA_HEADER_DIR}/math/Constants.h
    ${PICA_HEADER_DIR}/math/FP.h
    ${PICA_HEADER_DIR}/math/Vectors.h
    ${PICA_HEADER_DIR}/particles/Ensemble.h
    ${PICA_HEADER_DIR}/particles/Particle.h
    ${PICA_HEADER_DIR}/particles/ParticleSystem.h
	${PICA_HEADER_DIR}/threading/OpenMPHelper.h
	${PICA_HEADER_DIR}/utility/Utility.h	
    ${PICA_HEADER_DIR}/Parameters.h
    ${PICA_HEADER_DIR}/ParticlePush.h	

	currentDeposition/CurrentDeposition.cpp
	fieldSolver/Fdtd.cpp
	fieldSolver/FieldSolver.cpp	
    grid/Grid.cpp	
    grid/ScalarField.cpp	
    particles/Ensemble.cpp
    particles/Particle.cpp
    particles/ParticleSystem.cpp	
	utility/Utility.cpp	
	ParticlePush.cpp)