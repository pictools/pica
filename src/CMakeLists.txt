set(PICA_INCLUDE_DIR ../include)	
set(PICA_HEADER_DIR ${PICA_INCLUDE_DIR}/pica)

include_directories(
    ${PICA_INCLUDE_DIR})
	
add_library(pica
    ${PICA_HEADER_DIR}/math/Constants.h
    ${PICA_HEADER_DIR}/math/FP.h
    ${PICA_HEADER_DIR}/math/Vectors.h
    ${PICA_HEADER_DIR}/CurrentDeposition.h	
    ${PICA_HEADER_DIR}/Ensemble.h
    ${PICA_HEADER_DIR}/Fdtd.h
    ${PICA_HEADER_DIR}/FieldSolver.h		
    ${PICA_HEADER_DIR}/Formfactor.h
    ${PICA_HEADER_DIR}/Grid.h
	${PICA_HEADER_DIR}/OpenMPHelper.h
    ${PICA_HEADER_DIR}/Parameters.h
    ${PICA_HEADER_DIR}/Particle.h
    ${PICA_HEADER_DIR}/ParticlePush.h	
    ${PICA_HEADER_DIR}/ParticleSystem.h
    ${PICA_HEADER_DIR}/ScalarField.h
    ${PICA_HEADER_DIR}/TileInterpolator.h	
	${PICA_HEADER_DIR}/Utility.h
	
	CurrentDeposition.cpp
    Ensemble.cpp
	Fdtd.cpp
	FieldSolver.cpp
    Grid.cpp
    Particle.cpp
	ParticlePush.cpp
    ParticleSystem.cpp
    ScalarField.cpp
	Utility.cpp)