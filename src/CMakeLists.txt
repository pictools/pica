set(PICA_INCLUDE_DIR ../include)	
set(PICA_HEADER_DIR ${PICA_INCLUDE_DIR}/pica)

include_directories(
    ${PICA_INCLUDE_DIR})
	
add_library(pica
    ${PICA_HEADER_DIR}/grid/Grid.h
    ${PICA_HEADER_DIR}/grid/ScalarField.h	
    ${PICA_HEADER_DIR}/math/Constants.h
    ${PICA_HEADER_DIR}/math/FP.h
    ${PICA_HEADER_DIR}/math/Vectors.h
    ${PICA_HEADER_DIR}/particle/Particle.h
	${PICA_HEADER_DIR}/threading/OpenMPHelper.h
	${PICA_HEADER_DIR}/utility/Utility.h	
    ${PICA_HEADER_DIR}/CurrentDeposition.h	
    ${PICA_HEADER_DIR}/Ensemble.h
    ${PICA_HEADER_DIR}/Fdtd.h
    ${PICA_HEADER_DIR}/FieldSolver.h		
    ${PICA_HEADER_DIR}/Formfactor.h
    ${PICA_HEADER_DIR}/Parameters.h
    ${PICA_HEADER_DIR}/ParticlePush.h	
    ${PICA_HEADER_DIR}/ParticleSystem.h
    ${PICA_HEADER_DIR}/TileInterpolator.h	

    grid/Grid.cpp	
    grid/ScalarField.cpp	
    particle/Particle.cpp	
	utility/Utility.cpp	
	CurrentDeposition.cpp
    Ensemble.cpp
	Fdtd.cpp
	FieldSolver.cpp
	ParticlePush.cpp
    ParticleSystem.cpp)