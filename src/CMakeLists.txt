set(PICA_INCLUDE_DIR ../include)	
set(PICA_HEADER_DIR ${PICA_INCLUDE_DIR}/pica)

include_directories(
    ${PICA_INCLUDE_DIR})
	
add_library(pica
    ${PICA_HEADER_DIR}/fieldInterpolation/FieldInterpolator.h	
    ${PICA_HEADER_DIR}/fieldInterpolation/Formfactor.h
    ${PICA_HEADER_DIR}/fieldInterpolation/TileInterpolator.h	
    ${PICA_HEADER_DIR}/fieldSolver/YeeSolver.h	
    ${PICA_HEADER_DIR}/grid/Grid.h
    ${PICA_HEADER_DIR}/grid/YeeGrid.h
    ${PICA_HEADER_DIR}/grid/ScalarField.h
    ${PICA_HEADER_DIR}/math/Constants.h
    ${PICA_HEADER_DIR}/math/Dimension.h	
    ${PICA_HEADER_DIR}/math/FP.h
    ${PICA_HEADER_DIR}/math/Vectors.h
    ${PICA_HEADER_DIR}/particles/Ensemble.h
    ${PICA_HEADER_DIR}/particles/EnsembleOrdered.h
    ${PICA_HEADER_DIR}/particles/EnsembleSupercells.h
    ${PICA_HEADER_DIR}/particles/EnsembleUnordered.h
    ${PICA_HEADER_DIR}/particles/Particle.h
    ${PICA_HEADER_DIR}/particles/ParticleArray.h
    ${PICA_HEADER_DIR}/particles/ParticleInterface.h
    ${PICA_HEADER_DIR}/particles/ParticleSystem.h
    ${PICA_HEADER_DIR}/particles/ParticleTraits.h
    ${PICA_HEADER_DIR}/particlePush/BorisPusher.h	
    ${PICA_HEADER_DIR}/particlePush/ParticlePush.h
	${PICA_HEADER_DIR}/threading/OpenMPHelper.h
	${PICA_HEADER_DIR}/utility/Array.h	
	${PICA_HEADER_DIR}/utility/Assert.h	
	${PICA_HEADER_DIR}/utility/Utility.h	
    ${PICA_HEADER_DIR}/Parameters.h

    grid/Grid.cpp	
    grid/ScalarField.cpp	
    particles/Ensemble.cpp
    particles/Particle.cpp
    particles/ParticleSystem.cpp
	particlePush/ParticlePush.cpp	
	utility/Utility.cpp)